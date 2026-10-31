#!/usr/bin/env python3
"""Compact exact prime-state catalogues and full-rho tiled searches."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, exp, log
from pathlib import Path
import random
from time import monotonic

import numpy as np
from scipy.optimize import brentq, differential_evolution
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure
from search_symmetric_tiled_filters import polynomial_hash
from symmetric_prime_filter_synthesis import local_endpoint_image, synthesize


def source_length(prime, modes):
    value = 1
    while value < modes:
        value *= prime
    return value


def degree_vectors(prime, modes, m, seed):
    size = modes * (prime - 1) + 1
    levels = sorted({0, 1, max(1, m // 2), m - 1, m, m + 1, 2 * m, 4 * m})
    generator = random.Random(seed)
    if size <= 5:
        supports = set(range(1 << size))
    else:
        supports = {0, (1 << size) - 1}
        for weight in range(size):
            supports.add(1 << weight)
            supports.add(sum(1 << w for w in range(weight, size)))
            supports.add(sum(1 << w for w in range(weight + 1)))
            supports.add(sum(1 << w for w in range(size) if comb(w, weight) % prime))
        supports.update(generator.randrange(1 << size) for _ in range(16))
    vectors = {tuple(level if (support >> weight) & 1 else 0 for weight in range(size))
               for support in supports for level in levels}
    vectors.update(tuple(generator.choice(levels) for _ in range(size)) for _ in range(16))
    return tuple(sorted(vectors)), levels, tuple(sorted(supports))


def build_catalogue(prime, modes, ramifications, seed):
    last = modes * (prime - 1)
    n0 = source_length(prime, modes)
    catalogue, sources = {}, []
    priorities = (tuple(range(last + 1)), tuple(range(last, -1, -1)),
                  tuple(sorted(range(last + 1), key=lambda w: (abs(w - last // 2), -w))))
    for P in ramifications:
        m = n0 * P
        vectors, levels, supports = degree_vectors(prime, modes, m, seed)
        calls = feasible = 0
        started = monotonic()
        for number, caps in enumerate(vectors):
            maxima = tuple(max(caps[h] + m * (weight - h) for h in range(weight + 1)
                                if comb(weight, h) % prime) + 1 for weight in range(last + 1))
            seen = set()
            for priority in priorities:
                target = [0] * (last + 1)
                for weight in priority:
                    low, high = target[weight], maxima[weight]
                    while high - low > 1:
                        target[weight] = (low + high) // 2
                        calls += 1
                        if synthesize(prime, modes, m, caps, target, construct=False):
                            low = target[weight]
                        else:
                            high = target[weight]
                    target[weight] = low
                if tuple(target) in seen:
                    continue
                seen.add(tuple(target))
                witness = synthesize(prime, modes, m, caps, target)
                assert witness and witness['actual_orders'] == tuple(target)
                feasible += 1
                key = (tuple(Fraction(value, m) for value in witness['actual_orders']),
                       tuple(Fraction(value, m) for value in witness['actual_degrees']))
                if key not in catalogue:
                    catalogue[key] = {
                        'source_interval_length': n0, 'source_m': m, 'source_ramification': P,
                        'source_degree_bounds': caps, 'source_order_bounds': tuple(target),
                        'beta_orders': witness['actual_orders'], 'degrees': witness['actual_degrees'],
                        'endpoint_values': witness['endpoint_values'],
                        'coefficient_field_order': witness['coefficient_field_order'],
                        'witness_polynomial_sha256': polynomial_hash(witness),
                    }
            if (number + 1) % 100 == 0:
                print('CATALOGUE_PROGRESS', prime, modes, P, number + 1, len(vectors),
                      len(catalogue), calls, monotonic() - started, flush=True)
        sources.append({'ramification': P, 'm': m, 'source_interval_length': n0,
                        'degree_vectors': len(vectors), 'degree_levels': levels,
                        'weight_support_masks': supports, 'row_priorities': priorities,
                        'exact_feasibility_queries': calls, 'feasible_target_pairs': feasible,
                        'elapsed_seconds': monotonic() - started})
        local_endpoint_image.cache_clear()
        print('CATALOGUE', prime, modes, P, len(catalogue), calls, monotonic() - started, flush=True)
    return tuple(catalogue.values()), sources


def partitions(n, modes, seed):
    count, remainder = divmod(n, modes)
    indices = list(range(n))
    interleaved = [i + j * count for i in range(count) for j in range(modes)]
    interleaved += list(range(count * modes, n))
    reflected = [value for i in range(n // 2) for value in (i, n - 1 - i)]
    if n % 2:
        reflected.append(n // 2)
    shuffled = indices[:]
    random.Random(seed).shuffle(shuffled)
    results = {}
    signatures = set()
    for name, order in (('contiguous', indices), ('interleaved', interleaved),
                        ('reflected', reflected), ('shuffled', shuffled)):
        blocks = tuple(tuple(sorted(2 * index + 1 for index in order[start:start + modes]))
                       for start in range(0, count * modes, modes))
        spectators = tuple(2 * index + 1 for index in order[count * modes:])
        assert len(spectators) == remainder
        assert sorted([value for block in blocks for value in block] + list(spectators)) == list(range(1, 2 * n, 2))
        signature = tuple(sorted(blocks)), tuple(sorted(spectators))
        if signature not in signatures:
            results[name] = blocks, spectators
            signatures.add(signature)
    return results


class TiledPressure:
    def __init__(self, prime, n, modes, catalogue, blocks):
        self.prime, self.n, self.modes, self.catalogue = prime, n, modes, catalogue
        self.blocks = blocks
        self.energies = np.asarray(blocks, dtype=float) / n
        self.weights = np.arange(modes * (prime - 1) + 1)
        self.orders = np.asarray([[value / row['source_m'] for value in row['beta_orders']] for row in catalogue])
        self.degrees = np.asarray([[value / row['source_m'] for value in row['degrees']] for row in catalogue])
        self.calls = 0

    def relative_cost(self, tau, charge):
        log_coefficients = np.zeros((len(self.blocks), 1))
        for index in range(self.modes):
            width = log_coefficients.shape[1] + self.prime - 1
            expanded = []
            for digit in range(self.prime):
                row = np.full((len(self.blocks), width), -np.inf)
                row[:, digit:digit + log_coefficients.shape[1]] = (
                    log_coefficients - tau * self.energies[:, index, None] * digit)
                expanded.append(row)
            log_coefficients = logsumexp(np.asarray(expanded), axis=0)
        physical = log_coefficients - logsumexp(log_coefficients, axis=1)[:, None]
        shifted = log_coefficients - tau * charge * self.weights
        shifted -= logsumexp(shifted, axis=1)[:, None]
        minimum = np.zeros(len(self.blocks))
        selected = np.full(len(self.blocks), -1, dtype=int)
        for start in range(0, len(self.catalogue), 256):
            stop = min(start + 256, len(self.catalogue))
            values = logsumexp(shifted[:, None, :] +
                2 * tau * (self.weights - self.orders[start:stop])[None, :, :], axis=2)
            values += logsumexp(physical[:, None, :] +
                2 * tau * self.degrees[start:stop][None, :, :], axis=2)
            indices = np.argmin(values, axis=1)
            best = values[np.arange(len(indices)), indices]
            improved = best < minimum - 1e-12
            minimum[improved] = best[improved]
            selected[improved] = start + indices[improved]
        return float(minimum.sum()), selected, minimum

    def root(self, parameters):
        tau, charge = exp(parameters[0]), parameters[1]
        relative, _, _ = self.relative_cost(tau, charge)
        self.calls += 1
        function = lambda rho: gaussian_pressure(self.prime, self.n, self.n * charge, tau / self.n, rho) + relative
        upper = 2.0
        while function(upper) > 0 and upper < 64:
            upper *= 2
        if function(0.4) <= 0:
            raise ArithmeticError('Unexpected strong gain requires a separate rigorous replay.')
        return brentq(function, 0.4, upper, xtol=2e-10)

    def optimize(self, seed, iterations):
        best = (float('inf'), None)
        for tau in (.15, .4, .8, 1.5, 3, 5, 7, 10, 15):
            for charge in (-4, -2, -1, 0, 1, 2, 4):
                point = (log(tau), charge)
                value = self.root(point)
                if value < best[0]:
                    best = value, point
        result = differential_evolution(self.root, ((-4, log(15)), (-6, 6)),
            seed=seed, popsize=8, maxiter=iterations, tol=2e-7, polish=True, x0=best[1])
        tau, charge = exp(result.x[0]), float(result.x[1])
        relative, selected, values = self.relative_cost(tau, charge)
        pressure = gaussian_pressure(self.prime, self.n, self.n * charge, tau / self.n) + relative
        return {'rho': float(result.fun), 'tau': tau, 'scaled_charge': charge,
                'scale': tau / self.n, 'charge': self.n * charge,
                'relative_log_cost': relative, 'full_critical_pressure': pressure,
                'selected_catalogue_indices': selected.tolist(),
                'block_relative_log_costs': values.tolist(),
                'actual_mode_energy_blocks': self.blocks,
                'nonbaseline_blocks': int((selected >= 0).sum()),
                'root_evaluations': self.calls, 'optimizer_success': bool(result.success)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=3)
    parser.add_argument('--block-size', type=int, choices=(2, 4), default=2)
    parser.add_argument('--ramifications', default='1')
    parser.add_argument('--lengths', default='27,81')
    parser.add_argument('--iterations', type=int, default=35)
    parser.add_argument('--seed', type=int, default=720493)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ramifications = tuple(map(int, args.ramifications.split(',')))
    if any(P <= 0 or P % 2 == 0 for P in ramifications):
        raise ValueError('Use positive odd source ramifications.')
    catalogue, sources = build_catalogue(args.prime, args.block_size, ramifications, args.seed)
    payload = json.dumps(catalogue, separators=(',', ':'))
    result = {'status': 'Exact catalogue complete; listed continuous optimizations are exploratory.',
              'prime': args.prime, 'block_size': args.block_size,
              'source_ramifications': ramifications, 'source_runs': sources,
              'catalogue_profiles': len(catalogue), 'catalogue': catalogue,
              'catalogue_sha256': sha256(payload.encode()).hexdigest(),
              'seed': args.seed, 'optimizer_iterations': args.iterations,
              'optimized_runs': [], 'strict_crossing_candidate_found': False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    for n in map(int, args.lengths.split(',')):
        remaining = n
        while remaining > 1 and remaining % args.prime == 0:
            remaining //= args.prime
        if remaining != 1 or n < source_length(args.prime, args.block_size):
            raise ValueError('Use a prime-power physical interval large enough for each block.')
        for name, (blocks, spectators) in partitions(n, args.block_size, args.seed).items():
            started = monotonic()
            run = TiledPressure(args.prime, n, args.block_size, catalogue, blocks)
            row = run.optimize(args.seed, args.iterations)
            row.update(n=n, partition=name, untouched_spectator_mode_energies=spectators,
                       elapsed_seconds=monotonic() - started)
            result['optimized_runs'].append(row)
            result['strict_crossing_candidate_found'] |= row['full_critical_pressure'] < -1e-9
            args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
            print('TILED', args.prime, n, args.block_size, name, row['rho'],
                  row['nonbaseline_blocks'], row['relative_log_cost'], row['elapsed_seconds'], flush=True)
    result['status'] = 'Exact catalogue and every requested genuine-rho optimization completed; no all-family exclusion is inferred.'
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('DONE', args.prime, args.block_size, len(catalogue),
          min(row['rho'] for row in result['optimized_runs']), flush=True)


if __name__ == '__main__':
    main()
