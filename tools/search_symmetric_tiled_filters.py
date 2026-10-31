#!/usr/bin/env python3
"""Build exact symmetric polynomial catalogues and optimize genuine tiled rho.

Feasibility and polynomial profiles are exact. Continuous optimization is
exploratory; no negative floating pressure is a certified crossing witness.
Every retained tile records its actual distinct physical mode energies.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import exp, log, sqrt
from pathlib import Path
import random
from time import monotonic

import numpy as np
from scipy.optimize import brentq, differential_evolution
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure
from symmetric_binary_filter_synthesis import local_endpoint_image, submasks, synthesize


def polynomial_hash(witness):
    payload = [sorted(value.items()) for value in witness['f']]
    return sha256(json.dumps(payload, separators=(',', ':')).encode()).hexdigest()


def degree_catalogue(modes, m, random_cases, seed):
    size = modes + 1
    # Refinement of the fractional degree offsets is through literal odd
    # ramification; no fractional polynomial exponent enters a witness.
    levels = sorted({0, 1, max(1, m // 4), max(1, m // 2), max(1, 3 * m // 4),
                     m - 1, m, m + 1, 2 * m, 2 * m + 1, 4 * m})
    supports = set()
    if modes <= 8:
        supports.update(range(1 << size))
    else:
        for weight in range(size):
            supports.add(1 << weight)
            supports.add(sum(1 << w for w in range(weight, size)))
            supports.add(sum(1 << w for w in range(size) if w & weight == weight))
    result = {tuple(level if (support >> w) & 1 else 0 for w in range(size))
              for support in supports for level in levels}
    # A small fixed multilevel family supplements the two-level supports.
    generator = random.Random(seed)
    for _ in range(random_cases):
        result.add(tuple(generator.choice(levels) for _ in range(size)))
    for threshold in range(size):
        for slope in (1, 2):
            result.add(tuple(max(0, (weight - threshold) * m // 2) * slope
                             for weight in range(size)))
    return tuple(sorted(result)), levels


def build_catalogue(modes, ramifications, random_cases, seed):
    catalogue = {}
    source_runs = []
    priorities = [tuple(range(modes + 1)), tuple(range(modes, -1, -1)),
                  tuple(sorted(range(modes + 1), key=lambda w: (abs(w - modes // 2), -w)))]
    for P in ramifications:
        m = modes * P
        degrees, levels = degree_catalogue(modes, m, random_cases, seed)
        calls = feasible = 0
        started = monotonic()
        for caps in degrees:
            maxima = tuple(max(caps[h] + m * (weight - h) for h in submasks(weight)) + 1
                           for weight in range(modes + 1))
            targets_seen = set()
            for priority in priorities:
                target = [0] * (modes + 1)
                for weight in priority:
                    lower, upper = target[weight], maxima[weight]
                    while upper - lower > 1:
                        target[weight] = (lower + upper) // 2
                        calls += 1
                        if synthesize(modes, m, caps, target, construct=False):
                            lower = target[weight]
                        else:
                            upper = target[weight]
                    target[weight] = lower
                if tuple(target) in targets_seen:
                    continue
                targets_seen.add(tuple(target))
                witness = synthesize(modes, m, caps, target)
                assert witness and witness['actual_orders'] == tuple(target)
                feasible += 1
                key = (tuple(Fraction(value, m) for value in witness['actual_orders']),
                       tuple(Fraction(value, m) for value in witness['actual_degrees']))
                if key not in catalogue:
                    catalogue[key] = {
                        'source_m': m, 'source_ramification': P,
                        'source_degree_bounds': caps, 'source_order_bounds': tuple(target),
                        'beta_orders': witness['actual_orders'], 'degrees': witness['actual_degrees'],
                        'endpoint_values': witness['endpoint_values'],
                        'coefficient_field_order': witness['coefficient_field_order'],
                        'witness_polynomial_sha256': polynomial_hash(witness),
                    }
        source_runs.append({'ramification': P, 'm': m, 'degree_vectors': len(degrees),
                            'degree_levels': levels, 'row_priorities': priorities,
                            'exact_feasibility_queries': calls, 'feasible_target_pairs': feasible,
                            'elapsed_seconds': monotonic() - started})
        print('CATALOGUE', modes, P, len(catalogue), calls, monotonic() - started, flush=True)
        local_endpoint_image.cache_clear()
    return tuple(catalogue.values()), source_runs


def partitions(n, modes, seed):
    if n % modes:
        raise ValueError('This tiled driver requires the block size to divide n.')
    indices = list(range(n))
    orders = {'contiguous': indices,
              'interleaved': [i + j * (n // modes) for i in range(n // modes) for j in range(modes)],
              'reflected': [value for i in range(n // 2) for value in (i, n - 1 - i)]}
    generator = random.Random(seed)
    shuffled = indices[:]
    generator.shuffle(shuffled)
    orders['shuffled'] = shuffled
    result = {}
    signatures = set()
    for name, order in orders.items():
        blocks = tuple(tuple(sorted(2 * index + 1 for index in order[start:start + modes]))
                       for start in range(0, n, modes))
        signature = tuple(sorted(blocks))
        assert sorted(value for block in blocks for value in block) == list(range(1, 2 * n, 2))
        if signature not in signatures:
            signatures.add(signature)
            result[name] = blocks
    return result


class TiledPressure:
    def __init__(self, n, modes, catalogue, blocks):
        self.n, self.modes, self.catalogue = n, modes, catalogue
        self.blocks = blocks
        self.energies = np.asarray(blocks, dtype=float) / n
        self.weights = np.arange(modes + 1)
        self.orders = np.asarray([[value / row['source_m'] for value in row['beta_orders']]
                                  for row in catalogue])
        self.degrees = np.asarray([[value / row['source_m'] for value in row['degrees']]
                                   for row in catalogue])
        self.calls = 0

    def relative_cost(self, tau, scaled_charge):
        log_coefficients = np.zeros((len(self.blocks), 1))
        for index in range(self.modes):
            padding = np.full((len(self.blocks), 1), -np.inf)
            log_coefficients = np.logaddexp(
                np.c_[log_coefficients, padding],
                np.c_[padding, log_coefficients - tau * self.energies[:, index, None]])
        physical = log_coefficients - logsumexp(log_coefficients, axis=1)[:, None]
        shifted = log_coefficients - tau * scaled_charge * self.weights
        shifted -= logsumexp(shifted, axis=1)[:, None]
        # The baseline is always admissible and has relative cost exactly 0.
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
        tau, scaled_charge = exp(parameters[0]), parameters[1]
        relative, _, _ = self.relative_cost(tau, scaled_charge)
        self.calls += 1
        function = lambda rho: gaussian_pressure(
            2, self.n, self.n * scaled_charge, tau / self.n, rho) + relative
        upper = 2.0
        while function(upper) > 0 and upper < 64:
            upper *= 2
        if function(0.4) <= 0:
            raise ArithmeticError('Unexpected strong gain: isolate and rigorously replay before proceeding.')
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
        relative, selected, pieces = self.relative_cost(tau, charge)
        pressure = gaussian_pressure(2, self.n, self.n * charge, tau / self.n) + relative
        return {
            'rho': float(result.fun), 'tau': tau, 'scaled_charge': charge,
            'scale': tau / self.n, 'charge': self.n * charge,
            'relative_log_cost': relative, 'full_critical_pressure': pressure,
            'selected_catalogue_indices': selected.tolist(),
            'block_relative_log_costs': pieces.tolist(), 'actual_mode_energy_blocks': self.blocks,
            'nonbaseline_blocks': int((selected >= 0).sum()),
            'root_evaluations': self.calls, 'optimizer_success': bool(result.success),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block-size', type=int, choices=(2, 4, 8, 16), default=4)
    parser.add_argument('--ramifications', default='1,3')
    parser.add_argument('--lengths', default='64,128')
    parser.add_argument('--random-degree-cases', type=int, default=32)
    parser.add_argument('--iterations', type=int, default=45)
    parser.add_argument('--seed', type=int, default=720491)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ramifications = tuple(map(int, args.ramifications.split(',')))
    if any(P <= 0 or P % 2 == 0 for P in ramifications):
        raise ValueError('Use positive odd ramifications.')
    catalogue, sources = build_catalogue(args.block_size, ramifications,
                                         args.random_degree_cases, args.seed)
    payload = json.dumps(catalogue, separators=(',', ':'))
    records = []
    for n in map(int, args.lengths.split(',')):
        if n < args.block_size or n & (n - 1):
            raise ValueError('Use power-of-two physical lengths at least the block size.')
        for name, blocks in partitions(n, args.block_size, args.seed).items():
            started = monotonic()
            run = TiledPressure(n, args.block_size, catalogue, blocks)
            row = run.optimize(args.seed, args.iterations)
            row.update(n=n, partition=name, elapsed_seconds=monotonic() - started)
            records.append(row)
            print('TILED', n, args.block_size, name, row['rho'], row['nonbaseline_blocks'],
                  row['relative_log_cost'], row['elapsed_seconds'], flush=True)
    result = {
        'status': 'Exact polynomial catalogue and exploratory genuine-rho optimization; no inference of an all-family floor.',
        'block_size': args.block_size, 'source_ramifications': ramifications,
        'catalogue_profiles': len(catalogue), 'source_runs': sources,
        'catalogue_sha256': sha256(payload.encode()).hexdigest(), 'catalogue': catalogue,
        'seed': args.seed, 'optimizer_iterations': args.iterations,
        'floating_gain_tolerance': 1e-9,
        'strict_crossing_candidate_found': any(row['full_critical_pressure'] < -1e-9 for row in records),
        'optimized_runs': records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('DONE', args.block_size, len(catalogue), min(row['rho'] for row in records), flush=True)


if __name__ == '__main__':
    main()
