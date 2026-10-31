#!/usr/bin/env python3
"""Reproducible coordinate and coherent physical-code gauge searches.

All splits and local polynomial coefficients are exact. The full weighted
Smith oracle uses exact F4 arithmetic with floating valuation weights.
Scale, Gibbs types, and gauge optimization are numerical search evidence.
"""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import time

from physical_split_crossings import SplitProfile, check_split, retract
from reflection_dual_rail_gaussian import physical_bit_energies, tensor


def odd_fiber_maps(groups, dimension):
    groups = [tuple(group) for group in groups]
    flattened = [symbol for group in groups for symbol in group]
    if (not groups or any(not group or len(group) % 2 == 0 for group in groups)
            or len(set(flattened)) != len(flattened)
            or any(not isinstance(symbol, int) or not 0 <= symbol < dimension
                   for symbol in flattened)):
        raise ValueError('Use nonempty, odd, disjoint fibers inside the ambient alphabet.')
    maps = [[int(symbol in group) for symbol in range(dimension)] for group in groups]
    check_split(maps, maps)
    return maps


def energies():
    weights = physical_bit_energies(8)
    return [sum(weight for bit, weight in enumerate(weights) if symbol >> bit & 1)
            for symbol in range(16)]


def start_params(energy, log_scale):
    """Input energies already contain log(s); letter zero is the reference."""
    import numpy as np

    size = len(energy)
    shifted = np.asarray(energy, float) - energy[0]
    params = np.zeros(8 * (size - 1) + 1)
    for end, legs in [(0, (2, 3)), (1, (0, 1))]:
        for leg in legs:
            offset = (4 * end + leg) * (size - 1)
            params[offset:offset + size - 1] = shifted[1:] / 2
    params[-1] = log_scale
    return params


def optimize(seed, starts, maximum_iterations):
    import math
    import numpy as np
    from scipy.optimize import minimize
    from native_rectangular_smith_f4 import RectangularOracle

    size = seed['q_x']
    oracle = RectangularOracle(seed['values'], size, size)
    records = []

    def objective(params):
        try:
            rho, px, py, *_ = oracle.types(params)
            entropy = -sum(value * math.log(value) for value in px + py if value > 0)
            if entropy < 1e-9:
                return 1., np.zeros(len(params))
            _, gradient = oracle.evaluate(params, px, py)
            return rho, np.asarray(gradient) / (2 * entropy)
        except (ArithmeticError, ValueError):
            return 2., np.zeros(len(params))

    for name, start in starts:
        initial = objective(start)[0]
        result = minimize(objective, start, jac=True, method='SLSQP',
                          bounds=[(-20, 20)] * (len(start) - 1) + [(.0001, 2)],
                          options={'maxiter': maximum_iterations, 'ftol': 1e-11})
        use_initial = not result.fun < initial
        params = start if use_initial else result.x
        rho, px, py, orders, gx, gy = oracle.types(params)
        records.append({
            'start': name, 'initial_rho': initial, 'rho': rho,
            'params': list(map(float, params)), 'type_x': px, 'type_y': py,
            'energy_x_times_log_s': gx, 'energy_y_times_log_s': gy,
            'orders_times_log_s': orders, 'evaluations': result.nfev,
            'optimizer_success': bool(result.success), 'optimizer_message': str(result.message),
            'used_initial_candidate': use_initial,
        })
    return min(records, key=lambda record: record['rho']), records


def coordinate_candidates(count, rng):
    candidates = [(0, 2, 4, 8, 12), (0, 1, 2, 4, 8), (0, 4, 8, 12, 14),
                  (0, 2, 6, 8, 10), (0, 3, 5, 6, 9)]
    while len(candidates) < count:
        alphabet = tuple(sorted([0] + rng.sample(range(1, 16), 4)))
        if alphabet not in candidates:
            candidates.append(alphabet)
    for alphabet in candidates[:count]:
        yield 'coordinate', tuple((symbol,) for symbol in alphabet)


def coherent_candidates(random_count, rng, energy):
    alphabet = (0, 2, 4, 8, 12)
    omitted = [symbol for symbol in range(16) if symbol not in alphabet]
    candidates = [('coordinate', [(symbol,) for symbol in alphabet])]
    for a, b in combinations(omitted, 2):
        candidates.append(('three-letter-vacuum', [(0, a, b)]
                           + [(symbol,) for symbol in alphabet[1:]]))
    ordered = sorted(range(16), key=lambda symbol: energy[symbol])
    candidates.append(('energy-triples', [tuple(ordered[3 * i:3 * i + 3]) for i in range(5)]))
    candidates.append(('consecutive-triples', [tuple(range(3 * i, 3 * i + 3)) for i in range(5)]))
    for _ in range(random_count):
        letters = list(range(16))
        rng.shuffle(letters)
        candidates.append(('random-triples', [tuple(sorted(letters[3 * i:3 * i + 3]))
                                               for i in range(5)]))
    seen = set()
    for name, groups in candidates:
        groups = tuple(sorted((tuple(sorted(group)) for group in groups),
                              key=lambda group: (min(energy[symbol] for symbol in group), group)))
        if groups not in seen:
            seen.add(groups)
            yield name, groups


def compressed_candidates(random_count, rng, energy):
    ordered = sorted(range(16), key=lambda symbol: energy[symbol])
    triples = [tuple(sorted(ordered[i:i + 3])) for i in range(14)]
    triples += list(combinations(ordered[-6:], 3))
    triples += [tuple(sorted(rng.sample(range(16), 3))) for _ in range(random_count)]
    for triple in dict.fromkeys(tuple(sorted(group)) for group in triples):
        groups = [triple] + [(symbol,) for symbol in range(16) if symbol not in triple]
        groups.sort(key=lambda group: min(energy[symbol] for symbol in group))
        yield 'one-triple', tuple(groups)


def main():
    import numpy as np

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', required=True, choices=('coordinate', 'coherent', 'compressed'))
    parser.add_argument('--count', type=int, default=30)
    parser.add_argument('--random-count', type=int)
    parser.add_argument('--refine', type=int, default=20)
    parser.add_argument('--seed', type=int, default=20261031)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.count < 1 or args.refine < 0 or (args.random_count is not None and args.random_count < 0):
        parser.error('Counts must be nonnegative, with --count positive.')
    start = time.time()
    base, energy = tensor(8), energies()
    rng = random.Random(args.seed)
    if args.suite == 'coordinate':
        candidates = coordinate_candidates(args.count, rng)
    elif args.suite == 'coherent':
        candidates = coherent_candidates(100 if args.random_count is None else args.random_count,
                                         rng, energy)
    else:
        candidates = compressed_candidates(30 if args.random_count is None else args.random_count,
                                           rng, energy)
    records, refinements, seeds = [], [], {}
    endpoint_entries = 0
    source_names = [Path(__file__).name, 'physical_split_crossings.py',
                    'reflection_dual_rail_gaussian.py', 'native_rectangular_smith_f4.py',
                    'native_rectangular_smith_f4.cpp', 'native_matrix_smith_f4.py',
                    'native_matrix_smith_f4.cpp']
    hashes = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in source_names}

    def save():
        best = min([record['best']['rho'] for record in records]
                   + [record['best']['rho'] for record in refinements])
        report = {'suite': args.suite, 'seed': args.seed, 'case_count': len(records),
                  'endpoint_entries': endpoint_entries, 'cases': records,
                  'refinements': refinements, 'best_rho': best,
                  'elapsed_seconds': time.time() - start, 'source_sha256': hashes,
                  'scope': 'Exact split maps and endpoint checks, complete polynomial profiles, numerical gauges and Gibbs roots only; no optimality or threshold certificate.'}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({key: report[key] for key in ('suite', 'case_count', 'best_rho',
                                                       'elapsed_seconds')}
                         | {'refined': len(refinements)}), flush=True)

    for name, groups in candidates:
        maps = odd_fiber_maps(groups, 16)
        seed = retract(base, maps, maps)
        endpoint_entries += len(seed['values'])
        if args.suite == 'coordinate':
            selected = np.asarray([energy[group[0]] for group in groups], float)
            target = np.zeros(5)
            target[np.argsort(selected)] = [0, 0, 2.0346278387877708,
                                            2.6290972704856848, 3.03102987564101]
            starts = [('natural', start_params(selected * .22, .22)),
                      ('target', start_params(target, .12)),
                      ('uniform', start_params(np.zeros(5), .08))]
            best, runs = optimize(seed, starts, 250)
            records.append({'name': name, 'groups': groups, 'best': best, 'runs': runs})
        else:
            oracle = SplitProfile(seed)
            height_cases = [
                ('minimum', [min(energy[symbol] for symbol in group) for group in groups]),
                ('mean', [sum(energy[symbol] for symbol in group) / len(group) for group in groups]),
            ]
            height_cases.append(('uniform', [0] * len(groups)) if args.suite == 'coherent' else
                                ('maximum', [max(energy[symbol] for symbol in group) for group in groups]))
            values = []
            for method, selected in height_cases:
                result = oracle.pressure(selected)
                values.append({'method': method, 'energy': selected,
                               'rho': result['rho'], 'log_scale': result['log_scale'],
                               'profiles': result['profiles']})
            records.append({'name': name, 'groups': groups, 'best': min(values, key=lambda item: item['rho']),
                            'height_cases': values})
            seeds[groups] = seed
        if len(records) % 10 == 0:
            save()

    if args.suite == 'coherent':
        records.sort(key=lambda record: record['best']['rho'])
        for record in records[:args.refine]:
            initial = record['best']
            selected = np.asarray(initial['energy'])
            target = np.zeros(5)
            target[np.argsort(selected)] = [0, 0, 2, 8 / 3, 3]
            starts = [('static-best', start_params(selected * initial['log_scale'], initial['log_scale'])),
                      ('target', start_params(target, .12))]
            best, runs = optimize(seeds[record['groups']], starts, 350)
            refinements.append({'name': record['name'], 'groups': record['groups'],
                                'best': best, 'runs': runs})
            save()
    save()


if __name__ == '__main__':
    main()
