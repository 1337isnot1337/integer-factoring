#!/usr/bin/env python3
"""Bounded full-basis and full-gauge search for reflection crossing seeds.

Every physical basis and endpoint is exact over F4. Complete native Smith
profiles use exact coefficient arithmetic and floating valuation weights.
Only an independently replayed rational certificate can establish a new
crossing parameter. Completed optimization records are saved atomically.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize

from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import (
    F4_PRODUCT, flatten, identity_split, retract,
)
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from search_physical_code_gauges import start_params


def inverse(matrix):
    n = len(matrix)
    work = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if work[i][j]), None)
        if pivot is None:
            raise ValueError('singular physical basis')
        work[j], work[pivot] = work[pivot], work[j]
        inv = (0, 1, 3, 2)[work[j][j]]
        work[j] = [F4_PRODUCT[inv][value] for value in work[j]]
        for i in range(n):
            if i != j:
                coefficient = work[i][j]
                work[i] = [a ^ F4_PRODUCT[coefficient][b] for a, b in zip(work[i], work[j])]
    assert all(work[i][j] == int(i == j) for i in range(n) for j in range(n))
    return [row[n:] for row in work]


def transposed_inverse(matrix):
    return [list(column) for column in zip(*inverse(matrix))]


def change_basis(base, left, right):
    return retract(base, left, transposed_inverse(left), right, transposed_inverse(right))


def row_shear(matrix, target, source, scalar):
    result = [list(row) for row in matrix]
    result[target] = [a ^ F4_PRODUCT[scalar][b]
                      for a, b in zip(result[target], result[source])]
    return result


def canonical_basis(matrix):
    """Remove independent nonzero row scalings and symbol permutations."""
    rows = []
    for row in matrix:
        first = next(value for value in row if value)
        inv = (0, 1, 3, 2)[first]
        rows.append(tuple(F4_PRODUCT[inv][value] for value in row))
    return tuple(sorted(rows))


def generate_basis_pair(index, size, rng, elites):
    eye = identity_split(size)
    elementary = [(i, j, a) for i in range(size) for j in range(size)
                  if i != j for a in (1, 2, 3)]
    if index == 0:
        return 'identity', eye, eye
    if size == 4 and index <= 6:
        evaluation = []
        for value in range(4):
            row, power = [], 1
            for _ in range(4):
                row.append(power)
                power = F4_PRODUCT[power][value]
            evaluation.append(row)
        pascal = [[math.comb(i, j) % 2 if j <= i else 0 for j in range(4)] for i in range(4)]
        special = [(evaluation, evaluation), (inverse(evaluation), inverse(evaluation)),
                   (pascal, pascal), (evaluation, eye), (eye, evaluation),
                   (inverse(evaluation), eye)]
        left, right = special[index - 1]
        return 'evaluation-or-pascal', left, right
    offset = 7 if size == 4 else 1
    if index < offset + 3 * len(elementary):
        i, j, a = elementary[(index - offset) // 3]
        changed = row_shear(eye, i, j, a)
        mode = (index - offset) % 3
        return ('single-shear', changed if mode != 1 else eye,
                changed if mode != 0 else eye)
    if elites and index % 3:
        chosen = rng.choice(elites)
        left, right = [list(r) for r in chosen['basis_x']], [list(r) for r in chosen['basis_y']]
        changes = rng.randint(1, 3)
        label = 'elite-neighbor'
    else:
        left, right = eye, eye
        changes = rng.randint(size, 3 * size)
        label = 'independent-shears'
    for _ in range(changes):
        i, j, a = rng.choice(elementary)
        if rng.randrange(2):
            left = row_shear(left, i, j, a)
        else:
            right = row_shear(right, i, j, a)
    if index % 11 == 0:
        right = [list(row) for row in left]
        label += '-paired'
    return label, left, right


def packed_polynomial(value):
    low, high = value
    return sum((((low >> i) & 1) + 2 * ((high >> i) & 1)) << (2 * i)
               for i in range(max(low.bit_length(), high.bit_length())))


def exact_controls(base):
    """Compare native spectra with independent rational fraction-free elimination."""
    from native_matrix_smith_f4 import MatrixOracle
    from weighted_smith_oracle_f4 import WeightedSmithOracle
    size = base['q_x']
    eye = identity_split(size)
    left = row_shear(row_shear(eye, 0, size - 1, 2), 1, 0, 3)
    right = row_shear(eye, size - 1, 0, 1)
    records = []
    for name, seed in [('identity', base), ('nontrivial-basis', change_basis(base, left, right))]:
        for cut in (0, 1):
            matrix, _, _ = flatten(seed, cut)
            n = len(matrix)
            native = MatrixOracle(matrix, n)
            rw = [Fraction((5 * i + 1) % 7 - 3, 2) for i in range(n)]
            cw = [Fraction((3 * i + 4) % 5 - 2, 2) for i in range(n)]
            for infinity in (False, True):
                oracle = WeightedSmithOracle([[packed_polynomial(p) for p in row] for row in matrix],
                                             infinity=infinity)
                exact, _ = oracle.profile(rw, cw, Fraction(3, 2))
                observed = native.profile(rw, cw, scale=1.5, infinity=infinity)['orders']
                assert sorted(observed) == sorted(map(float, exact))
                records.append({'basis': name, 'cut': cut, 'infinity': infinity,
                                'orders': [[q.numerator, q.denominator] for q in exact]})
    return {'complete_endpoint_controls': 2 * size ** 4, 'exact_profile_controls': records}


def optimize(seed, starts, iterations, cut, counters):
    oracle = RectangularOracle(seed['values'], seed['q_x'], seed['q_u'])
    valid_best = [math.inf, None]

    def objective(params):
        counters['objective_calls'] += 1
        try:
            rho, px, py, *_ = oracle.types(params, cut)
            entropy = -sum(value * math.log(value) for value in px + py if value > 0)
            if entropy < 1e-9:
                counters['degenerate_type_calls'] += 1
                return 1000000.0, np.zeros(len(params))
            _, gradient = oracle.evaluate(params, px, py, cut)
            if rho < valid_best[0]:
                valid_best[:] = [rho, np.asarray(params, float).copy()]
            return rho, np.asarray(gradient) / (2 * entropy)
        except (ArithmeticError, ValueError):
            counters['unresolved_calls'] += 1
            return 1000000.0, np.zeros(len(params))

    records = []
    for name, start in starts:
        valid_best[:] = [math.inf, None]
        initial = objective(start)[0]
        result = minimize(objective, start, jac=True, method='SLSQP',
                          bounds=[(-24, 24)] * (len(start) - 1) + [(0.0002, 2)],
                          options={'maxiter': iterations, 'ftol': 2e-11})
        if valid_best[1] is None:
            records.append({'start': name, 'status': 'no resolved complete profile',
                            'rho': None, 'evaluations': int(result.nfev), 'cut': cut})
            continue
        selected = valid_best[1]
        rho, px, py, orders, gx, gy = oracle.types(selected, cut)
        records.append({'start': name, 'initial_rho': float(initial), 'rho': float(rho),
                        'params': list(map(float, selected)), 'type_x': px, 'type_y': py,
                        'scaled_orders': orders, 'scaled_energy_x': gx, 'scaled_energy_y': gy,
                        'iterations': int(result.nit), 'evaluations': int(result.nfev),
                        'success': bool(result.success), 'message': str(result.message),
                        'retained_initial': bool(initial <= result.fun), 'cut': cut})
    resolved = [record for record in records if record['rho'] is not None]
    if not resolved:
        raise ArithmeticError('every start failed to produce a complete profile')
    return min(resolved, key=lambda item: item['rho']), records


def write_atomic(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=4)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--iterations', type=int, default=240)
    parser.add_argument('--starts', type=int, default=3)
    parser.add_argument('--seed', type=int, default=311704)
    parser.add_argument('--cases', type=int, default=1000000)
    parser.add_argument('--controls-only', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.interval not in (4, 8) or args.seconds <= 0 or min(args.iterations, args.starts, args.cases) < 1:
        parser.error('positive budgets and interval4 or8 required')
    base = tensor(args.interval)
    size = base['q_x']
    controls = exact_controls(base) if args.interval == 4 else {'scope': 'large run uses previously checked engine'}
    if args.controls_only:
        write_atomic(args.output, controls)
        print(json.dumps({'exact_profile_controls': len(controls['exact_profile_controls'])}), flush=True)
        return
    energy_bits = physical_bit_energies(args.interval)
    energy = [sum(value for bit, value in enumerate(energy_bits) if symbol >> bit & 1)
              for symbol in range(size)]
    scale = 0.2 if args.interval == 4 else 0.09
    natural = start_params(np.asarray(energy) * scale, scale)
    rng = random.Random(args.seed)
    numerical = np.random.default_rng(args.seed)
    names = [Path(__file__).name, 'physical_split_crossings.py', 'reflection_dual_rail_gaussian.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py']
    report = {'schema_version': 1, 'scope': 'bounded numerical full physical-basis search; no exact pressure or global optimality claim',
              'interval': args.interval, 'alphabet': size, 'seed': args.seed,
              'requested_seconds': args.seconds, 'iteration_limit': args.iterations,
              'starts_per_basis': args.starts, 'controls': controls, 'cases': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'counters': {'objective_calls': 0, 'unresolved_calls': 0, 'degenerate_type_calls': 0}}
    seen, elites, best = set(), [], None
    started, index = time.monotonic(), 0
    while time.monotonic() - started < args.seconds and len(report['cases']) < args.cases:
        label, left, right = generate_basis_pair(index, size, rng, elites)
        index += 1
        key = (canonical_basis(left), canonical_basis(right))
        if key in seen:
            continue
        seen.add(key)
        seed = change_basis(base, left, right)
        cut = 0 if report['cases'] and len(report['cases']) % 7 == 0 else 1
        starts = [('natural', natural.copy())]
        if best is not None and args.starts >= 2:
            starts.append(('previous-best', np.asarray(best['best']['params'], float)))
        while len(starts) < args.starts:
            start = natural + numerical.normal(0, 0.12 if len(starts) % 2 else 0.45, len(natural))
            start[-1] = min(1.2, max(0.005, scale * math.exp(numerical.normal(0, 0.3))))
            starts.append(('perturbed-natural', start))
        local, runs = optimize(seed, starts, args.iterations, cut, report['counters'])
        record = {'index': index - 1, 'kind': label, 'basis_x': left, 'basis_y': right,
                  'degree': seed['degree'], 'endpoint_entries': size ** 4, 'best': local,
                  'runs': runs, 'elapsed_seconds': time.monotonic() - started}
        report['cases'].append(record)
        elites.append(record)
        elites.sort(key=lambda item: item['best']['rho'])
        del elites[12:]
        if best is None or local['rho'] < best['best']['rho']:
            best = record
        report.update(best_rho=best['best']['rho'], best_case_index=best['index'],
                      completed_cases=len(report['cases']), elapsed_seconds=time.monotonic() - started,
                      completed=False)
        write_atomic(args.output, report)
        print(json.dumps({'cases': len(report['cases']), 'best_rho': report['best_rho'],
                          'current_rho': local['rho'], 'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        if report['best_rho'] < 1 / math.sqrt(2) - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate requires exact independent replay'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'requested time or case budget completed'))
    write_atomic(args.output, report)


if __name__ == '__main__':
    main()
