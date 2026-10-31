#!/usr/bin/env python3
"""A bounded complete F3 projector catalogue in two near-power directions."""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from coupled_identity_mode_filters import add
from coupled_identity_operator_filters import identity_kernel, states
from prime_deficit_identity_operators import zero_profiles, infinity_profiles, Pressure


@lru_cache(maxsize=2)
def projectors(size):
    identity = tuple(tuple(int(i == j) for j in range(size)) for i in range(size))
    zero = ((0,)*size,)*size
    result = {zero, identity}
    for left in product(range(3), repeat=size):
        if next((value for value in left if value), 0) != 1:
            continue
        for right in product(range(3), repeat=size):
            if sum(a*b for a, b in zip(left, right)) % 3 != 1:
                continue
            first = tuple(tuple(a*b % 3 for b in right) for a in left)
            result.add(first)
            result.add(tuple(tuple((identity[i][j]-first[i][j]) % 3 for j in range(size)) for i in range(size)))
    assert len(result) == (14 if size == 2 else 236)
    assert all(all(sum(matrix[i][k]*matrix[k][j] for k in range(size)) % 3 == matrix[i][j]
                       for i in range(size) for j in range(size)) for matrix in result)
    return tuple(sorted(result))


def paired_catalogue(height=1):
    alphabet = states(3, 2)
    sectors = {level: tuple(state for state in alphabet if sum(state) == level) for level in (1, 2, 3)}
    choices = [((sectors[2], matrix),) for matrix in projectors(3)]
    choices += [((sectors[1], matrix),
                 (sectors[3], tuple(tuple(matrix[1-i][1-j] for j in range(2)) for i in range(2))))
                for matrix in projectors(2)]
    result = {}
    for blocks, top_pole in product(choices, (False, True)):
        kernel = identity_kernel(3, 2)
        for labels, matrix in blocks:
            for i, row in enumerate(labels):
                for j, column in enumerate(labels):
                    if matrix[i][j]:
                        index = (row, column)
                        kernel[index] = add(kernel.get(index, {}), {0: 2, -height: 1}, 3, matrix[i][j])
        if top_pole:
            kernel[(2, 2), (2, 2)] = {-height: 1}
        kernel = {index: value for index, value in kernel.items() if value}
        key = tuple((a, b, tuple(sorted(value.items()))) for (a, b), value in sorted(kernel.items()))
        result[key] = kernel
    assert len(result) == 498
    return tuple(result.values())


def rank_one_catalogue(n, P):
    alphabet = states(3, 2)
    result = [identity_kernel(3, 2)]
    counts = {}
    for occupation in (1, 2, 3):
        labels = tuple(state for state in alphabet if sum(state) == occupation)
        rank_one = tuple(matrix for matrix in projectors(len(labels))
                         if sum(matrix[i][i] for i in range(len(labels))) % 3 == 1)
        counts[occupation] = len(rank_one)
        for matrix, coefficient, exponent in product(rank_one, (1, 2), (-n*P, -1, 1, n*P)):
            kernel = identity_kernel(3, 2)
            for i, row in enumerate(labels):
                for j, column in enumerate(labels):
                    if matrix[i][j]:
                        index = (row, column)
                        kernel[index] = add(kernel.get(index, {}), {0: 1, exponent: 2}, 3,
                                            coefficient*matrix[i][j])
            result.append({index: value for index, value in kernel.items() if value})
    assert counts == {1: 12, 2: 117, 3: 12} and len(result) == 1129
    keys = {tuple((a, b, tuple(sorted(value.items()))) for (a, b), value in sorted(kernel.items()))
            for kernel in result}
    assert len(keys) == len(result)
    return tuple(result)


def encode(kernel):
    return [{'row': row, 'column': column, 'polynomial': value}
            for (row, column), value in sorted(kernel.items())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=int, default=27)
    parser.add_argument('--initial-modes', type=int, choices=(0, 1), default=0)
    parser.add_argument('--P', type=int, default=1)
    parser.add_argument('--height', type=int, default=1)
    parser.add_argument('--family', choices=('rank-one', 'paired'), default='rank-one')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.height < 1:
        raise ValueError('Use a positive Laurent pole height.')
    n, initial, terminal = args.parent-1, args.initial_modes, 2-args.initial_modes
    geometry = (3, args.parent, 1, initial, terminal, args.P)
    start = time.monotonic()
    grid = tuple(product(math.log(2*n)/n*np.array((.25, .5, .75, 1., 1.5, 2.)),
                         sorted(set((-2*n, -n, -4, -2, -1, 0, 1, 2, 4, n/2, n, 2*n)))))
    scales = np.array([value[0] for value in grid])
    charges = np.array([value[1] for value in grid])
    classes, best = {}, math.inf
    candidates = rank_one_catalogue(n, args.P) if args.family == 'rank-one' else paired_catalogue(args.height)
    for completed, kernel in enumerate(candidates, 1):
        blocks = zero_profiles(*geometry, kernel)
        infinity = infinity_profiles(*geometry, kernel)
        key = (tuple((block['charge'], block['orders']) for block in blocks), infinity)
        if key in classes:
            classes[key]['count'] += 1
        else:
            pressure = Pressure(*geometry, blocks, infinity)
            value, relative = pressure.grid(scales, charges), pressure.grid(scales, charges, True)
            i, j = int(np.argmin(value)), int(np.argmin(relative))
            record = {'count': 1, 'kernel': encode(kernel), 'zero_blocks': blocks, 'infinity_active': infinity,
                      'critical_grid': {'pressure': float(value[i]), 'scale': float(scales[i]), 'charge': float(charges[i])},
                      'relative_grid': {'change': float(relative[j]), 'scale': float(scales[j]), 'charge': float(charges[j])}}
            lower = min(degree for value in kernel.values() for degree in value)
            upper = max(degree for value in kernel.values() for degree in value)
            record['degree_accounting'] = {'common_Laurent_clearing': max(0, -lower),
                'coefficient_degree_after_clearing': upper+max(0, -lower),
                'Gaussian_degree_upper_bound': 8*n*n,
                'full_polynomial_degree_upper_bound': args.P*(8*n*n+8*n)+upper+max(0, -lower)}
            classes[key] = record
            if value[i] < best-1e-10:
                best = float(value[i])
                print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                                  'best_grid': record['critical_grid']}), flush=True)
        if completed % 50 == 0:
            print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                              'seconds': time.monotonic()-start}), flush=True)
    records = sorted(classes.values(), key=lambda row: row['critical_grid']['pressure'])
    for row in records[:30]:
        pressure = Pressure(*geometry, row['zero_blocks'], row['infinity_active'])
        point = row['critical_grid']
        fit = minimize(lambda value: float(pressure.grid([math.exp(float(value[0]))], [float(value[1])])[0]),
                       [math.log(point['scale']), point['charge']], method='Nelder-Mead',
                       bounds=[(-math.log(n)-4, math.log(1.5)), (-2*n, 2*n)],
                       options={'maxiter': 200, 'fatol': 1e-10})
        refined = {'pressure': float(fit.fun), 'scale': math.exp(float(fit.x[0])),
                   'charge': float(fit.x[1]), 'success': bool(fit.success)}
        row['refined'] = refined if refined['pressure'] <= point['pressure'] else dict(point)
    assert completed == len(candidates) == sum(row['count'] for row in records)
    result = {'status': 'Complete bounded F3 projector family; exact profiles and numerical pressure search, no global obstruction or improved bound claim.',
              'prime': 3, 'parent': args.parent, 'n': n, 'deficit': 1,
              'initial_modes': initial, 'terminal_modes': terminal, 'ramification': args.P, 'height': args.height,
              'family': args.family, 'rank_one_projectors_by_occupation': {1: 12, 2: 117, 3: 12},
              'rank_one_coefficients': (1, 2), 'rank_one_exponents': (-n*args.P, -1, 1, n*args.P),
              'completed_candidates': completed, 'complete': True, 'profile_classes': len(records),
              'grid_size': len(grid),
              'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in records),
              'best_relative_grid_change': min(row['relative_grid']['change'] for row in records),
              'records': records, 'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'records'}), flush=True)


if __name__ == '__main__':
    main()
