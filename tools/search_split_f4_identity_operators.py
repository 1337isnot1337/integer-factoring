#!/usr/bin/env python3
"""Search explicitly finite F4 projector kernels in three split identity modes.

The profiles are exact. Pressure minimization is numerical and is search
evidence only. The catalogue includes every rank-one and rank-two idempotent
in each three-dimensional occupation sector, with four prescribed couplings
between sectors, all four top-state colors, and all four endpoint flags.
"""
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

from compiled_f4_identity_operators import compile_family, combine_kernel
from coupled_identity_operator_filters import _geometry, identity_kernel, states
from f4_identity_operator_filters import zero_profiles, infinity_profiles
from f4_repaired_identity_filters import MUL, INV
from search_coupled_identity_operator_filters import Pressure, encode_kernel


def dot(first, second):
    result = 0
    for a, b in zip(first, second):
        result ^= MUL[a][b]
    return result


def matrix_product(first, second):
    return tuple(tuple(dot(row, column) for column in zip(*second)) for row in first)


@lru_cache(maxsize=1)
def projectors():
    """All 674 idempotent 3x3 matrices over F4, independently count checked."""
    zero = ((0, 0, 0),) * 3
    identity = tuple(tuple(int(i == j) for j in range(3)) for i in range(3))
    rank_one = set()
    for left in product(range(4), repeat=3):
        leading = next((value for value in left if value), None)
        if leading != 1:
            continue
        for right in product(range(4), repeat=3):
            if dot(left, right) == 1:
                rank_one.add(tuple(tuple(MUL[a][b] for b in right) for a in left))
    assert len(rank_one) == 336
    result = {zero, identity} | rank_one
    result |= {tuple(tuple(value ^ identity[i][j] for j, value in enumerate(row))
                     for i, row in enumerate(matrix)) for matrix in rank_one}
    assert len(result) == 674
    assert all(matrix_product(matrix, matrix) == matrix for matrix in result)
    return tuple(sorted(result))


def reverse(matrix, transpose=False):
    return tuple(tuple(matrix[2-(j if transpose else i)][2-(i if transpose else j)]
                       for j in range(3)) for i in range(3))


@lru_cache(maxsize=1)
def sector_pairs():
    zero = ((0, 0, 0),) * 3
    result = set()
    for matrix in projectors():
        result.add((matrix, zero))
        result.add((zero, matrix))
        result.add((matrix, reverse(matrix)))
        result.add((matrix, reverse(matrix, True)))
    return tuple(sorted(result))


def affine_sources(n, ramification, geometry, one_height, two_height, top_height):
    selected = _geometry(2, n, 3, geometry)['selected']
    alphabet = states(2, 3)
    sectors = tuple(tuple(state for state in alphabet if sum(state) == level) for level in (1, 2))
    perturbations = []
    for labels, height in zip(sectors, (one_height, two_height)):
        perturbations.extend({(row, column): {0: 1, -height: 1}}
                             for row, column in product(labels, repeat=2))
    top = (1, 1, 1)
    perturbations.append({(top, top): {0: 1, -top_height: 1}})
    flag = {}
    for labels in sectors:
        for row, column in product(labels, repeat=2):
            if sum(j * digit for j, digit in zip(selected, row)) > sum(
                    j * digit for j, digit in zip(selected, column)):
                flag[row, column] = {0: 1}
    perturbations.append(flag)
    constant = identity_kernel(2, 3)
    return constant, tuple(perturbations)


def catalogue():
    for first, second in sector_pairs():
        flat = tuple(value for matrix in (first, second) for row in matrix for value in row)
        for top_color, endpoint_flag in product(range(4), repeat=2):
            yield flat + (top_color, endpoint_flag)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--P', type=int, default=1)
    parser.add_argument('--geometry', default='split:1:2')
    parser.add_argument('--one-height', type=int, default=1)
    parser.add_argument('--two-height', type=int, default=1)
    parser.add_argument('--top-height', type=int, default=1)
    parser.add_argument('--refinements', type=int, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.one_height, args.two_height, args.top_height) <= 0:
        raise ValueError('The three slow-pole heights must be positive.')
    start = time.monotonic()
    constant, perturbations = affine_sources(args.n, args.P, args.geometry,
                                             args.one_height, args.two_height, args.top_height)
    compiled = compile_family(args.n, args.P, 3, constant, perturbations, args.geometry)
    grid = tuple(product(math.log(2 * args.n) / args.n * np.array((.25, .5, .75, 1., 1.5, 2.)),
                         sorted(set((-2*args.n, -args.n, -4, -2, -1, 0, 1, 2, 4,
                                     args.n/2, args.n, 2*args.n)))))
    scales = np.array([point[0] for point in grid])
    charges = np.array([point[1] for point in grid])
    classes, completed, best = {}, 0, math.inf
    for colors in catalogue():
        blocks, infinity = compiled(colors)
        key = (tuple((row['charge'], row['orders']) for row in blocks), infinity)
        completed += 1
        if key in classes:
            classes[key]['rows'] += 1
        else:
            pressure = Pressure(args.n, args.P, 3, args.geometry, blocks, infinity)
            values = pressure.grid(scales, charges)
            relative = pressure.grid(scales, charges, True)
            i, j = int(np.argmin(values)), int(np.argmin(relative))
            record = {'rows': 1, 'colors': colors, 'zero_blocks': blocks, 'infinity_active': infinity,
                      'critical_grid': {'pressure': float(values[i]), 'scale': float(scales[i]),
                                        'charge': float(charges[i])},
                      'relative_grid': {'change': float(relative[j]), 'scale': float(scales[j]),
                                        'charge': float(charges[j])}}
            classes[key] = record
            if values[i] < best - 1e-10:
                best = float(values[i])
                print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                                  'best_critical_grid': record['critical_grid']}), flush=True)
        if completed % 4096 == 0:
            print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                              'seconds': time.monotonic()-start}), flush=True)
    rows = sorted(classes.values(), key=lambda row: row['critical_grid']['pressure'])
    for row in rows[:args.refinements]:
        pressure = Pressure(args.n, args.P, 3, args.geometry, row['zero_blocks'], row['infinity_active'])
        point = row['critical_grid']
        optimum = minimize(lambda value: pressure(*value), [math.log(point['scale']), point['charge']],
                           method='Nelder-Mead', bounds=[(-math.log(args.n)-4, math.log(1.5)),
                                                        (-2*args.n, 2*args.n)],
                           options={'maxiter': 200, 'fatol': 1e-10})
        refined = {'pressure': float(optimum.fun), 'scale': math.exp(float(optimum.x[0])),
                   'charge': float(optimum.x[1]), 'success': bool(optimum.success)}
        row['refined'] = refined if refined['pressure'] <= point['pressure'] else dict(point)
    exact_replays = []
    # These reconstruct the original Laurent kernel and bypass the bit-plane compiler.
    for row in rows[:5]:
        kernel = combine_kernel(constant, perturbations, row['colors'])
        assert zero_profiles(args.n, args.P, 3, kernel, args.geometry) == row['zero_blocks']
        assert infinity_profiles(args.n, args.P, 3, kernel, args.geometry) == row['infinity_active']
        exact_replays.append({'colors': row['colors'], 'kernel': encode_kernel(kernel)})
    expected = len(sector_pairs()) * 16
    assert completed == expected == sum(row['rows'] for row in rows)
    result = {'status': 'Complete specified F4 projector catalogue with exact profiles and numerical pressure optimization; not a global obstruction.',
              'n': args.n, 'ramification': args.P, 'geometry': args.geometry,
              'one_height': args.one_height, 'two_height': args.two_height, 'top_height': args.top_height,
              'idempotents': len(projectors()), 'sector_pair_count': len(sector_pairs()),
              'completed_candidates': completed, 'complete': True, 'profile_class_count': len(rows),
              'grid_size': len(grid),
              'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in rows),
              'best_grid_relative_change': min(row['relative_grid']['change'] for row in rows),
              'exact_best_profile_replays': exact_replays,
              'profile_classes': rows, 'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('profile_classes', 'exact_best_profile_replays')}), flush=True)


if __name__ == '__main__':
    main()
