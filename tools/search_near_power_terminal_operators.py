#!/usr/bin/env python3
"""Finite complete F4 projector screen at n=2^k-1 terminal modes."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from compiled_f4_identity_operators import combine_kernel
from coupled_identity_operator_filters import identity_kernel, states
from near_power_terminal_operators import (
    geometry_data, infinity_profiles, zero_profiles as full_zero_profiles, TerminalPressure,
)
from near_power_terminal_weyl import compile_family, zero_profiles
from search_coupled_identity_operator_filters import encode_kernel
from search_split_f4_identity_operators import catalogue, projectors, sector_pairs


def affine_sources(n, one_height=1, two_height=1, top_height=1, initial_modes=0):
    selected = geometry_data(n, 3, initial_modes)['selected']
    alphabet = states(2, 3)
    sectors = tuple(tuple(state for state in alphabet if sum(state) == level) for level in (1, 2))
    perturbations = []
    for labels, height in zip(sectors, (one_height, two_height)):
        perturbations.extend({(row, column): {0: 1, -height: 1}}
                             for row, column in product(labels, repeat=2))
    perturbations.append({((1, 1, 1), (1, 1, 1)): {0: 1, -top_height: 1}})
    flag = {}
    for labels in sectors:
        for row, column in product(labels, repeat=2):
            if sum(j*digit for j, digit in zip(selected, row)) > sum(
                    j*digit for j, digit in zip(selected, column)):
                flag[row, column] = {0: 1}
    perturbations.append(flag)
    return identity_kernel(2, 3), tuple(perturbations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=31)
    parser.add_argument('--P', type=int, default=1)
    parser.add_argument('--initial-modes', type=int, choices=(0, 1), default=0)
    parser.add_argument('--one-height', type=int, default=1)
    parser.add_argument('--two-height', type=int, default=1)
    parser.add_argument('--top-height', type=int, default=1)
    parser.add_argument('--refinements', type=int, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.one_height, args.two_height, args.top_height) <= 0:
        raise ValueError('Use positive slow-pole heights.')
    start = time.monotonic()
    constant, perturbations = affine_sources(args.n, args.one_height, args.two_height, args.top_height, args.initial_modes)
    compiled = compile_family(args.n, args.P, 3, constant, perturbations, args.initial_modes)
    grid = tuple(product(math.log(2*args.n)/args.n*np.array((.25, .5, .75, 1., 1.5, 2.)),
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
            pressure = TerminalPressure(args.n, args.P, 3, blocks, infinity, args.initial_modes)
            values, relative = pressure.grid(scales, charges), pressure.grid(scales, charges, True)
            i, j = int(np.argmin(values)), int(np.argmin(relative))
            row = {'rows': 1, 'colors': colors, 'zero_blocks': blocks, 'infinity_active': infinity,
                   'critical_grid': {'pressure': float(values[i]), 'scale': float(scales[i]), 'charge': float(charges[i])},
                   'relative_grid': {'change': float(relative[j]), 'scale': float(scales[j]), 'charge': float(charges[j])}}
            classes[key] = row
            if values[i] < best-1e-10:
                best = float(values[i])
                print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                                  'best_critical_grid': row['critical_grid']}), flush=True)
        if completed % 4096 == 0:
            print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                              'seconds': time.monotonic()-start}), flush=True)
    rows = sorted(classes.values(), key=lambda row: row['critical_grid']['pressure'])
    for row in rows[:args.refinements]:
        pressure = TerminalPressure(args.n, args.P, 3, row['zero_blocks'], row['infinity_active'], args.initial_modes)
        point = row['critical_grid']
        optimum = minimize(lambda value: pressure(*value), [math.log(point['scale']), point['charge']],
                           method='Nelder-Mead', bounds=[(-math.log(args.n)-4, math.log(1.5)),
                                                        (-2*args.n, 2*args.n)],
                           options={'maxiter': 200, 'fatol': 1e-10})
        refined = {'pressure': float(optimum.fun), 'scale': math.exp(float(optimum.x[0])),
                   'charge': float(optimum.x[1]), 'success': bool(optimum.success)}
        row['refined'] = refined if refined['pressure'] <= point['pressure'] else dict(point)
    replays = []
    for row in rows[:5]:
        kernel = combine_kernel(constant, perturbations, row['colors'])
        assert zero_profiles(args.n, args.P, 3, kernel, args.initial_modes) == row['zero_blocks']
        assert infinity_profiles(args.n, args.P, 3, kernel, args.initial_modes) == row['infinity_active']
        replays.append({'colors': row['colors'], 'kernel': encode_kernel(kernel)})
    # Independently replay the complete rectangular core, without the dummy transform.
    assert full_zero_profiles(args.n, args.P, 3, combine_kernel(constant, perturbations, rows[0]['colors']), args.initial_modes) == rows[0]['zero_blocks']
    assert completed == len(sector_pairs())*16 == sum(row['rows'] for row in rows)
    result = {'status': 'Complete specified F4 projector catalogue, exact profiles, numerical pressure search; no global obstruction.',
              'n': args.n, 'ramification': args.P, 'initial_modes': args.initial_modes, 'selected': geometry_data(args.n, 3, args.initial_modes)['selected'],
              'retained': geometry_data(args.n, 3, args.initial_modes)['retained'],
              'one_height': args.one_height, 'two_height': args.two_height, 'top_height': args.top_height,
              'idempotents': len(projectors()), 'sector_pairs': len(sector_pairs()),
              'completed_candidates': completed, 'complete': True, 'profile_class_count': len(rows),
              'grid_size': len(grid),
              'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in rows),
              'best_relative_grid_change': min(row['relative_grid']['change'] for row in rows),
              'exact_best_profile_replays': replays, 'profile_classes': rows,
              'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('profile_classes', 'exact_best_profile_replays')}), flush=True)


if __name__ == '__main__':
    main()
