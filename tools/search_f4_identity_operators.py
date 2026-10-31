#!/usr/bin/env python3
"""Exhaust finite F4 two-mode operator blocks of a prescribed slow degree."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from compiled_f4_identity_operators import compile_family, combine_kernel
from coupled_identity_operator_filters import identity_kernel
from search_coupled_identity_operator_filters import Pressure, encode_kernel


def endpoint_blocks(policy):
    if policy == 'identity':
        return (((1, 0), (0, 1)),)
    return tuple(((a, 0), (b, d)) for a, b, d in product((1, 2, 3), range(4), (1, 2, 3)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--P', type=int, default=1)
    parser.add_argument('--degree', type=int, choices=(1, 2), default=1)
    parser.add_argument('--endpoint', choices=('identity', 'filtered'), default='filtered')
    parser.add_argument('--geometry', choices=('initial', 'mixed', 'terminal'), default='mixed')
    parser.add_argument('--top-exponent', type=int, default=0)
    parser.add_argument('--vacuum-color', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--refinements', type=int, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    middle = ((1, 0), (0, 1))
    perturbations = [{(a, b): {0: 1, -degree: 1}}
                     for degree in range(1, args.degree+1) for a, b in product(middle, repeat=2)]
    grid = tuple(product(math.log(2*args.n)/args.n*np.array((.25, .5, .75, 1., 1.5, 2.)),
                         sorted(set((-args.n, -4, -2, -1, 0, 1, 2, 4, args.n/2, args.n, 2*args.n)))))
    scales, charges = np.array([point[0] for point in grid]), np.array([point[1] for point in grid])
    classes, completed, best = {}, 0, math.inf
    for endpoint in endpoint_blocks(args.endpoint):
        constant = identity_kernel(2, 2)
        constant[(0, 0), (0, 0)] = {0: args.vacuum_color}
        constant[(1, 1), (1, 1)] = {-args.top_exponent: 1}
        for i, a in enumerate(middle):
            for j, b in enumerate(middle):
                if endpoint[i][j]:
                    constant[a, b] = {0: endpoint[i][j]}
                else:
                    constant.pop((a, b), None)
        compiled = compile_family(args.n, args.P, 2, constant, perturbations, args.geometry)
        for colors in product(range(4), repeat=len(perturbations)):
            blocks, infinity = compiled(colors)
            key = (tuple((row['charge'], row['orders']) for row in blocks), infinity)
            completed += 1
            if key in classes:
                classes[key]['rows'] += 1
                continue
            pressure = Pressure(args.n, args.P, 2, args.geometry, blocks, infinity)
            values = pressure.grid(scales, charges)
            relative = pressure.grid(scales, charges, True)
            i, j = int(np.argmin(values)), int(np.argmin(relative))
            record = {'rows': 1, 'endpoint': endpoint, 'colors': colors,
                      'zero_blocks': blocks, 'infinity_active': infinity,
                      'critical_grid': {'pressure': float(values[i]), 'scale': float(scales[i]),
                                        'charge': float(charges[i])},
                      'relative_grid': {'change': float(relative[j]), 'scale': float(scales[j]),
                                        'charge': float(charges[j])}}
            classes[key] = record
            if values[i] < best-1e-10:
                best = float(values[i])
                print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                                  'best_critical_grid': record['critical_grid']}), flush=True)
        print(json.dumps({'completed': completed, 'profile_classes': len(classes),
                          'seconds': time.monotonic()-start}), flush=True)
    rows = sorted(classes.values(), key=lambda row: row['critical_grid']['pressure'])
    for row in rows[:args.refinements]:
        pressure = Pressure(args.n, args.P, 2, args.geometry, row['zero_blocks'], row['infinity_active'])
        point = row['critical_grid']
        optimum = minimize(lambda value: pressure(*value), [math.log(point['scale']), point['charge']],
                           method='Nelder-Mead', bounds=[(-math.log(args.n)-4, math.log(1.5)),
                                                        (-2*args.n, 2*args.n)],
                           options={'maxiter': 200, 'fatol': 1e-10})
        refined = {'pressure': float(optimum.fun), 'scale': math.exp(float(optimum.x[0])),
                   'charge': float(optimum.x[1]), 'success': bool(optimum.success)}
        row['refined'] = refined if refined['pressure'] <= point['pressure'] else dict(point)
    expected = len(endpoint_blocks(args.endpoint))*4**len(perturbations)
    assert completed == expected == sum(row['rows'] for row in rows)
    result = {'status': 'Complete finite F4 operator catalogue with exact profiles and numerical pressure screen; not a global obstruction.',
              'n': args.n, 'ramification': args.P, 'geometry': args.geometry, 'slow_degree': args.degree,
              'endpoint_policy': args.endpoint, 'vacuum_color': args.vacuum_color,
              'top_exponent': args.top_exponent, 'completed_candidates': completed,
              'complete': True, 'profile_class_count': len(rows), 'grid_size': len(grid),
              'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in rows),
              'best_grid_relative_change': min(row['relative_grid']['change'] for row in rows),
              'profile_classes': rows, 'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'profile_classes'}), flush=True)


if __name__ == '__main__':
    main()
