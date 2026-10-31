#!/usr/bin/env python3
"""Exhaust a Laurent F4 resonance catalogue across selected resonant quadruples."""
from __future__ import annotations

import argparse
from itertools import combinations, product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from f4_repaired_identity_filters import MUL, add
from resonant_identity_operators import FIRST, SECOND, STATES, pressure, signature
from prime_state_filter_pressure import gaussian_pressure


def endpoints():
    return tuple(((a, b), (c, d)) for a, b, c, d in product(range(4), repeat=4)
                 if MUL[a][d] ^ MUL[b][c])


def block(endpoint, correction, first_power=0, second_power=-1):
    return tuple(tuple(add({0: endpoint[i][j]},
                           {first_power: correction[2*i+j], second_power: correction[2*i+j]})
                       for j in range(2)) for i in range(2))


def catalogue(n, P, first_power=0, second_power=-1):
    if first_power == second_power:
        raise ValueError('Use distinct perturbation powers.')
    classes, count = {}, 0
    for endpoint in endpoints():
        for correction in product(range(4), repeat=4):
            values = signature(n, P, block(endpoint, correction, first_power, second_power))
            count += 1
            if values not in classes:
                classes[values] = {'signature': values, 'count': 0,
                                   'endpoint': endpoint, 'correction': correction}
            classes[values]['count'] += 1
    assert len(endpoints()) == 180 and count == 46080
    return tuple(classes.values())


def selections(n, family):
    if family == 'split':
        return ((0, 1, n-2, n-1),)
    if family == 'terminal':
        return ((n-4, n-3, n-2, n-1),)
    return tuple(row for row in combinations(range(n), 4) if row[0]+row[3] == row[1]+row[2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--P', type=int, default=1)
    parser.add_argument('--first-power', type=int, default=0)
    parser.add_argument('--second-power', type=int, default=-1)
    parser.add_argument('--window-lower', type=int, choices=(-2, -1),
                        help='Exhaust every entry in a complete consecutive three-power Laurent window.')
    parser.add_argument('--selections', choices=('all', 'split', 'terminal'), default='all')
    parser.add_argument('--refinements', type=int, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.n < 4 or args.n & (args.n-1):
        raise ValueError('Use a power-of-two interval of length at least four.')
    start = time.monotonic()
    if args.window_lower is None:
        classes = catalogue(args.n, args.P, args.first_power, args.second_power)
        exact_kernels = 46080
    else:
        from resonant_laurent_catalogue import complete_catalogue
        classes = complete_catalogue(args.n, args.P, args.window_lower)
        exact_kernels = 11796480
    selected = selections(args.n, args.selections)
    print(json.dumps({'exact_kernels': exact_kernels, 'profile_signatures': len(classes),
                      'selected_quadruples': len(selected), 'seconds': time.monotonic()-start}), flush=True)
    occupations = np.array([sum(state) for state in STATES])
    state_array = np.array(STATES)
    grid = tuple(product(math.log(2*args.n)/args.n*np.array((.25, .5, .75, 1., 1.5, 2.)),
                         sorted(set((-2*args.n, -args.n, -8, -4, -2, -1, 0, 1, 2, 4, 8, args.n/2, args.n, 2*args.n)))))
    scales = np.array([point[0] for point in grid])
    charges = np.array([point[1] for point in grid])
    gaussian = np.array([gaussian_pressure(2, args.n, float(J), float(t)) for t, J in grid])
    state_charge = occupations[:, None]*(1-charges)
    chunks = []
    for offset in range(0, len(selected), 128):
        values = selected[offset:offset+128]
        weights = np.array(values) @ state_array.T
        baseline_zero = 2*(weights+(4-occupations)*args.n)
        baseline_infinity = 2*weights+(1-2*args.n)*occupations
        reference = logsumexp(-(baseline_zero[:, :, None]+state_charge)*scales, axis=1)
        reference += logsumexp(-baseline_infinity[:, :, None]*scales, axis=1)
        chunks.append((values, weights, baseline_infinity, reference))
    best = math.inf
    for index, record in enumerate(classes):
        vacuum, single_a, single_b, pair_a, pair_b, infinity_a, infinity_b = record['signature']
        extra = []
        for state in STATES:
            if sum(state) == 0:
                value = vacuum
            elif sum(state) == 1:
                value = single_a if state[0] or state[3] else single_b
            elif state == FIRST:
                value = pair_a
            elif state == SECOND:
                value = pair_b
            else:
                value = (4-sum(state))*args.n*args.P
            extra.append(value/args.P)
        extra = np.array(extra)
        record['critical_grid'] = {'pressure': math.inf}
        record['relative_grid'] = {'change': math.inf}
        for values, weights, baseline_infinity, reference in chunks:
            zero = 2*(weights+extra)
            infinity = baseline_infinity.astype(float)
            infinity[:, STATES.index(FIRST)] += infinity_a/args.P
            infinity[:, STATES.index(SECOND)] += infinity_b/args.P
            relative = logsumexp(-(zero[:, :, None]+state_charge)*scales, axis=1)
            relative += logsumexp(-infinity[:, :, None]*scales, axis=1)
            relative -= reference
            actual = relative+gaussian
            i, j = np.unravel_index(np.argmin(actual), actual.shape)
            if actual[i, j] < record['critical_grid']['pressure']:
                record['critical_grid'] = {'pressure': float(actual[i, j]), 'selected': values[i],
                                           'scale': float(scales[j]), 'charge': float(charges[j])}
            i, j = np.unravel_index(np.argmin(relative), relative.shape)
            if relative[i, j] < record['relative_grid']['change']:
                record['relative_grid'] = {'change': float(relative[i, j]), 'selected': values[i],
                                           'scale': float(scales[j]), 'charge': float(charges[j])}
        if record['critical_grid']['pressure'] < best-1e-10:
            best = record['critical_grid']['pressure']
            print(json.dumps({'completed_signatures': index+1, 'best_grid': record['critical_grid']}), flush=True)
    rows = sorted(classes, key=lambda row: row['critical_grid']['pressure'])
    for record in rows[:args.refinements]:
        point = record['critical_grid']
        fit = minimize(lambda x: pressure(args.n, args.P, point['selected'], record['signature'],
                                          [math.exp(float(x[0]))], [float(x[1])])[0],
                       [math.log(point['scale']), point['charge']], method='Nelder-Mead',
                       bounds=[(-math.log(args.n)-4, math.log(2)), (-2*args.n, 2*args.n)],
                       options={'maxiter': 250, 'fatol': 1e-10})
        record['refined'] = {'pressure': float(fit.fun), 'selected': point['selected'],
                             'scale': math.exp(float(fit.x[0])), 'charge': float(fit.x[1]),
                             'success': bool(fit.success)}
    result = {'status': 'Complete specified Laurent resonance catalogue with exact signatures and numerical pressure search; not a global obstruction.',
              'n': args.n, 'ramification': args.P, 'first_power': args.first_power,
              'second_power': args.second_power, 'selection_family': args.selections,
              'window_lower': args.window_lower,
              'selected_quadruples': selected, 'exact_kernels': exact_kernels,
              'kernel_selection_pairs': exact_kernels*len(selected), 'profile_signatures': len(rows),
              'grid_size': len(grid), 'complete': True,
              'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in rows),
              'best_relative_grid_change': min(row['relative_grid']['change'] for row in rows),
              'records': rows, 'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('selected_quadruples', 'records')}), flush=True)


if __name__ == '__main__':
    main()
