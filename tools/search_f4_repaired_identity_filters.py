#!/usr/bin/env python3
"""Finite pressure screens of exact repaired F4 two-direction filters.

The characteristic remains two; each binary occupation state receives
an arbitrary nonzero F4 coefficient. The diagonal endpoint repair has
zero profile cost. Floating pressure values only guide discovery.
"""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import exp, log
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from coupled_identity_operator_filters import _geometry
from f4_repaired_identity_filters import compiled_monomial_profiles, monomial_filters, zero_profiles
from prime_state_filter_pressure import gaussian_pressure
from search_coupled_identity_filters import candidate_exponents as pair_candidates


def candidate_exponents(n, ramification, modes, grid=None):
    if modes == 2:
        return pair_candidates(2, n, ramification, grid)
    if modes != 3 or grid is not None:
        raise ValueError('Use the pair grid or the fixed three-mode nonlinear exponent catalog.')
    states = tuple(product(range(2), repeat=modes))
    tables = {(0,) * len(states)}
    for chosen in states:
        for value in (-n * ramification, -n * ramification // 2,
                       -2 * ramification, -ramification, ramification,
                       2 * ramification, n * ramification // 2, n * ramification):
            values = tuple(value * int(state == chosen) for state in states)
            tables.add(tuple(value - values[0] for value in values))
    for a, b, c in product((-1, 0, 1), (-1, 0, 1), (-2, -1, 0, 1, 2)):
        tables.add(tuple(ramification * (a * sum(state) + b * sum(state) * (sum(state) - 1) // 2
                                         + c * int(sum(state) == 3)) for state in states))
    for pair_coefficients in product((-1, 0, 1), repeat=3):
        for cubic in (-1, 0, 1):
            tables.add(tuple(ramification * (
                pair_coefficients[0] * state[0] * state[1]
                + pair_coefficients[1] * state[0] * state[2]
                + pair_coefficients[2] * state[1] * state[2]
                + cubic * state[0] * state[1] * state[2]) for state in states))
    identity = (0,) * len(states)
    return (identity,) + tuple(sorted(tables - {identity}))


def pressure_grid(n, ramification, geometry, scales, charges, modes=2):
    data = _geometry(2, n, modes, geometry)
    selected, retained = data['selected'], data['retained']
    points = tuple(product(scales, charges))
    scale = np.asarray([point[0] for point in points])
    charge = np.asarray([point[1] for point in points])
    gaussian = np.asarray([2 * j + 1 - charge for j in selected]
                          + [2 * j + 1 + charge for j in retained]
                          + [np.full(len(points), 2 * j + 1 - 2 * n) for j in selected]).T
    gaussian_ground = np.minimum(gaussian, 0).sum(axis=1)
    gaussian_thermal = np.logaddexp(0, -scale[:, None] * np.abs(gaussian)).sum(axis=1)
    background = np.asarray([gaussian_pressure(2, n, float(j), float(t)) for t, j in points])
    states = np.asarray(tuple(product(range(2), repeat=modes)))
    infinity_base = states @ np.asarray([2 * j + 1 - 2 * n for j in selected])

    def evaluate(exponents, blocks):
        zero_base = np.asarray([value / ramification for block in blocks for value in block['orders']])
        block_charges = np.asarray([block['charge'] for block in blocks for value in block['orders']])
        zero = zero_base + (1 - charge[:, None]) * block_charges
        infinity = infinity_base + 2 * np.asarray(exponents) / ramification
        zero_ground = np.min(zero, axis=1)
        infinity_ground = np.min(infinity)
        zero_thermal = logsumexp(-scale[:, None] * (zero - zero_ground[:, None]), axis=1)
        infinity_thermal = logsumexp(-scale[:, None] * (infinity - infinity_ground), axis=1)
        relative = (scale * (-zero_ground - infinity_ground + gaussian_ground)
                    + zero_thermal + infinity_thermal - gaussian_thermal)
        critical = background + relative
        i, j = int(np.argmin(critical)), int(np.argmin(relative))
        return (float(critical[i]), *points[i]), (float(relative[j]), *points[j])

    return evaluate


def pressure_function(n, ramification, modes, geometry, exponents, blocks):
    data = _geometry(2, n, modes, geometry)
    selected, retained = data['selected'], data['retained']
    states = np.asarray(tuple(product(range(2), repeat=modes)))
    infinity = (states @ np.asarray([2 * j + 1 - 2 * n for j in selected])
                + 2 * np.asarray(exponents) / ramification)
    zero_base = np.asarray([value / ramification for block in blocks for value in block['orders']])
    block_charges = np.asarray([block['charge'] for block in blocks for value in block['orders']])

    def evaluate(log_scale, charge, relative=False):
        scale = exp(log_scale)
        zero = zero_base + (1 - charge) * block_charges
        gaussian = np.asarray([2 * j + 1 - charge for j in selected]
                              + [2 * j + 1 + charge for j in retained]
                              + [2 * j + 1 - 2 * n for j in selected])
        low_zero, low_infinity = np.min(zero), np.min(infinity)
        ground = -low_zero - low_infinity + np.minimum(gaussian, 0).sum()
        thermal = (logsumexp(-scale * (zero - low_zero))
                   + logsumexp(-scale * (infinity - low_infinity))
                   - np.logaddexp(0, -scale * np.abs(gaussian)).sum())
        change = float(scale * ground + thermal)
        return change if relative else gaussian_pressure(2, n, charge, scale) + change

    return evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=8)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--geometry', choices=('initial', 'terminal', 'mixed'), default='mixed')
    parser.add_argument('--exponent-grid')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--refine', type=int, default=16)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    grid = tuple(map(int, args.exponent_grid.split(','))) if args.exponent_grid else None
    candidates = candidate_exponents(args.n, args.ramification, args.modes, grid)
    if args.limit:
        candidates = candidates[:args.limit]
    colors = tuple((1,) + values for values in product((1, 2, 3), repeat=2 ** args.modes - 1))
    scales = sorted({value / args.n for value in (.1, .3, 1, 2, 3, 5, 8)} | {.125, .25, .5, 1., 2., 4.})
    charges = sorted({1 - 2 * args.n, -args.n, -3, -1, 0, 1, 2, 3,
                      args.n / 2, args.n, 2 * args.n})
    evaluate_grid = pressure_grid(args.n, args.ramification, args.geometry, scales, charges, args.modes)
    rows, seen = [], {}
    result = {'status': 'Finite floating screen of exact F4 profiles; no certified pressure witness.',
              'characteristic': 2, 'coefficient_field': 'F4', 'occupation_states': 2 ** args.modes,
              'modes': args.modes,
              'n': args.n, 'ramification': args.ramification, 'geometry': args.geometry,
              'rho': '1/sqrt(2)', 'evaluation_scale': 't=P*log(s)',
              'exponent_vectors': len(candidates), 'nonzero_coefficient_vectors_per_exponent': len(colors),
              'scale_grid': scales, 'charge_grid': charges, 'gain_detection_tolerance': 1e-9,
              'refinement_bounds': {'log_scale': [-log(args.n) - 5, log(4)],
                                    'charge': [1 - 2 * args.n, 2 * args.n]},
              'rows': rows}
    start, count, best, best_relative = monotonic(), 0, float('inf'), 0.
    for i, exponents in enumerate(candidates):
        compiled = (compiled_monomial_profiles(args.n, args.ramification, args.modes,
                                               exponents, args.geometry)
                    if args.modes > 2 else None)
        for coefficients in colors:
            if compiled is not None:
                blocks = compiled(coefficients)
            else:
                alpha = monomial_filters(args.modes, exponents, coefficients)
                blocks = zero_profiles(args.n, args.ramification, args.modes, alpha, args.geometry)
            count += 1
            key = exponents, tuple((block['charge'], block['orders']) for block in blocks)
            if key in seen:
                rows[seen[key]]['coefficient_multiplicity'] += 1
                continue
            seen[key] = len(rows)
            critical, relative = evaluate_grid(exponents, blocks)
            row = {'exponents': exponents, 'coefficients': coefficients,
                   'critical_grid': critical, 'relative_grid': relative,
                   'zero_charge_blocks_at_J1': blocks, 'coefficient_multiplicity': 1}
            rows.append(row)
            if critical[0] < best - 1e-9:
                best = critical[0]
                print('BEST', i, exponents, coefficients, critical, relative, flush=True)
            if relative[0] < best_relative - 1e-9:
                best_relative = relative[0]
                print('RELATIVE_GAIN', i, exponents, coefficients, critical, relative, flush=True)
        if (i + 1) % (5 if args.modes > 2 else 40) == 0:
            print('PROGRESS', i + 1, len(candidates), count, len(rows), monotonic() - start, flush=True)
    for row in sorted(rows, key=lambda value: value['critical_grid'][0])[:args.refine]:
        evaluate = pressure_function(args.n, args.ramification, args.modes, args.geometry,
                                     row['exponents'], row['zero_charge_blocks_at_J1'])
        value, scale, charge = row['critical_grid']
        assert abs(evaluate(log(scale), charge) - value) < 1e-9
        optimized = minimize(lambda point: evaluate(*point), [log(scale), charge],
                             method='Nelder-Mead',
                             bounds=[(-log(args.n) - 5, log(4)),
                                     (1 - 2 * args.n, 2 * args.n)],
                             options={'maxiter': 250, 'xatol': 1e-8, 'fatol': 1e-10})
        row['refined'] = {'critical_pressure': float(optimized.fun),
                          'scale': exp(float(optimized.x[0])), 'charge': float(optimized.x[1]),
                          'relative_log_cost': evaluate(*optimized.x, relative=True),
                          'scale_at_search_boundary': bool(
                              abs(optimized.x[0] - log(4)) < 1e-6
                              or abs(optimized.x[0] + log(args.n) + 5) < 1e-6)}
    result['elapsed_seconds'] = monotonic() - start
    result['completed_coefficient_cases'] = count
    result['distinct_profile_count'] = len(rows)
    result['minimum_sampled_relative_log_cost'] = min(row['relative_grid'][0] for row in rows)
    result['minimum_sampled_critical_pressure'] = min(row['critical_grid'][0] for row in rows)
    result['minimum_refined_critical_pressure'] = min(
        (row['refined']['critical_pressure'] for row in rows if 'refined' in row), default=None)
    result['same_gauge_gain_below_tolerance'] = any(
        row['relative_grid'][0] < -1e-9
        or row.get('refined', {}).get('relative_log_cost', 0) < -1e-9 for row in rows)
    result['negative_critical_pressure_found'] = any(
        row['critical_grid'][0] < -1e-9
        or row.get('refined', {}).get('critical_pressure', 0) < -1e-9 for row in rows)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('DONE', count, len(rows), result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
