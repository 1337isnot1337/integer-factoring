#!/usr/bin/env python3
"""Finite floating screens of exact correlated two-mode polynomial filters.

Both profiles and the crossing endpoint come from the exact coupled
core theorems. This driver only guides discovery: a negative pressure
must receive an independent rigorous inequality certificate before it
is used as a factoring bound. A finite no-gain screen is not a barrier.
"""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import exp, fsum, log
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import coupled_identity_extreme_filters as extreme
import coupled_identity_mode_filters as initial
from prime_state_filter_pressure import gaussian_pressure


def prepare_pressure(prime, n, ramification, geometry, exponents, blocks):
    if geometry == 'initial':
        selected, retained = (0, 1), (n - 2, n - 1)
    else:
        data = extreme.geometry_data(prime, n, 2, geometry)
        selected, retained = data['selected'], data['retained']
    states = np.asarray(tuple(product(range(prime), repeat=2)), dtype=float)
    infinity = states @ np.asarray([2 * j + 1 - 2 * n for j in selected])
    infinity += 2 * np.asarray(exponents, dtype=float) / ramification
    zero_base = np.asarray([value / ramification for block in blocks for value in block['orders']])
    charges = np.asarray([block['charge'] for block in blocks for value in block['orders']])
    digits = np.arange(prime)

    def evaluate(log_scale, charge, relative=False):
        scale = exp(log_scale)
        zero = zero_base + (1 - charge) * charges
        gaussian = [digits * (2 * j + 1 - charge) for j in selected]
        gaussian += [digits * (2 * j + 1 + charge) for j in retained]
        gaussian += [digits * (2 * j + 1 - 2 * n) for j in selected]
        values = [zero, infinity] + gaussian
        signs = [1, 1] + [-1] * len(gaussian)
        minima = [float(np.min(value)) for value in values]
        ground = fsum(-sign * minimum for sign, minimum in zip(signs, minima))
        thermal = fsum(sign * float(logsumexp(-(value - minimum) * scale))
                       for sign, value, minimum in zip(signs, values, minima))
        change = ground * scale + thermal
        return change if relative else gaussian_pressure(prime, n, charge, scale) + change

    return evaluate


def candidate_exponents(prime, n, ramification, custom_grid=None):
    if prime == 2:
        grid = (tuple(custom_grid) if custom_grid is not None else
                tuple(sorted({-n * ramification, -(n * ramification // 2),
                              -2 * ramification, -ramification, 0,
                              ramification, 2 * ramification,
                              n * ramification // 2, n * ramification})))
        candidates = {(0,) + values for values in product(grid, repeat=3)}
        candidates.add((0, 0, 0, 0))
    else:
        states = tuple(product(range(prime), repeat=2))
        candidates = {(0,) * len(states)}
        for a, b, c in product((-1, 0, 1), (-1, 0, 1), (-2, -1, 1, 2)):
            candidates.add(tuple(ramification * (a * r + b * s + c * r * s)
                                 for r, s in states))
        for chosen in states:
            for value in (-n * ramification, -2 * ramification, -ramification,
                          ramification, 2 * ramification, n * ramification):
                exponents = tuple(value * int(state == chosen) for state in states)
                candidates.add(tuple(x - exponents[0] for x in exponents))
        for value in (-2, -1, 1, 2):
            for function in (lambda r, s: r == s,
                             lambda r, s: r + s == prime - 1,
                             lambda r, s: (r + s) % prime,
                             lambda r, s: (r - s) ** 2):
                exponents = tuple(ramification * value * int(function(r, s)) for r, s in states)
                candidates.add(tuple(x - exponents[0] for x in exponents))
    identity = (0,) * (prime * prime)
    return (identity,) + tuple(sorted(candidates - {identity}))


def exact_blocks(prime, n, ramification, geometry, exponents, python_backend=False):
    alpha = initial.monomial_filters(prime, 2, exponents)
    if prime == 2 and not python_backend:
        from binary_coupled_filter_profiles import zero_profiles
        if geometry == 'initial':
            return zero_profiles(n, ramification, 2, alpha)
        return zero_profiles(n, ramification, 2, alpha,
                             tensor_builder=extreme.active_tensor, geometry=geometry)
    if geometry == 'initial':
        return initial.zero_profiles(prime, n, ramification, 2, alpha)
    return extreme.zero_profiles(prime, n, ramification, 2, alpha, geometry)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=2)
    parser.add_argument('--n', type=int, default=8)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--geometry', choices=('initial', 'terminal', 'mixed'), default='mixed')
    parser.add_argument('--exponent-grid', help='Optional comma-separated integer grid for binary exponents.')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--refine', type=int, default=20)
    parser.add_argument('--python-backend', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    p, n, ramification = args.prime, args.n, args.ramification
    grid = tuple(map(int, args.exponent_grid.split(','))) if args.exponent_grid else None
    candidates = candidate_exponents(p, n, ramification, grid)
    if args.limit:
        candidates = candidates[:args.limit]
    scales = sorted({value / n for value in (.1, .3, 1, 2, 3, 5, 8)} | {.125, .25, .5, 1.})
    charges = sorted({1 - 2 * n, -n, -3, -1, 0, 1, 2, 3, n / 2, n, 2 * n})
    started = monotonic()
    rows = []
    result = {
        'status': 'Finite floating screen of exact profiles; no numerical value is a certified witness.',
        'prime': p, 'n': n, 'ramification': ramification, 'geometry': args.geometry,
        'rho': '1/sqrt(2)', 'evaluation_scale': 't=P*log(s)',
        'candidate_count': len(candidates), 'scale_grid': scales, 'charge_grid': charges,
        'rows': rows,
    }
    best = float('inf')
    best_relative = 0.
    for index, exponents in enumerate(candidates):
        blocks = exact_blocks(p, n, ramification, args.geometry, exponents, args.python_backend)
        evaluate = prepare_pressure(p, n, ramification, args.geometry, exponents, blocks)
        critical = min((evaluate(log(scale), charge), scale, charge)
                       for scale, charge in product(scales, charges))
        relative = min((evaluate(log(scale), charge, relative=True), scale, charge)
                       for scale, charge in product(scales, charges))
        row = {'exponents': exponents, 'critical_grid': critical,
               'relative_grid': relative, 'zero_charge_blocks_at_J1': blocks}
        rows.append(row)
        if critical[0] < best - 1e-9:
            best = critical[0]
            print('BEST', index, exponents, critical, relative, flush=True)
        if relative[0] < best_relative - 1e-9:
            best_relative = relative[0]
            print('RELATIVE_GAIN', index, exponents, critical, relative, flush=True)
        if (index + 1) % 20 == 0:
            result['elapsed_seconds'] = monotonic() - started
            args.output.write_text(json.dumps(result, indent=2) + '\n')
            print('PROGRESS', index + 1, len(candidates), result['elapsed_seconds'], flush=True)
    for row in sorted(rows, key=lambda value: value['critical_grid'][0])[:args.refine]:
        evaluate = prepare_pressure(p, n, ramification, args.geometry,
                                    row['exponents'], row['zero_charge_blocks_at_J1'])
        _, scale, charge = row['critical_grid']
        optimized = minimize(lambda point: evaluate(*point), [log(scale), charge],
                             method='Nelder-Mead',
                             bounds=[(-log(n) - 5, log(4)), (1 - 2 * n, 2 * n)],
                             options={'maxiter': 250, 'xatol': 1e-8, 'fatol': 1e-10})
        row['refined'] = {'critical_pressure': float(optimized.fun),
                          'scale': exp(float(optimized.x[0])), 'charge': float(optimized.x[1])}
        print('REFINED', row['exponents'], row['refined'], flush=True)
    result['elapsed_seconds'] = monotonic() - started
    result['completed_count'] = len(rows)
    result['minimum_sampled_relative_log_cost'] = min(row['relative_grid'][0] for row in rows)
    result['minimum_refined_critical_pressure'] = min(
        (row['refined']['critical_pressure'] for row in rows if 'refined' in row), default=None)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('DONE', len(rows), result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
