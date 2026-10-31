#!/usr/bin/env python3
"""Finite monomial and dual-monomial prime-state filter searches.

Coefficient cancellations and the resulting profiles are exact. Pressure
values and bounded continuous refinements are exploratory floating point
calculations. A finite negative search does not prove a family obstruction.
"""
from __future__ import annotations
import argparse
from itertools import product
import json
from math import comb, exp, log, sqrt
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from prime_state_filter_profiles import beta_polynomials
from prime_state_filter_pressure import filter_pressure, gaussian_pressure
from polynomial_state_filter_synthesis import transform


def monomial_profiles(prime, n, ramification, exponents):
    beta = beta_polynomials(prime, n, ramification, exponents)
    return tuple(min(value) for value in beta), tuple(
        (1 - 2 * n) * ramification * digit + 2 * e
        for digit, e in enumerate(exponents))


def dual_polynomials(prime, n, ramification, exponents):
    """Inverse transform of beta_k=(-1)^k q^f_k, with q=-v."""
    alpha = []
    for digit in range(prime):
        m = prime - 1 - digit
        value = {}
        for k in range(m + 1):
            for h in range(m - k + 1):
                degree = exponents[k] + n * ramification * h
                coefficient = comb(m, k) * comb(m - k, h) * (-1)**h
                value[degree] = (value.get(degree, 0) + coefficient) % prime
        alpha.append({degree: coefficient for degree, coefficient in value.items() if coefficient})
    return tuple(alpha)


def dual_profiles(prime, n, ramification, exponents):
    alpha = dual_polynomials(prime, n, ramification, exponents)
    return tuple(2 * e for e in exponents), tuple(
        (1 - 2 * n) * ramification * digit - 2 * max(value)
        for digit, value in enumerate(alpha))


def digit_log(prime, energy):
    return float(logsumexp(-np.arange(prime) * energy))


def equality_partition(values):
    seen = {}
    return tuple(seen.setdefault(value, len(seen)) for value in values)


def run_search(family, prime, n, ramification, refine=60):
    start = time.monotonic()
    scale = n * ramification
    levels = sorted(set((-2 * scale, -scale, -2, -1, 0, 1, 2, scale, 2 * scale)))
    parameters = [(0,) + rest for rest in product(levels, repeat=prime - 1)]
    if family == 'dual':
        parameters = [tuple(value + k * scale for k, value in enumerate(row)) for row in parameters]
    profile = monomial_profiles if family == 'monomial' else dual_profiles
    pairs = [profile(prime, n, ramification, row) for row in parameters]
    zero = np.array([pair[0] for pair in pairs], dtype=float)
    infinity = np.array([pair[1] for pair in pairs], dtype=float)
    charges = sorted(set((-2 * n, -n, -4, -2, -1, 0, .5, 1, 1.5, 2, 4, n, 2 * n)))
    scales = log(2 * n) / n * np.array((.125, .25, .5, .75, 1, 1.5, 2, 3, 4))
    digits = np.arange(prime)
    best = np.full(len(parameters), np.inf)
    best_charge, best_scale = np.zeros(len(parameters)), np.zeros(len(parameters))
    for charge, evaluation in product(charges, scales):
        constant = (gaussian_pressure(prime, n, charge, evaluation)
                    - digit_log(prime, (2 * n - 1 + charge) * evaluation)
                    - digit_log(prime, (1 - 2 * n) * evaluation))
        values = (constant + logsumexp(-(zero / ramification - (1 - charge) * digits) * evaluation, axis=1)
                  + logsumexp(-infinity / ramification * evaluation, axis=1))
        changed = values < best
        best[changed] = values[changed]
        best_charge[changed], best_scale[changed] = charge, evaluation
    records = []
    lower, upper = -log(n) - 4, log(min(4, 4 * log(2 * n) / n))
    for index in np.argsort(best)[:min(refine, len(parameters))]:
        beta_orders, infinity_orders = pairs[index]
        def objective(value):
            return filter_pressure(prime, n, ramification, float(value[0]), exp(float(value[1])),
                                   beta_orders, infinity_orders)
        initial = [max(-2 * n, min(2 * n, float(best_charge[index]))),
                   max(lower, min(upper, log(best_scale[index])))]
        optimized = minimize(objective, initial, method='Nelder-Mead',
                             bounds=[(-2 * n, 2 * n), (lower, upper)],
                             options={'maxiter': 250, 'xatol': 1e-8, 'fatol': 1e-10})
        grid_value = filter_pressure(prime, n, ramification, float(best_charge[index]),
                                     float(best_scale[index]), beta_orders, infinity_orders)
        use_grid = grid_value < float(optimized.fun)
        entry = {'parameters': parameters[index], 'beta_orders': beta_orders,
                 'infinity_active': infinity_orders,
                 'charge': float(best_charge[index]) if use_grid else float(optimized.x[0]),
                 'scale': float(best_scale[index]) if use_grid else exp(float(optimized.x[1])),
                 'pressure': grid_value if use_grid else float(optimized.fun),
                 'point_source': 'grid' if use_grid else 'bounded refinement',
                 'grid_pressure': float(best[index]), 'optimizer_success': bool(optimized.success)}
        if family == 'dual':
            alpha = dual_polynomials(prime, n, ramification, parameters[index])
            beta = transform(alpha, n * ramification, prime)
            assert beta == tuple({e: (-1)**k % prime} for k, e in enumerate(parameters[index]))
            assert all(sum(value.values()) % prime == 1 for value in alpha)
            entry['alpha'] = alpha
        records.append(entry)
    records.sort(key=lambda value: value['pressure'])
    return {'status': 'Finite exploratory search; no global obstruction or improved bound asserted.',
            'family': family, 'prime': prime, 'n': n, 'ramification': ramification,
            'levels': levels, 'candidates': len(parameters),
            'equality_partitions': len({equality_partition(row) for row in parameters}) if family == 'monomial' else None,
            'charge_grid': charges, 'scale_grid': scales.tolist(), 'grid_pairs': len(charges) * len(scales),
            'refinement_bounds': {'charge': [-2 * n, 2 * n], 'log_scale': [lower, upper]},
            'refined_candidates': len(records), 'best': records[0], 'top': records[:10],
            'seconds': time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=('monomial', 'dual'), default='monomial')
    parser.add_argument('--prime', type=int, default=5)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--refine', type=int, default=60)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.refine < 1:
        raise ValueError('Refine at least one candidate.')
    beta_polynomials(args.prime, args.n, args.ramification, (0,) * args.prime)
    result = run_search(args.family, args.prime, args.n, args.ramification, args.refine)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key in ('family', 'prime', 'n', 'ramification', 'candidates', 'best', 'seconds')}))


if __name__ == '__main__':
    main()
