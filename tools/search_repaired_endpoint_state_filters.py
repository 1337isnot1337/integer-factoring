#!/usr/bin/env python3
"""Finite exact-profile screens of signed scalar occupation filters.

Every nonzero prime-field coefficient vector is included, up to common
scaling. Exponent tables are a documented finite family. Pressure is
evaluated in floating point on a documented finite charge/scale grid.
"""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure
from repaired_endpoint_state_filters import beta_polynomials, profile_data


def candidate_exponents(prime, n, ramification):
    if prime == 3:
        grid = sorted({-2 * n * ramification, -n * ramification, -n * ramification // 2,
                       -3 * ramification, -2 * ramification, -ramification, -1, 0, 1,
                       ramification, 2 * ramification, 3 * ramification,
                       n * ramification // 2, n * ramification, 2 * n * ramification})
        return tuple((0,) + values for values in product(grid, repeat=2))
    candidates = {(0,) + values for values in product((-ramification, 0, ramification),
                                                      repeat=prime - 1)}
    for value in (-n * ramification, -n * ramification // 2, -2 * ramification, -1,
                  1, 2 * ramification, n * ramification // 2, n * ramification):
        for function in (lambda digit: digit,
                         lambda digit: digit * (digit - 1) // 2,
                         lambda digit: digit % 2,
                         lambda digit: int(digit > 0),
                         lambda digit: int(digit == prime - 1)):
            candidates.add(tuple(value * function(digit) for digit in range(prime)))
    return tuple(sorted(candidates))


def screen(prime, n, ramification, geometry):
    if prime not in (3, 5):
        raise ValueError('This finite exponent catalog is defined for p=3 or p=5.')
    profile_data(prime, n, ramification, 0, (0,) * prime, (1,) * prime, geometry)
    selected = 0 if geometry == 'initial' else n - 1
    retained = n - 1 - selected
    scales = sorted({value / n for value in (.1, .3, 1, 2, 3, 5, 8)}
                    | {.125, .25, .5, 1., 2., 4.})
    charges = sorted({1 - 2 * n, -n, -3, -1, 0, 1, 2, 3, n / 2, n, 2 * n})
    points = tuple(product(scales, charges))
    scale = np.asarray([point[0] for point in points])
    charge = np.asarray([point[1] for point in points])
    a, b = 2 * selected + 1 - charge, 2 * retained + 1 + charge
    kappa = 2 * selected + 1 - 2 * n
    digits = np.arange(prime)
    gaussian_zero = b[:, None] * digits
    gaussian_infinity = kappa * digits
    zero_ground = np.min(gaussian_zero, axis=1)
    zero_thermal = logsumexp(-scale[:, None] * (gaussian_zero - zero_ground[:, None]), axis=1)
    infinity_ground = float(np.min(gaussian_infinity))
    infinity_thermal = logsumexp(-scale[:, None] * (gaussian_infinity - infinity_ground), axis=1)
    background = np.asarray([gaussian_pressure(prime, n, float(j), float(t)) for t, j in points])
    exponents = candidate_exponents(prime, n, ramification)
    coefficients = tuple((1,) + values for values in product(range(1, prime), repeat=prime - 1))
    seen, rows, count = {}, [], 0
    start = monotonic()
    for powers in exponents:
        infinity = kappa * digits + 2 * np.asarray(powers) / ramification
        new_infinity_ground = float(np.min(infinity))
        new_infinity_thermal = logsumexp(-scale[:, None] * (infinity - new_infinity_ground), axis=1)
        for colors in coefficients:
            orders = tuple(2 * min(value) for value in beta_polynomials(
                prime, n, ramification, powers, colors))
            count += 1
            key = powers, orders
            if key in seen:
                rows[seen[key]]['coefficient_multiplicity'] += 1
                continue
            seen[key] = len(rows)
            zero = np.asarray(orders) / ramification - a[:, None] * digits
            new_zero_ground = np.min(zero, axis=1)
            new_zero_thermal = logsumexp(-scale[:, None] * (zero - new_zero_ground[:, None]), axis=1)
            relative = (scale * (-new_zero_ground - new_infinity_ground
                                 + zero_ground + infinity_ground)
                        + new_zero_thermal + new_infinity_thermal - zero_thermal - infinity_thermal)
            critical = background + relative
            i, j = int(np.argmin(critical)), int(np.argmin(relative))
            rows.append({'exponents': powers, 'coefficients': colors, 'beta_orders': orders,
                         'critical_grid': (float(critical[i]), *points[i]),
                         'relative_grid': (float(relative[j]), *points[j]),
                         'coefficient_multiplicity': 1})
    return {
        'status': 'Finite floating screen of exact signed profiles; no certified pressure witness.',
        'prime': prime, 'n': n, 'ramification': ramification, 'geometry': geometry,
        'rho': '1/sqrt(2)', 'evaluation_scale': 't=P*log(s)',
        'exponent_vectors': len(exponents), 'nonzero_coefficient_vectors_per_exponent': len(coefficients),
        'completed_coefficient_cases': count, 'distinct_profile_count': len(rows),
        'scale_grid': scales, 'charge_grid': charges, 'gain_detection_tolerance': 1e-9,
        'minimum_sampled_critical_pressure': min(row['critical_grid'][0] for row in rows),
        'minimum_sampled_relative_log_cost': min(row['relative_grid'][0] for row in rows),
        'negative_critical_pressure_found': any(row['critical_grid'][0] < -1e-9 for row in rows),
        'same_gauge_gain_below_tolerance': any(row['relative_grid'][0] < -1e-9 for row in rows),
        'elapsed_seconds': monotonic() - start, 'rows': rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=5)
    parser.add_argument('--n', type=int, default=25)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--geometry', choices=('initial', 'terminal'), default='initial')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = screen(args.prime, args.n, args.ramification, args.geometry)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
