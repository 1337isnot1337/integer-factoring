#!/usr/bin/env python3
"""Exact regular-concentration, endpoint, and both-cut joint-filter checks."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations, product
import json
from math import factorial
from pathlib import Path

from coupled_identity_mode_filters import (
    add, concentration_data, monomial_filters, multiply,
    polynomial_seed, profile_data, shift,
)


def determinant(matrix, prime):
    result = {}
    for permutation in permutations(range(len(matrix))):
        value = {0: 1}
        for i, j in enumerate(permutation):
            value = multiply(value, matrix[i][j], prime)
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(len(matrix)) for j in range(i + 1, len(matrix)))
        result = add(result, value, prime, (-1) ** inversions)
    return result


def concentration_controls():
    cases = []
    identities = 0
    for prime, modes, n in product((2, 3, 5), (1, 2, 3), range(3, 11)):
        data = concentration_data(prime, n, modes)
        forms, top, interpolation = data['all_forms'], data['top'], data['interpolation']
        f = data['active_x']
        # Construct (1-tau)^(r(r-1)/2).
        expected_determinant = {0: 1}
        for _ in range(modes * (modes - 1) // 2):
            expected_determinant = multiply(expected_determinant, {0: 1, 1: prime - 1}, prime)
        assert determinant(f, prime) == expected_determinant
        for i in range(n - modes):
            for j in range(modes):
                value = {}
                for h in range(modes):
                    value = add(value, f[h][j], prime, interpolation[i][h])
                assert value == forms[i][j]
                identities += 1
        w = [[{0: int(h == k)} if h == k else {} for k in range(modes)]
             for h in range(modes)]
        for h, t in enumerate(top):
            for k in range(modes):
                for i in range(n - modes):
                    w[h][k] = add(w[h][k], {t - i: 1}, prime,
                                   interpolation[i][h] * interpolation[i][k])
                assert w[h][k].get(0, 0) == int(h == k)
                assert all(power >= 0 for power in w[h][k])
        assert determinant(w, prime).get(0) == 1
        for h in range(modes):
            for j in range(modes):
                value = {}
                for k in range(modes):
                    value = add(value, multiply(w[h][k], f[k][j], prime), prime)
                assert value == data['active_v'][h][j]
                identities += 1
        # Every nonconstant entry of W and every weighted first shear has
        # positive order. The theorem's inverse formulas then use a unit
        # determinant and regular adjugates, including when n=0 or1 mod p.
        assert all(t - i >= 1 for t in top for i in range(n - modes))
        cases.append({'prime': prime, 'n': n, 'modes': modes,
                      'n_mod_prime': n % prime, 'active_zero_modes': 2 * modes})
    return {'exact_polynomial_form_identities': identities,
            'interpolation_and_unit_gram_controls': cases}


def expand_modes(prime, orders):
    values = [0]
    for order in orders:
        values = [value + digit * order for value in values for digit in range(prime)]
    return values


def direct_profiles(prime, kernel, n, ramification, charge):
    from check_all_field_pascal_crossings import ramify_negative, scale, weighted_orders
    from check_truncated_symmetric_crossings import basis

    digits = basis(n, prime)
    lookup = {word: i for i, word in enumerate(digits)}
    dimension = prime ** n
    tensor = {}
    for powers, coefficient in kernel.items():
        labels = tuple(lookup[powers[i * n:(i + 1) * n]] for i in range(4))
        factor = 1
        for digit in powers[2 * n:]:
            factor *= factorial(digit)
        tensor[labels] = ramify_negative(scale(coefficient, factor, prime), prime)
    p = ramification
    heights = (
        [p * charge] * n, [p * (2 * i + 1) for i in range(n)],
        [-p * charge] * n, [p * (2 * i - 4 * n + 1) for i in range(n)],
        [p * (2 * i - n + 2) for i in range(n)], [-p * n] * n,
        [p * (2 * i - n + 1) for i in range(n)], [p * (3 * n - 1)] * n,
    )
    profiles = []
    precision = 0
    for endpoint, axes in enumerate(((3, 1, 0, 2), (3, 2, 0, 1))):
        full = [[()] * (dimension ** 2) for _ in range(dimension ** 2)]
        for labels, value in tensor.items():
            i = labels[axes[0]] * dimension + labels[axes[1]]
            j = labels[axes[2]] * dimension + labels[axes[3]]
            full[i][j] = value
        rows, columns, sectors = [], [], []
        for a, b in product(range(dimension), repeat=2):
            sectors.append(sum(digits[a]) + (1 if endpoint else -1) * sum(digits[b]))
            rows.append(sum(heights[4 * endpoint + axes[k]][i] * digits[label][i]
                            for k, label in enumerate((a, b)) for i in range(n)))
            columns.append(sum(heights[4 * endpoint + axes[k + 2]][i] * digits[label][i]
                               for k, label in enumerate((a, b)) for i in range(n)))
        orders = []
        for sector in sorted(set(sectors)):
            indices = [i for i, value in enumerate(sectors) if value == sector]
            profile, used = weighted_orders(
                [[full[i][j] for j in indices] for i in indices],
                [rows[i] for i in indices], [columns[j] for j in indices],
                prime, infinity=bool(endpoint))
            orders.extend(profile)
            precision = max(precision, used)
        profiles.append(tuple(sorted(orders)))
    return profiles, precision


def full_tensor_controls():
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel

    controls = []
    endpoint_entries = 0
    choices = []
    for n in (2, 3):
        for exponents in ((0, 0, 0, 0), (0, 0, 0, 1), (0, 1, 1, 0), (-1, 0, 1, 2)):
            choices.append((2, n, 1, monomial_filters(2, 2, exponents)))
    for exponents in ((0,) * 9, (0, 0, 0, 0, 1, 0, 0, 0, 0),
                      (0, 1, 0, 1, 0, -1, 0, -1, 0)):
        choices.append((3, 2, 1, monomial_filters(3, 2, exponents)))
    nonmonomial = monomial_filters(3, 2, (0,) * 9)
    nonmonomial[(1, 1)] = {0: 1, 1: 1, 2: 2}
    choices.append((3, 2, 1, nonmonomial))
    choices.append((2, 3, 3, monomial_filters(2, 2, (0, 0, 0, 1))))
    for prime, n, ramification, alpha in choices:
        kernel, clearing = polynomial_seed(prime, n, ramification, 2, alpha)
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, prime), prime)
        evaluate = lambda value: sum(coefficient * (-1) ** i
                                     for i, coefficient in enumerate(value)) % prime
        for powers in product(range(prime), repeat=4 * n):
            assert evaluate(kernel.get(powers, ())) == evaluate(gaussian.get(powers, ()))
            endpoint_entries += 1
        for charge in (0, 2):
            actual, precision = direct_profiles(prime, kernel, n, ramification, charge)
            data = profile_data(prime, n, ramification, charge, 2, alpha)
            expected = []
            for cut, scalar in (('zero', 2 * clearing), ('infinity', -2 * clearing)):
                spectators = expand_modes(prime, data[cut + '_spectators'])
                expected.append(tuple(sorted(a + b + scalar for a in data[cut + '_active']
                                             for b in spectators)))
            assert actual == expected, (prime, n, ramification, charge, actual, expected)
            controls.append({'prime': prime, 'n': n, 'ramification': ramification,
                             'charge': charge, 'clearing_order': clearing,
                             'matrix_order': prime ** (2 * n),
                             'maximum_series_precision': precision})
    return {'physical_endpoint_coefficient_equalities': endpoint_entries,
            'complete_two_cut_controls': controls}


def single_mode_reduction_controls():
    from prime_state_filter_profiles import active_profiles

    count = 0
    for prime, n, ramification, charge in product((2, 3), (2, 4), (1, 3), (0, 2)):
        exponents = tuple(digit * (digit - 1) // 2 for digit in range(prime))
        data = profile_data(prime, n, ramification, charge, 1,
                            monomial_filters(prime, 1, exponents))
        expected = active_profiles(prime, n, ramification, charge, exponents)
        for key in ('zero_active', 'infinity_active', 'zero_spectators', 'infinity_spectators'):
            assert data[key] == expected[key]
        count += 1
    return {'independent_prime_scalar_profile_comparisons': count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {
        'status': 'Exact bounded cores and coupled polynomial crossings; no pressure gain claimed.',
        'simultaneous_regular_concentration': concentration_controls(),
        'single_mode_regression': single_mode_reduction_controls(),
        'literal_polynomial_tensors': full_tensor_controls(),
    }
    encoded = json.dumps(result, indent=2, default=str) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
