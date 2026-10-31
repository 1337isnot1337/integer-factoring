#!/usr/bin/env python3
"""Exact prime-power terminal/mixed interpolation and crossing controls."""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import comb
from pathlib import Path

from coupled_identity_extreme_filters import geometry_data, profile_data, reverse_form
from coupled_identity_mode_filters import add, monomial_filters, multiply, polynomial_seed, shift
from check_coupled_identity_mode_filters import determinant, direct_profiles, expand_modes


def raw_form(index, mode, prime):
    result = {}
    for h in range(min(index, mode) + 1):
        for power in range(h + 1):
            coefficient = comb(mode, h) * comb(index, h) * comb(h, power) * (-1) ** power
            result[power] = (result.get(power, 0) + coefficient) % prime
    return {power: value for power, value in result.items() if value}


def interpolation_controls():
    reversals = form_checks = 0
    cases = []
    for prime, n in ((2, 4), (2, 8), (2, 16), (3, 3), (3, 9), (5, 5), (5, 25)):
        for geometry, modes in (('terminal', 1), ('terminal', 2), ('terminal', 3), ('mixed', 2)):
            if modes > n:
                continue
            data = geometry_data(prime, n, modes, geometry)
            selected, retained, discarded = data['selected'], data['retained'], data['discarded']
            forms, f = data['all_forms'], data['active_x']
            for i in range(n):
                for j, mode in enumerate(selected):
                    assert forms[i][j] == raw_form(i, mode, prime)
                    reversals += 1
            unit = data['gaussian_v_unit']
            if geometry == 'terminal':
                b = tuple(tuple(shift({0: coefficient}, i - retained[h], prime)
                                if coefficient else {} for h, coefficient in enumerate(row))
                          for i, row in zip(discarded, data['interpolation_constants']))
                denominator = {0: 1}
            else:
                b = data['interpolation_numerators']
                denominator = unit
            for k, i in enumerate(discarded):
                for h, retained_index in enumerate(retained):
                    if b[k][h]:
                        assert min(b[k][h]) >= max(0, i - retained_index)
                for j in range(modes):
                    value = {}
                    for h in range(modes):
                        value = add(value, multiply(b[k][h], f[h][j], prime), prime)
                    assert value == multiply(denominator, forms[i][j], prime)
                    form_checks += 1
            denominator_squared = multiply(denominator, denominator, prime)
            w_numerator = []
            for h, retained_index in enumerate(retained):
                row = []
                for j in range(modes):
                    value = denominator_squared if h == j else {}
                    for k, i in enumerate(discarded):
                        value = add(value, shift(multiply(b[k][h], b[k][j], prime),
                                                  retained_index - i, prime), prime)
                    assert all(power >= 0 for power in value)
                    if h >= j:
                        assert value.get(0, 0) == int(h == j)
                    row.append(value)
                w_numerator.append(row)
            assert determinant(w_numerator, prime).get(0) == 1
            for h, retained_index in enumerate(retained):
                for j in range(modes):
                    value = multiply(denominator, f[h][j], prime)
                    for k, i in enumerate(discarded):
                        value = add(value, shift(multiply(b[k][h], forms[i][j], prime),
                                                  retained_index - i, prime), prime)
                    assert value == data['active_v'][h][j]
                    form_checks += 1
            cases.append({'prime': prime, 'n': n, 'modes': modes, 'geometry': geometry,
                          'selected': selected, 'retained': retained,
                          'zero_unit_degree': max(unit)})
    return {'direct_Lucas_polynomial_identities': reversals,
            'exact_interpolation_and_active_form_identities': form_checks,
            'weighted_regular_unit_controls': cases}


def scalar_terminal_controls():
    from prime_state_filter_profiles import beta_polynomials

    count = 0
    for prime, n in ((2, 4), (2, 16), (3, 3), (3, 9), (5, 5)):
        for ramification, charge in product((1, 3), (0, 2)):
            exponents = tuple((j * j + j) // 2 - 1 for j in range(prime))
            alpha = monomial_filters(prime, 1, exponents)
            actual = profile_data(prime, n, ramification, charge, 1, alpha, 'terminal')
            beta = beta_polynomials(prime, n, ramification, exponents)
            a = ramification * (2 * n - 1 - charge)
            expected = tuple(sorted(min(beta[k]) + a * (r - k)
                                    for r, k in product(range(prime), repeat=2)))
            assert actual['zero_active'] == expected
            assert actual['infinity_active'] == tuple(sorted(-ramification * j + 2 * e
                                                              for j, e in enumerate(exponents)))
            count += 1
    return {'independent_prime_scalar_profile_controls': count}


def complete_binary_controls():
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel

    records = []
    endpoints = 0
    for geometry, n, exponents in (
            ('terminal', 2, (0, 0, 0, 1)),
            ('mixed', 2, (0, 1, 1, 0)),
            ('terminal', 4, (0, 0, 0, 1)),
            ('mixed', 4, (0, 0, 0, 1))):
        prime, ramification, charge = 2, 1, 0
        alpha = monomial_filters(prime, 2, exponents)
        data = geometry_data(prime, n, 2, geometry)
        tensor, clearing = polynomial_seed(prime, n, ramification, 2, alpha, data['selected'])
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, prime), prime)
        for powers in product(range(prime), repeat=4 * n):
            assert sum(tensor.get(powers, ())) % prime == sum(gaussian.get(powers, ())) % prime
            endpoints += 1
        full, precision = direct_profiles(prime, tensor, n, ramification, charge)
        active = profile_data(prime, n, ramification, charge, 2, alpha, geometry)
        expected = []
        for cut, scalar in (('zero', 2 * clearing), ('infinity', -2 * clearing)):
            spectators = expand_modes(prime, active[cut + '_spectators'])
            expected.append(tuple(sorted(a + b + scalar for a in active[cut + '_active']
                                         for b in spectators)))
        assert full == expected
        records.append({'geometry': geometry, 'n': n, 'exponents': exponents,
                        'clearing_order': clearing, 'full_matrix_order': prime ** (2 * n),
                        'maximum_exact_jet_precision': precision})
    return {'physical_endpoint_coefficient_equalities': endpoints,
            'complete_both_cut_replays': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {
        'status': 'Exact prime-power extreme-mode cores and polynomial endpoints; no pressure gain claimed.',
        'regular_interpolation': interpolation_controls(),
        'single_terminal_mode': scalar_terminal_controls(),
        'full_polynomial_tensor_replays': complete_binary_controls(),
    }
    encoded = json.dumps(result, indent=2, default=str) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
