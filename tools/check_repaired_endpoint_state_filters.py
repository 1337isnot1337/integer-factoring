#!/usr/bin/env python3
"""Exact diagonal-Y repairs, fixed-mode-total endpoints, and signed profiles."""
from __future__ import annotations

import argparse
from itertools import permutations, product
import json
from math import factorial
from pathlib import Path

import coupled_identity_extreme_filters as extreme
import coupled_identity_mode_filters as coupled
from check_coupled_identity_mode_filters import direct_profiles, expand_modes
from repaired_endpoint_state_filters import profile_data, scalar_filters


def scalar_profile_controls():
    count = 0
    for prime, n, geometry in ((3, 2, 'initial'), (3, 9, 'initial'), (3, 9, 'terminal'),
                               (5, 1, 'initial'), (5, 5, 'terminal')):
        choices = ((1,) * prime, tuple(1 if i % 2 == 0 else prime - 1 for i in range(prime)),
                   tuple(i + 1 for i in range(prime - 1)) + (1,))
        for ramification, coefficients in product((1, 3), choices):
            exponents = tuple((i * i + i) // 2 - 1 for i in range(prime))
            alpha = scalar_filters(exponents, coefficients)
            for charge in (0, 2):
                scalar = profile_data(prime, n, ramification, charge, exponents, coefficients, geometry)
                if geometry == 'initial':
                    matrix = coupled.profile_data(prime, n, ramification, charge, 1, alpha,
                                                   strict_endpoint=False)
                else:
                    matrix = extreme.profile_data(prime, n, ramification, charge, 1, alpha,
                                                   geometry='terminal', strict_endpoint=False)
                for key in ('zero_active', 'infinity_active', 'zero_spectators', 'infinity_spectators'):
                    assert matrix[key] == scalar[key]
                count += 1
    return {'signed_scalar_against_literal_active_matrix_profiles': count}


def endpoint_and_full_cut_controls():
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel

    cases = []
    equalities = 0
    choices = [(3, n, 1, 1, (0,), scalar_filters(exponents, coefficients))
               for n in (1, 2)
               for coefficients in ((1, 2, 1), (1, 1, 2), (1, 2, 2))
               for exponents in ((0, 1, 1), (0, -1, 2))]
    choices.append((5, 1, 1, 1, (0,), scalar_filters((0, 1, 1, -1, 2), (1, 2, 3, 4, 1))))
    joint = coupled.monomial_filters(3, 2, (0, 0, 0, 0, 1, 0, 0, 0, 0))
    joint[(1, 1)] = {-1: 2}
    choices.append((3, 2, 1, 2, (0, 1), joint))
    for prime, n, ramification, modes, indices, alpha in choices:
        raw, clearing = coupled.polynomial_seed(prime, n, ramification, modes, alpha, indices,
                                                 strict_endpoint=False)
        repaired, second_clearing = coupled.polynomial_seed(
            prime, n, ramification, modes, alpha, indices,
            strict_endpoint=False, repair_endpoint=True)
        assert clearing == second_clearing
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, prime), prime)
        evaluate = lambda value: sum(coefficient * (-1) ** i
                                     for i, coefficient in enumerate(value)) % prime
        endpoints = {state: sum(value.values()) % prime for state, value in alpha.items()}
        for powers in product(range(prime), repeat=4 * n):
            state = tuple(powers[n + j] for j in indices)
            expected = evaluate(gaussian.get(powers, ()))
            assert evaluate(raw.get(powers, ())) == endpoints[state] * expected % prime
            assert evaluate(repaired.get(powers, ())) == expected
            equalities += 2
        raw_orders, precision = direct_profiles(prime, raw, n, ramification, 0)
        repaired_orders, second_precision = direct_profiles(prime, repaired, n, ramification, 0)
        assert raw_orders == repaired_orders
        data = coupled.profile_data(prime, n, ramification, 0, modes, alpha, strict_endpoint=False)
        for cut, orders, scalar_shift in (('zero', raw_orders[0], 2 * clearing),
                                          ('infinity', raw_orders[1], -2 * clearing)):
            spectators = expand_modes(prime, data[cut + '_spectators'])
            expected = tuple(sorted(a + b + scalar_shift for a in data[cut + '_active']
                                    for b in spectators))
            assert orders == expected
        cases.append({'prime': prime, 'n': n, 'modes': modes,
                      'endpoint_values': [endpoints[state] for state in sorted(endpoints)],
                      'matrix_order': prime ** (2 * n), 'clearing_order': clearing,
                      'maximum_exact_jet_precision': max(precision, second_precision)})
    return {'raw_label_factor_and_repaired_endpoint_equalities': equalities,
            'complete_both_cut_replays': cases}


def fixed_mode_total_controls():
    from check_truncated_symmetric_crossings import basis

    records = []
    for prime, n, copies, modes, coefficients in ((3, 1, 3, 1, (1, 1, 2)),
                                                (5, 1, 2, 1, (1, 2, 3, 4, 1)),
                                                (3, 2, 2, 2, (1, 1, 1, 1, 2, 1, 1, 1, 1))):
        states = tuple(product(range(prime), repeat=modes))
        alpha = {state: {0: coefficient} for state, coefficient in zip(states, coefficients)}
        raw, _ = coupled.polynomial_seed(prime, n, 1, modes, alpha,
                                          strict_endpoint=False)
        repaired, _ = coupled.polynomial_seed(prime, n, 1, modes, alpha,
                                               strict_endpoint=False, repair_endpoint=True)
        digits = basis(n, prime)
        lookup = {word: i for i, word in enumerate(digits)}
        dimension = len(digits)
        entries = []
        for tensor in (raw, repaired):
            values = {}
            for powers, polynomial in tensor.items():
                labels = tuple(lookup[powers[j * n:(j + 1) * n]] for j in range(4))
                value = sum(coefficient * (-1) ** i for i, coefficient in enumerate(polynomial))
                for digit in powers[2 * n:]:
                    value *= factorial(digit)
                values[labels] = value % prime
            entries.append(values)
        words = tuple(product(range(dimension), repeat=copies))
        totals = {word: tuple(sum(digits[letter][j] for letter in word) for j in range(n))
                  for word in words}
        pairs = [(x, v) for x in words for v in words if totals[x] == totals[v]]
        scalar = {}
        for word in words:
            value = 1
            for label in word:
                value *= alpha[digits[label][:modes]][0]
            scalar[word] = value % prime
        count = 0
        for x, v in pairs:
            for y, u in pairs:
                raw_value = repaired_value = 1
                for i in range(copies):
                    labels = (v[i], y[i], u[i], x[i])
                    raw_value = raw_value * entries[0].get(labels, 0) % prime
                    repaired_value = repaired_value * entries[1].get(labels, 0) % prime
                expected = int(x == v and y == u)
                assert raw_value == scalar[y] * expected % prime
                assert repaired_value == expected
                count += 1
        records.append({'prime': prime, 'n': n, 'copies': copies,
                        'raw_and_repaired_fixed_mode_total_entries': count})
    # A full letter type gives a scalar. Mode totals alone need not do so.
    assert (0 + 2) == (1 + 1)
    assert (1 * 2) % 3 != (1 * 1) % 3
    return {'fixed_mode_total_crossings': records,
            'mode_totals_do_not_determine_nonlinear_type_scalar': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {
        'status': 'Exact nonzero-endpoint diagonal repairs and signed profiles; no pressure gain claimed.',
        'scalar_profile_controls': scalar_profile_controls(),
        'original_polynomial_tensor_controls': endpoint_and_full_cut_controls(),
        'typed_normalization': fixed_mode_total_controls(),
    }
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
