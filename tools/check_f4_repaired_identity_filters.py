#!/usr/bin/env python3
"""Original physical endpoint and profile checks for repaired F4 filters."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

from f4_repaired_identity_filters import (
    MUL, active_tensor, evaluate_one, laurent_orders, monomial_filters,
    polynomial_seed, profile_data, reduced_zero_profiles, zero_profiles,
)


def expanded(orders):
    result = [0]
    for order in orders:
        result += [value + order for value in result]
    return result


def physical_profiles(tensor, n, ramification, charge):
    """Build the original physical matrices without the active reduction."""
    digits = tuple(product(range(2), repeat=n))
    lookup = {state: i for i, state in enumerate(digits)}
    dimension = 2 ** n
    heights = (
        [ramification * charge] * n,
        [ramification * (2 * i + 1) for i in range(n)],
        [-ramification * charge] * n,
        [ramification * (2 * i - 4 * n + 1) for i in range(n)],
        [ramification * (2 * i - n + 2) for i in range(n)],
        [-ramification * n] * n,
        [ramification * (2 * i - n + 1) for i in range(n)],
        [ramification * (3 * n - 1)] * n,
    )
    profiles, precision = [], 0
    for place, axes in enumerate(((3, 1, 0, 2), (3, 2, 0, 1))):
        matrix = [[{} for _ in range(dimension ** 2)] for _ in range(dimension ** 2)]
        for powers, value in tensor.items():
            labels = [lookup[powers[j * n:(j + 1) * n]] for j in range(4)]
            row = labels[axes[0]] * dimension + labels[axes[1]]
            column = labels[axes[2]] * dimension + labels[axes[3]]
            weight = sum(heights[4 * place + j][i] * powers[j * n + i]
                         for j in range(4) for i in range(n))
            matrix[row][column] = {weight + (1 if place == 0 else -1) * 2 * degree: coefficient
                                   for degree, coefficient in value.items()}
        sectors = [sum(digits[a]) + (1 if place else -1) * sum(digits[b])
                   for a, b in product(range(dimension), repeat=2)]
        orders = []
        for sector in sorted(set(sectors)):
            indices = [i for i, value in enumerate(sectors) if value == sector]
            values, used = laurent_orders([[matrix[i][j] for j in indices] for i in indices])
            orders.extend(values)
            precision = max(precision, used)
        profiles.append(tuple(sorted(orders)))
    return tuple(profiles), precision


def binary_controls():
    import coupled_identity_mode_filters as initial
    import coupled_identity_extreme_filters as extreme
    from binary_coupled_filter_profiles import zero_profiles as binary_profiles

    count = 0
    for n, ramification, geometry in ((3, 1, 'initial'), (8, 1, 'initial'),
                                       (8, 3, 'terminal'), (8, 1, 'mixed')):
        for exponents in ((0, 0, 0, 0), (0, 1, -1, 2), (0, -n, n, 1)):
            alpha = monomial_filters(2, exponents, (1, 1, 1, 1))
            actual = zero_profiles(n, ramification, 2, alpha, geometry)
            builder = initial.active_tensor if geometry == 'initial' else extreme.active_tensor
            options = {} if geometry == 'initial' else {'geometry': geometry}
            expected = binary_profiles(n, ramification, 2, alpha,
                                       tensor_builder=builder, **options)
            assert [(row['charge'], row['orders']) for row in actual] == [
                (row['charge'], row['orders']) for row in expected]
            count += 1
    return {'complete_F2_profile_regressions': count}


def original_controls():
    from coupled_identity_operator_filters import _geometry

    records, endpoints = [], 0
    cases = ((2, 1, 0, 2, 'initial'), (2, 3, 2, 2, 'terminal'),
             (3, 1, -1, 2, 'initial'), (3, 1, 1, 3, 'initial'),
             (4, 1, 2, 2, 'mixed'), (4, 1, -1, 2, 'terminal'))
    for n, ramification, charge, modes, geometry in cases:
        states = tuple(product(range(2), repeat=modes))
        exponents = tuple(sum(state) ** 2 - 2 * sum(state) for state in states)
        coefficients = tuple(1 + (i % 3) for i in range(len(states)))
        alpha = monomial_filters(modes, exponents, coefficients)
        if n == 3 and modes == 2:
            alpha[states[0]] = {0: 2, 1: 1}
        raw, clearing, gaussian = polynomial_seed(n, ramification, modes, alpha, geometry, False)
        repaired, other, _ = polynomial_seed(n, ramification, modes, alpha, geometry, True)
        assert clearing == other
        selected = _geometry(2, n, modes, geometry)['selected']
        for powers in product(range(2), repeat=4 * n):
            endpoint = evaluate_one(gaussian.get(powers, {}))
            state = tuple(powers[n + j] for j in selected)
            assert evaluate_one(raw.get(powers, {})) == MUL[endpoint][evaluate_one(alpha[state])]
            assert evaluate_one(repaired.get(powers, {})) == endpoint
            endpoints += 2
        data = profile_data(n, ramification, charge, modes, alpha, geometry)
        expected = (
            tuple(sorted(a + b + 2 * clearing for a in data['zero_active']
                         for b in expanded(data['zero_spectators']))),
            tuple(sorted(a + b - 2 * clearing for a in data['infinity_active']
                         for b in expanded(data['infinity_spectators']))),
        )
        actual, precision = physical_profiles(repaired, n, ramification, charge)
        assert actual == expected
        raw_actual, raw_precision = physical_profiles(raw, n, ramification, charge)
        assert raw_actual == actual
        conjugate = {state: {degree: MUL[value][value] for degree, value in polynomial.items()}
                     for state, polynomial in alpha.items()}
        conjugate_blocks = zero_profiles(n, ramification, modes, conjugate, geometry)
        assert [(row['charge'], row['orders']) for row in conjugate_blocks] == [
            (row['charge'], row['orders']) for row in data['zero_charge_blocks_at_J1']]
        reduced_blocks = reduced_zero_profiles(n, ramification, modes, alpha, geometry)
        assert [(row['charge'], row['orders']) for row in reduced_blocks] == [
            (row['charge'], row['orders']) for row in data['zero_charge_blocks_at_J1']]
        records.append({'n': n, 'ramification': ramification, 'charge': charge,
                        'modes': modes, 'geometry': geometry, 'clearing': clearing,
                        'filter_coefficients': [sorted(alpha[state].items()) for state in states],
                        'endpoint_values': [evaluate_one(alpha[state]) for state in states],
                        'full_physical_matrix_order': 2 ** (2 * n),
                        'maximum_exact_jet_precision': max(precision, raw_precision)})
    return {'raw_and_repaired_endpoint_equalities': endpoints,
            'full_physical_both_cut_cases': records,
            'conjugate_field_profile_controls': len(records),
            'complete_Weyl_reduced_profile_controls': len(records)}


def native_field_controls():
    for omega in (2, 3):
        square = MUL[omega][omega]
        actual, precision = laurent_orders([[{0: 1}, {0: omega}],
                                            [{0: omega}, {0: square, 600: 1}]])
        assert actual == (0, 600) and precision >= 1024
    try:
        laurent_orders([[{0: 2, 8192: 3}]])
    except ArithmeticError:
        pass
    else:
        raise AssertionError('A complete source beyond the supported window was truncated.')
    return {'late_second_plane_cancellation_controls': 2,
            'unsupported_complete_source_rejected': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'Exact repaired F4 endpoints and both physical profiles; no pressure claim.',
              'binary_regressions': binary_controls(),
              'original_tensor_controls': original_controls(),
              'native_field_controls': native_field_controls()}
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
