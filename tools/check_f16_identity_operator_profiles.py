#!/usr/bin/env python3
"""Independent F16 minor, base-field, and original physical operator checks."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random

from coupled_identity_operator_filters import _geometry, states
from f16_identity_operator_profiles import (
    evaluate_one, infinity_profiles, laurent_orders, polynomial_seed,
    profile_data, zero_profiles,
)
from f16_diagonal_filter_profiles import encode, smith_orders, zero_profiles as minor_diagonal
from f4_identity_operator_filters import zero_profiles as f4_zero, infinity_profiles as f4_infinity
from check_f4_identity_operator_filters import colored_kernel
from check_f4_repaired_identity_filters import expanded


def genuine_kernel(modes, selected):
    alphabet = states(2, modes)
    kernel = {(state, state): {-(sum(state)**2-2*sum(state)): 1+(3*i)%15}
              for i, state in enumerate(alphabet)}
    for occupation in range(1, modes):
        sector = sorted((state for state in alphabet if sum(state) == occupation),
                        key=lambda state: sum(j*digit for j, digit in zip(selected, state)))
        low, high = sector[0], sector[-1]
        kernel[high, low] = {0: 4, -1: 7}
        kernel[low, high] = {0: 12, -1: 12}
    return kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(163033)
    minors = 0
    for size in (1, 2, 3):
        for case in range(10):
            matrix = []
            for i in range(size):
                row = []
                for j in range(size):
                    a, b = rng.randrange(16), rng.randrange(16)
                    value = {0: int(i == j) ^ a ^ b, 1: a, 2: b}
                    row.append({degree: coefficient for degree, coefficient in value.items() if coefficient})
                matrix.append(row)
            expected = smith_orders([[encode(value) for value in row] for row in matrix])
            shift = case % 3-1
            shifted = [[{degree+shift: value for degree, value in entry.items()} for entry in row] for row in matrix]
            actual, precision = laurent_orders(shifted)
            assert actual == tuple(value+shift for value in expected)
            minors += 1
    late, precision = laurent_orders([[{0: 1}, {0: 1}], [{0: 1}, {0: 1, 600: 4}]])
    assert late == (0, 600) and precision >= 1024
    base_field, diagonal = 0, 0
    for n, P, r, geometry in ((2, 1, 2, 'mixed'), (3, 3, 3, 'initial'),
                               (4, 1, 3, 'split:1:2'), (8, 1, 3, 'terminal'),
                               (32, 1, 3, 'split:1:2')):
        selected = _geometry(2, n, r, geometry)['selected']
        kernel = colored_kernel(r, selected)
        assert zero_profiles(n, P, r, kernel, geometry) == f4_zero(n, P, r, kernel, geometry)
        assert infinity_profiles(n, P, r, kernel, geometry) == f4_infinity(n, P, r, kernel, geometry)
        base_field += 2
        for case in range(4):
            alpha = {}
            for state in states(2, r):
                while True:
                    value = {0: rng.randrange(16), 1: rng.randrange(16), 3: rng.randrange(16)}
                    if evaluate_one(value):
                        break
                alpha[state] = {degree: coefficient for degree, coefficient in value.items() if coefficient}
            kernel = {(state, state): value for state, value in alpha.items()}
            assert zero_profiles(n, P, r, kernel, geometry) == minor_diagonal(n, P, r, alpha, geometry)
            diagonal += 1
    # The original physical matrix builder is field-independent. Only its
    # exact scalar-field Smith oracle is replaced for this independent replay.
    import check_f4_repaired_identity_filters as physical
    previous = physical.laurent_orders
    physical.laurent_orders = laurent_orders
    cases = [(2, 1, 0, 2, 'mixed'), (3, 1, 1, 3, 'split:1:2')]
    if args.extended:
        cases.append((3, 3, -1, 2, 'initial'))
    records = []
    try:
        for n, P, J, r, geometry in cases:
            selected = _geometry(2, n, r, geometry)['selected']
            kernel = genuine_kernel(r, selected)
            raw, clearing, gaussian = polynomial_seed(n, P, r, kernel, geometry, False)
            repaired, other, _ = polynomial_seed(n, P, r, kernel, geometry, True)
            assert clearing == other
            data = profile_data(n, P, J, r, kernel, geometry)
            expected = (
                tuple(sorted(a+b+2*clearing for a in data['zero_active'] for b in expanded(data['zero_spectators']))),
                tuple(sorted(a+b-2*clearing for a in data['infinity_active'] for b in expanded(data['infinity_spectators']))),
            )
            first, p1 = physical.physical_profiles(raw, n, P, J)
            second, p2 = physical.physical_profiles(repaired, n, P, J)
            assert first == second == expected
            support = set(raw)|set(repaired)|set(gaussian)
            assert all(evaluate_one(repaired.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
                       for powers in support)
            records.append({'n': n, 'ramification': P, 'charge': J, 'modes': r,
                            'geometry': geometry, 'original_matrix_order': 2**(2*n),
                            'endpoint_coefficients_checked': len(support),
                            'maximum_native_precision': max(p1, p2)})
    finally:
        physical.laurent_orders = previous
    result = {'status': 'Exact F16 restriction-of-scalars and original physical controls; no pressure claim.',
              'independent_all_minor_profiles': minors, 'late_extension_coefficient_degree': 600,
              'F4_two_cut_regressions': base_field, 'independent_diagonal_minor_profiles': diagonal,
              'original_physical_cases': records, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'independent_all_minor_profiles': minors,
                      'independent_diagonal_minor_profiles': diagonal,
                      'F4_two_cut_regressions': base_field,
                      'raw_and_repaired_original_cut_profiles': 4*len(records),
                      'endpoint_coefficients_checked': sum(row['endpoint_coefficients_checked'] for row in records)}))


if __name__ == '__main__':
    main()
