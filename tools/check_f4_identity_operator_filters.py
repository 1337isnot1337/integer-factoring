#!/usr/bin/env python3
"""Original physical replays of filtered F4 occupation-transfer kernels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from coupled_identity_operator_filters import _geometry, identity_kernel, states, zero_profiles as binary_zero
from check_coupled_identity_operator_filters import example_kernel
from check_f4_repaired_identity_filters import physical_profiles, expanded
from f4_identity_operator_filters import (
    MUL, normalized_kernel, full_zero_profiles, zero_profiles,
    infinity_profiles, profile_data, polynomial_seed, evaluate_one,
)


def colored_kernel(modes, selected):
    alphabet = states(2, modes)
    kernel = {(state, state): {-(sum(state)**2-2*sum(state)): 1+i%3}
              for i, state in enumerate(alphabet)}
    for occupation in range(1, modes):
        sector = sorted((state for state in alphabet if sum(state) == occupation),
                        key=lambda state: sum(j*digit for j, digit in zip(selected, state)))
        low, high = sector[0], sector[-1]
        kernel[high, low] = {0: 2, -1: 1}
        kernel[low, high] = {0: 3, -1: 3}
    return kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [(2, 1, 0, 2, 'initial'), (2, 3, 2, 2, 'terminal'),
             (3, 1, -1, 2, 'initial'), (3, 1, 1, 3, 'initial')]
    if args.extended:
        cases.append((4, 1, 2, 2, 'mixed'))
    records = []
    for n, P, J, modes, geometry in cases:
        start = time.monotonic()
        selected = _geometry(2, n, modes, geometry)['selected']
        kernel = colored_kernel(modes, selected)
        actual = zero_profiles(n, P, modes, kernel, geometry)
        expected = full_zero_profiles(n, P, modes, kernel, geometry)
        assert actual == expected
        raw, clearing, gaussian = polynomial_seed(n, P, modes, kernel, geometry, False)
        repaired, other, _ = polynomial_seed(n, P, modes, kernel, geometry, True)
        assert clearing == other
        support = set(raw) | set(repaired) | set(gaussian)
        assert all(evaluate_one(repaired.get(state, {})) == evaluate_one(gaussian.get(state, {}))
                   for state in support)
        data = profile_data(n, P, J, modes, kernel, geometry, actual)
        expected_physical = (
            tuple(sorted(a+b+2*clearing for a in data['zero_active'] for b in expanded(data['zero_spectators']))),
            tuple(sorted(a+b-2*clearing for a in data['infinity_active'] for b in expanded(data['infinity_spectators']))),
        )
        raw_profiles, raw_precision = physical_profiles(raw, n, P, J)
        repaired_profiles, repaired_precision = physical_profiles(repaired, n, P, J)
        assert raw_profiles == repaired_profiles == expected_physical
        conjugate = {index: {degree: MUL[value][value] for degree, value in polynomial.items()}
                     for index, polynomial in kernel.items()}
        assert zero_profiles(n, P, modes, conjugate, geometry) == actual
        assert infinity_profiles(n, P, modes, conjugate, geometry) == data['infinity_active']
        records.append({'n': n, 'ramification': P, 'charge': J, 'modes': modes, 'geometry': geometry,
                        'physical_matrix_order': 2**(2*n), 'endpoint_coefficients_checked': len(support),
                        'clearing': clearing, 'maximum_native_precision': max(raw_precision, repaired_precision),
                        'zero_blocks': actual, 'infinity_active': data['infinity_active'],
                        'seconds': time.monotonic()-start})
    binary = 0
    for n, modes, geometry in ((3, 2, 'initial'), (3, 3, 'initial'), (4, 2, 'mixed'), (4, 2, 'terminal')):
        kernel = example_kernel(2, modes, True, 'filtered')
        assert zero_profiles(n, 1, modes, kernel, geometry) == binary_zero(2, n, 1, modes, kernel, geometry, 'filtered')
        binary += 1
    invalid = identity_kernel(2, 2)
    invalid[(1, 0), (0, 1)] = {0: 2}
    try:
        normalized_kernel(2, invalid, (0, 1))
    except ValueError:
        pass
    else:
        raise AssertionError('A nonregular physical endpoint repair was accepted.')
    result = {'status': 'Exact F4 off-diagonal operator and physical endpoint controls; no pressure claim.',
              'original_physical_cases': records, 'binary_regressions': binary,
              'invalid_filtered_endpoint_rejected': True, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'raw_and_repaired_cut_profiles': 4*len(records),
                      'full_vs_reduced_profiles': len(records), 'binary_regressions': binary,
                      'endpoint_coefficients_checked': sum(row['endpoint_coefficients_checked'] for row in records)}))


if __name__ == '__main__':
    main()
