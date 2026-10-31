#!/usr/bin/env python3
"""Original off-diagonal crossings for arbitrary prime-power mode selections."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

from coupled_identity_operator_filters import (
    identity_kernel, states, zero_profiles as full_prime_profiles,
)
from coupled_identity_weyl_profiles import zero_profiles as prime_profiles, weyl_data
from selected_identity_mode_filters import geometry_data
from check_f4_identity_operator_filters import colored_kernel
from check_f4_repaired_identity_filters import expanded, physical_profiles
from f4_identity_operator_filters import polynomial_seed, profile_data, zero_profiles, full_zero_profiles


def kernel(prime, modes, selected):
    result = identity_kernel(prime, modes)
    alphabet = states(prime, modes)
    for state in alphabet:
        exponent = sum(state)**2-2*sum(state)
        result[state, state] = {-exponent: 1}
    for occupation in range(1, modes*(prime-1)):
        sector = sorted((state for state in alphabet if sum(state) == occupation),
                        key=lambda state: sum(j*digit for j, digit in zip(selected, state)))
        if len(sector) > 1:
            first, last = sector[0], sector[-1]
            result[first, last] = {0: 1, -1: prime-1}
            result[last, first] = {0: 1, -1: prime-1}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls = []
    for prime, n, P, selected in ((2, 4, 1, (1, 2)), (2, 8, 3, (1, 4, 6)),
                                   (2, 32, 1, (3, 17, 29)), (3, 9, 1, (2, 6)),
                                   (5, 5, 1, (2,))):
        modes, geometry = len(selected), 'selected:'+','.join(map(str, selected))
        data = geometry_data(prime, n, selected)
        assert data['retained'] == tuple(n-1-j for j in selected)
        value = kernel(prime, modes, selected)
        compact = prime_profiles(prime, n, P, modes, value, geometry)
        assert compact == full_prime_profiles(prime, n, P, modes, value, geometry)
        gaussian = prime_profiles(prime, n, P, modes, identity_kernel(prime, modes), geometry)
        expected = {}
        alphabet = states(prime, modes)
        for first, second in product(alphabet, repeat=2):
            charge = sum(first)-sum(second)
            order = 2*P*sum(j*a+(n-j)*b for j, a, b in zip(selected, first, second))
            expected.setdefault(charge, []).append(order)
        assert gaussian == tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                                 for charge, orders in sorted(expected.items()))
        contraction = weyl_data(prime, n, modes, geometry)
        numerator = contraction['inverse_parameter_numerator']
        assert all(not numerator[i][j] for i in range(modes) for j in range(modes) if i != j)
        controls.append({'prime': prime, 'n': n, 'ramification': P, 'selected': selected,
                         'full_active_orders': prime**(2*modes),
                         'reduced_matrix_order': prime**modes})
    physical_cases = [(4, 1, 0, (1, 2))]
    if args.extended:
        physical_cases.append((4, 3, -1, (0, 1, 3)))
    physical = []
    for n, P, J, selected in physical_cases:
        modes, geometry = len(selected), 'selected:'+','.join(map(str, selected))
        value = colored_kernel(modes, selected)
        blocks = zero_profiles(n, P, modes, value, geometry)
        assert blocks == full_zero_profiles(n, P, modes, value, geometry)
        tensor, clearing, gaussian = polynomial_seed(n, P, modes, value, geometry)
        profile = profile_data(n, P, J, modes, value, geometry, blocks)
        expected = (
            tuple(sorted(a+b+2*clearing for a in profile['zero_active'] for b in expanded(profile['zero_spectators']))),
            tuple(sorted(a+b-2*clearing for a in profile['infinity_active'] for b in expanded(profile['infinity_spectators']))),
        )
        actual, precision = physical_profiles(tensor, n, P, J)
        assert actual == expected
        physical.append({'n': n, 'ramification': P, 'charge': J, 'selected': selected,
                         'original_matrix_order': 2**(2*n),
                         'endpoint_coefficients_checked': len(set(tensor)|set(gaussian)),
                         'maximum_native_precision': precision})
    rejected = 0
    for prime, n, selected in ((2, 6, (1, 2)), (2, 4, (2, 1)), (3, 9, (2, 2)), (2, 4, (0, 4))):
        try:
            geometry_data(prime, n, selected)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('An unsupported interval or invalid selection was accepted.')
    result = {'status': 'Exact arbitrary-selection operator and original physical controls; no pressure claim.',
              'prime_active_controls': controls, 'original_physical_controls': physical,
              'invalid_geometry_rejections': rejected, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'full_active_vs_reduced_profiles': len(controls)+len(physical),
                      'original_cut_profiles': 2*len(physical),
                      'endpoint_coefficients_checked': sum(row['endpoint_coefficients_checked'] for row in physical),
                      'invalid_geometry_rejections': rejected}))


if __name__ == '__main__':
    main()
