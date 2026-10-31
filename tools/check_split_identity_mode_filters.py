#!/usr/bin/env python3
"""Weighted interpolation, Gaussian recovery, and physical split-core checks."""
from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
import json
from pathlib import Path
import time

from coupled_identity_operator_filters import identity_kernel, states
from coupled_identity_weyl_profiles import zero_profiles as prime_profiles
from split_identity_mode_filters import geometry_data
from check_f4_identity_operator_filters import colored_kernel
from check_f4_repaired_identity_filters import physical_profiles, expanded
from f4_identity_operator_filters import polynomial_seed, profile_data, zero_profiles, full_zero_profiles


def gaussian_control(prime, n, a, b, ramification=1):
    data = geometry_data(prime, n, a, b)
    modes = a+b
    geometry = f'split:{a}:{b}'
    actual = prime_profiles(prime, n, ramification, modes, identity_kernel(prime, modes), geometry)
    expected = defaultdict(list)
    alphabet = states(prime, modes)
    for first, second in product(alphabet, repeat=2):
        order = 2*ramification*(sum(j*d for j, d in zip(data['selected'], first))
                               +sum((i+1)*d for i, d in zip(data['retained'], second)))
        expected[sum(first)-sum(second)].append(order)
    expected = tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                     for charge, orders in sorted(expected.items()))
    assert actual == expected
    margins = [margin for margin in data['weighted_order_margins'] if margin is not None]
    return {'prime': prime, 'n': n, 'initial_modes': a, 'terminal_modes': b,
            'ramification': ramification, 'selected': data['selected'], 'retained': data['retained'],
            'retained_determinant_constant': data['interpolation_denominator_unit'][0],
            'weighted_coefficients_checked': len(data['weighted_order_margins']),
            'minimum_nonzero_order_margin': min(margins, default=None),
            'complete_Gaussian_zero_orders': prime**(2*modes)}


def physical_control(n, ramification, charge, a, b):
    start = time.monotonic()
    modes, geometry = a+b, f'split:{a}:{b}'
    data = geometry_data(2, n, a, b)
    kernel = colored_kernel(modes, data['selected'])
    compact = zero_profiles(n, ramification, modes, kernel, geometry)
    assert compact == full_zero_profiles(n, ramification, modes, kernel, geometry)
    tensor, clearing, gaussian = polynomial_seed(n, ramification, modes, kernel, geometry)
    profile = profile_data(n, ramification, charge, modes, kernel, geometry, compact)
    expected = (
        tuple(sorted(a+b+2*clearing for a in profile['zero_active'] for b in expanded(profile['zero_spectators']))),
        tuple(sorted(a+b-2*clearing for a in profile['infinity_active'] for b in expanded(profile['infinity_spectators']))),
    )
    actual, precision = physical_profiles(tensor, n, ramification, charge)
    assert actual == expected
    return {'n': n, 'ramification': ramification, 'charge': charge,
            'initial_modes': a, 'terminal_modes': b,
            'original_matrix_order': 2**(2*n), 'active_matrix_order': 2**(2*modes),
            'reduced_matrix_order': 2**modes,
            'endpoint_coefficients_checked': len(set(tensor)|set(gaussian)),
            'maximum_native_precision': precision, 'seconds': time.monotonic()-start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [(2, 4, 1, 2), (2, 8, 1, 2), (2, 16, 1, 2), (2, 32, 1, 2),
             (2, 64, 1, 2), (2, 32, 2, 1), (2, 32, 2, 2),
             (3, 9, 1, 2), (3, 9, 2, 1), (5, 5, 1, 1)]
    gaussian = [gaussian_control(*case) for case in cases]
    physical_cases = [(3, 1, 0, 1, 2), (4, 1, 1, 1, 2)]
    if args.extended:
        physical_cases.append((4, 3, -1, 2, 1))
    physical = [physical_control(*case) for case in physical_cases]
    result = {'status': 'Exact mixed-end interpolation and original physical controls; no pressure claim.',
              'Gaussian_and_weighted_controls': gaussian,
              'original_physical_controls': physical, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'weighted_interpolation_entries': sum(row['weighted_coefficients_checked'] for row in gaussian),
                      'Gaussian_zero_orders': sum(row['complete_Gaussian_zero_orders'] for row in gaussian),
                      'original_cut_profiles': 2*len(physical),
                      'endpoint_coefficients_checked': sum(row['endpoint_coefficients_checked'] for row in physical)}))


if __name__ == '__main__':
    main()
