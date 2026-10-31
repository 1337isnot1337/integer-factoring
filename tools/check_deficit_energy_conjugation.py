#!/usr/bin/env python3
"""Exact endpoint and original physical controls for monomial kernel conjugation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_selected_deficit_operators import example
from check_prime_state_filter_profiles import physical_cut_profiles, expand_spectators
from coupled_identity_operator_filters import normalized_kernel
from prime_deficit_identity_operators import polynomial_seed, profile_data
from search_qutrit_deficit_operators import rank_one_catalogue, energy_conjugate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    endpoints = 0
    kernels = rank_one_catalogue(26, 1, 27)
    for selected in ((7, 18), (0, 13)):
        for twist in (-1, 1):
            keys = set()
            for kernel in kernels:
                conjugate = energy_conjugate(kernel, selected, 1, twist)
                normalized_kernel(3, 2, conjugate, 'identity', selected)
                assert energy_conjugate(conjugate, selected, 1, -twist) == kernel
                keys.add(tuple((row, column, tuple(sorted(value.items())))
                               for (row, column), value in sorted(conjugate.items())))
                endpoints += 1
            assert len(keys) == 1129
    physical = []
    for twist, P, J in ((-1, 1, 0), (1, 3, -2)):
        selected = (0, 2)
        kernel = energy_conjugate(example(2, 2), selected, P, twist, prime=2)
        tensor, clearing, gaussian = polynomial_seed(2, 4, 1, 0, 2, P, kernel,
                                                     selected_indices=selected)
        dense = {powers: tuple(value.get(i, 0) for i in range(max(value)+1))
                 for powers, value in tensor.items()}
        actual = physical_cut_profiles(dense, 2, 3, P, J)
        data = profile_data(2, 4, 1, 0, 2, P, J, kernel, selected_indices=selected)
        expected = (
            tuple(sorted(a+b+2*clearing for a in data['zero_active']
                         for b in expand_spectators(data['zero_spectators'], 2))),
            tuple(sorted(a+b-2*clearing for a in data['infinity_active']
                         for b in expand_spectators(data['infinity_spectators'], 2))),
        )
        assert actual == expected
        physical.append({'selected': selected, 'energy_twist': twist, 'ramification': P,
                         'aggregate_charge': J, 'matrix_order': 64,
                         'common_clearing': clearing,
                         'physical_polynomial_degree': max(degree for value in tensor.values() for degree in value),
                         'endpoint_coefficients': len(set(tensor)|set(gaussian))})
    result = {'status': 'Exact energy-conjugated endpoint and complete original two-cut controls.',
              'qutrit_endpoint_and_inverse_conjugation_controls': endpoints,
              'distinct_kernels_per_selection_and_twist': 1129,
              'original_physical_controls': physical}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
