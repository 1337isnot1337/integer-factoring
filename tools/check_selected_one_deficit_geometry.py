#!/usr/bin/env python3
"""Exact arbitrary-selection and weighted eigenbasis controls below a prime power."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

from coupled_identity_mode_filters import add, multiply, shift
from selected_one_deficit_geometry import geometry_data, compact_parameter, adjugate_parameter, contraction_eigenbasis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    for prime, parent in ((2, 4), (2, 8), (3, 3), (3, 9), (5, 5)):
        for modes in range(1, min(3, parent-2)+1):
            cases.extend((prime, parent, selected) for selected in combinations(range(parent-1), modes))
    cases.extend((3, 27, selected) for selected in combinations(range(26), 2))
    cases.extend((3, 27, selected) for selected in ((0, 1, 17), (1, 7, 13), (2, 11, 25)))
    entries, interpolation, eigen_entries, endpoint_singular = 0, 0, 0, 0
    records = []
    for prime, parent, selected in cases:
        data = geometry_data(prime, parent, selected)
        local, original = compact_parameter(data), adjugate_parameter(data)
        size = data['augmented_modes']
        for i in range(size):
            for j in range(size):
                assert multiply(local['inverse_parameter_numerator'][i][j], original['denominator_unit'], prime) == multiply(
                    original['inverse_parameter_numerator'][i][j], local['denominator_unit'], prime)
                entries += 1
        for row, i in enumerate(data['discarded']):
            for column in range(size):
                value = {}
                for k in range(size):
                    value = add(value, multiply(data['interpolation_numerators'][row][k],
                                                data['augmented_x'][k][column], prime), prime)
                assert value == multiply(data['unit'], data['forms'][i][column], prime)
                interpolation += 1
        eigen = contraction_eigenbasis(data)
        basis = eigen['basis']
        e0, h = local['missing_polynomial_unit'], local['missing_order']
        diagonal = [multiply({0: 1, parent: prime-1}, add(e0, {parent-h: 1}, prime), prime)]+[e0]*(size-1)
        for i in range(size):
            for j in range(size):
                value = {}
                for k in range(size):
                    value = add(value, multiply(local['inverse_parameter_numerator'][i][k], basis[k][j], prime), prime)
                assert value == multiply(basis[i][j], diagonal[j], prime)
                eigen_entries += 1
        assert eigen['endpoint_determinant'] == size % prime
        endpoint_singular += eigen['endpoint_determinant'] == 0
        if parent == 27:
            records.append({'selected': selected, 'retained': data['retained'],
                            'initial_prefix_length': data['initial_prefix_length'],
                            'eigenbasis_endpoint_singular': not eigen['endpoint_determinant']})
    result = {'status': 'Uniform arbitrary-selection reduction and local eigenbasis; exact polynomial checks, no pressure claim.',
              'geometry_count': len(cases), 'compact_parameter_entry_equalities': entries,
              'interpolation_equalities': interpolation, 'weighted_eigenbasis_entry_equalities': eigen_entries,
              'endpoint_singular_local_eigenbases': endpoint_singular,
              'all_n26_two_mode_selections': 325, 'n26_records': records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'n26_records'}))


if __name__ == '__main__':
    main()
