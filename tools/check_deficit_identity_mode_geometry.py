#!/usr/bin/env python3
"""Polynomial and low-rank certificates for fixed-deficit identity cores."""
from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import combinations, product
import json
from math import comb
from pathlib import Path

from coupled_identity_mode_filters import add, multiply, shift
from coupled_identity_weyl_profiles import determinant
from deficit_identity_mode_geometry import geometry_data, weyl_data, remove_dummy_factor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [(2, 8, 0, 1, 1), (2, 8, 1, 0, 2), (2, 8, 1, 1, 1),
             (2, 8, 2, 0, 2), (2, 8, 2, 1, 1), (2, 8, 3, 1, 1),
             (2, 16, 2, 2, 1), (3, 3, 1, 0, 1), (3, 9, 1, 1, 1),
             (3, 9, 2, 0, 2), (3, 9, 2, 1, 1), (5, 5, 1, 1, 1),
             (5, 25, 2, 1, 1)]
    records, direct_count, gram_count, inverse_count, minor_count = [], 0, 0, 0, 0
    for p, parent, deficit, initial, terminal in cases:
        data = geometry_data(p, parent, deficit, initial, terminal)
        n, size = data['n'], data['augmented_modes']
        forms = data['forms']+data['tail_forms']
        for i in range(parent):
            for k, j in enumerate(data['augmented']):
                direct = {}
                for h in range(min(i, j)+1):
                    for degree in range(h+1):
                        value = comb(i, h)*comb(j, h)*comb(h, degree)*(-1)**degree % p
                        direct[degree] = (direct.get(degree, 0)+value) % p
                assert {degree: value for degree, value in direct.items() if value} == forms[i][k]
                direct_count += 1
        gram, full = [], []
        unit_parent = {degree: 1 for degree in range(parent)}
        for i in range(size):
            row, full_row = [], []
            for j in range(size):
                value, expected_full = {}, {}
                for a in range(parent):
                    term = shift(multiply(forms[a][i], forms[a][j], p), -a-1, p)
                    expected_full = add(expected_full, term, p)
                    if a < n:
                        value = add(value, term, p)
                target = shift(unit_parent, data['augmented'][i]-parent, p) if i == j else {}
                assert expected_full == target
                row.append(value)
                full_row.append(expected_full)
                gram_count += 1
            gram.append(row)
            full.append(full_row)
        local = weyl_data(p, parent, deficit, initial, terminal)
        delta, numerator = local['denominator_unit'], local['inverse_parameter_numerator']
        c = {0: 1, 1: p-1}
        for i in range(size):
            for j in range(size):
                value = {}
                for k in range(size):
                    first = add(numerator[i][k], delta if i == k else {}, p, -1)
                    value = add(value, multiply(first, multiply(c, gram[k][j], p), p), p)
                expected = shift(delta, data['augmented'][i], p) if i == j else {}
                assert value == expected
                inverse_count += 1
        scalar_denominator = {0: 1, parent: p-1}
        correction = [[add(multiply(scalar_denominator, numerator[i][j], p),
                           delta if i == j else {}, p, -1) for j in range(size)] for i in range(size)]
        if deficit < size:
            subsets = tuple(combinations(range(size), deficit+1))
            for rows, columns in product(subsets, repeat=2):
                assert not determinant([[correction[i][j] for j in columns] for i in rows], p)
                minor_count += 1
        records.append({'prime': p, 'parent': parent, 'n': n, 'deficit': deficit,
                        'initial_modes': initial, 'terminal_modes': terminal,
                        'selected': data['selected'], 'dummy': data['dummy'], 'retained': data['retained'],
                        'physical_active_modes': 2*data['modes']+deficit,
                        'reduced_matrix_order': p**size, 'rank_correction_bound': deficit})
    divisions = 0
    for p, dummy_orders in ((2, (6,)), (2, (8, 10)), (3, (4,)), (3, (10, 12)), (5, (8, 10, 12))):
        original = {-2: (-5, 1), -1: (-3, -3, 7), 0: (0, 4), 1: (2,)}
        current = original
        for order in dummy_orders:
            augmented = defaultdict(list)
            for charge, values in current.items():
                for digit in range(p):
                    augmented[charge+digit].extend(value+digit*order for value in values)
            current = {charge: tuple(sorted(values)) for charge, values in augmented.items()}
        for order in reversed(dummy_orders):
            current = remove_dummy_factor(current, p, order, min(current), max(current))
            divisions += 1
        assert current == original
    result = {'status': 'Uniform fixed-deficit construction with exact polynomial geometry and low-rank checks; no pressure or improved-bound claim.',
              'geometries': records, 'direct_form_identities': direct_count,
              'parent_Gram_entry_identities': gram_count, 'contraction_inverse_entry_identities': inverse_count,
              'vanishing_rank_bound_minors': minor_count, 'exact_dummy_factor_divisions': divisions}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'geometries'}))


if __name__ == '__main__':
    main()
