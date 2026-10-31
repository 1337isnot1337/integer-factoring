#!/usr/bin/env python3
"""Complete prime-field controls for arbitrary surplus-interval identity cores."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from coupled_identity_mode_filters import add, multiply, shift
from coupled_identity_weyl_profiles import determinant
from surplus_identity_mode_geometry import form, geometry_data, selected_geometry_data, selected_weyl_data


def c_power(prime, exponent):
    return {j: comb(exponent, j)*(-1)**j % prime for j in range(exponent+1)
            if comb(exponent, j) % prime}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    terminal = ((2, 4, 1, 2), (2, 8, 2, 3), (2, 8, 3, 3), (2, 32, 2, 4),
                (3, 3, 1, 1), (3, 9, 2, 3), (3, 27, 1, 2), (3, 9, 3, 3),
                (5, 5, 1, 2), (5, 5, 2, 2), (5, 25, 2, 3), (5, 5, 3, 3))
    cases = [(p, N, d, geometry_data(p, N, d, r)['selected']) for p, N, d, r in terminal]
    cases += [(2, 8, 2, (2, 6)), (2, 16, 1, (5, 11)),
              (3, 9, 2, (2, 6, 10)), (5, 5, 1, (1, 3))]
    lucas_count, determinants, grams, contractions, weighted = 0, 0, 0, 0, 0
    records = []
    for p, N, d, selected in cases:
        data = selected_geometry_data(p, N, d, selected)
        n, R = data['n'], data['augmented_modes']
        factor = add({0: 1}, c_power(p, N), p)
        for h in range(d):
            for i in range(N):
                assert form(p, i, N+h) == form(p, i, h)
                lucas_count += 1
            for k in range(d):
                assert form(p, N+k, N+h) == multiply(factor, form(p, k, h), p)
                lucas_count += 1
        for k in range(d):
            for j in data['old_augmented']:
                assert form(p, N+k, j) == form(p, k, j)
                lucas_count += 1
        small = [[form(p, i, j) for j in range(d)] for i in range(d)]
        assert determinant(small, p) == c_power(p, d*(d-1)//2)
        dx, dv = data['determinant_x'], data['determinant_v']
        expected = shift(c_power(p, N*d+d*(d-1)//2), 0, p, dx[0])
        assert dx == expected
        assert dv == multiply(dx, c_power(p, len(data['old_augmented'])*(N-1)), p)
        determinants += 3
        gram = []
        for a, j in enumerate(data['augmented']):
            row = []
            for b, k in enumerate(data['augmented']):
                original, reduced = {}, {}
                for i in range(n):
                    original = add(original, shift(multiply(form(p, i, j), form(p, i, k), p), -i-1, p), p)
                for h, i in enumerate(data['retained']):
                    reduced = add(reduced, shift(multiply(
                        data['augmented_v'][h][a], data['augmented_x'][h][b], p), -i-1, p), p)
                assert original == reduced
                row.append(original)
                grams += 1
            gram.append(row)
        parameter = selected_weyl_data(p, N, d, selected)
        delta, numerator = parameter['denominator_unit'], parameter['inverse_parameter_numerator']
        correction = [[add(value, delta if i == j else {}, p, -1)
                       for j, value in enumerate(row)] for i, row in enumerate(numerator)]
        c = {0: 1, 1: p-1}
        for i in range(R):
            for j in range(R):
                value = {}
                for k in range(R):
                    value = add(value, multiply(correction[i][k], multiply(c, gram[k][j], p), p), p)
                expected = shift(delta, data['augmented'][i], p) if i == j else {}
                assert value == expected
                contractions += 1
                value = correction[i][j]
                if value:
                    assert min(value) >= data['augmented'][i]+min(data['retained'])+1
                    assert min(value)-data['augmented'][i]+data['augmented'][j] >= 1
                weighted += 2
        assert set(data['augmented']) == set(selected)|set(range(d))|set(range(N, n))
        assert R <= len(selected)+2*d
        records.append({'prime': p, 'parent': N, 'surplus': d, 'n': n,
                        'selected': selected, 'augmented': data['augmented'],
                        'unchanged_actual_identities': data['unchanged_identity_indices'],
                        'retained': data['retained'], 'Weyl_matrix_order': p**R,
                        'actual_zero_active_modes': 2*R})
    from above_power_terminal_operators import geometry_data as binary_geometry, weyl_data as binary_weyl
    binary = 0
    for N, r in ((2, 1), (2, 2), (4, 2), (8, 3), (32, 3)):
        current, previous = geometry_data(2, N, 1, r), binary_geometry(N+1, r)
        for key in ('selected', 'augmented', 'retained', 'augmented_x', 'augmented_v'):
            assert current[key] == previous[key]
        a, b = selected_weyl_data(2, N, 1, current['selected']), binary_weyl(N+1, r)
        assert a['denominator_unit'] == b['denominator_unit']
        assert a['inverse_parameter_numerator'] == b['inverse_parameter_numerator']
        binary += 1
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/surplus_identity_mode_geometry.py', 'tools/check_surplus_identity_mode_geometry.py']
    result = {'status': 'Exact arbitrary-selection surplus geometry, unit blocks, Gram identities, and weighted Weyl parameters; no new full-tensor or pressure-search claim.',
              'geometries': records, 'Lucas_block_identities': lucas_count,
              'unit_determinant_identities': determinants, 'full_vs_reduced_Gram_entries': grams,
              'inverse_contraction_identities': contractions, 'weighted_regularity_entries': weighted,
              'binary_one_surplus_regressions': binary,
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in sources}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('geometries', 'source_sha256')}))


if __name__ == '__main__':
    main()
