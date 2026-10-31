#!/usr/bin/env python3
"""Exact rational matrix controls for cached core updates and safe cutoffs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_rectangular_defect_core_basis import data_for
from check_rectangular_identity_transports import IDENTITY, conjugation_spec
from check_rectangular_refined_repairs import single_shear
from check_rectangular_f4_repair import binary_divide, binary_mul, times_binary
from check_rectangular_defect_support import pgcd
from rectangular_defect_coordinates import coordinate_matrix
from rectangular_coordinate_transforms import transform_coordinate_matrix, normalize


def compare_exact(first, second, below=None):
    A, SA, MA = first
    B, SB, MB = second
    assert len(A) == len(B)
    da, db = int(MA['unit_denominator_hex'], 16), int(MB['unit_denominator_hex'], 16)
    common = pgcd(da, db)
    fa, fb = binary_divide(db, common), binary_divide(da, common)
    for i in range(len(A)):
        for j in range(len(A)):
            a, b = times_binary(A[i][j], fa), times_binary(B[i][j], fb)
            shift = min(SA[i][j], SB[i][j])
            difference = tuple((a[p] << (SA[i][j]-shift)) ^ (b[p] << (SB[i][j]-shift)) for p in (0, 1))
            if below is None:
                assert difference == (0, 0), (i, j)
            elif below > shift:
                mask = (1 << (below-shift))-1
                assert not ((difference[0] | difference[1]) & mask), (i, j, below)
    return len(A)**2


def truncate_fixture(compiled, width=12):
    matrix, shifts, metadata = compiled
    mask = (1 << width)-1
    observed = [[(value[0]&mask, value[1]&mask) for value in row] for row in matrix]
    copied_shifts = [row[:] for row in shifts]
    for i, row in enumerate(observed):
        for j, value in enumerate(row):
            observed[i][j], copied_shifts[i][j] = normalize(value, copied_shifts[i][j])
    meta = dict(metadata, exact_polynomial_matrix=False, profile_precision=width)
    # Deliberately leave the native common cutoff high. The updater must also
    # use min_nonzero_shift+width, as required for individually trimmed cells.
    return observed, copied_shifts, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = []
    post, _, _ = single_shear(1, 11, 2)
    pre, wx, wy = single_shear(3, 8, 3)
    V = [row[:] for row in IDENTITY]
    V[10][1], V[0][0] = 2, 2
    C, W = conjugation_spec(1, 2, V)
    for place in ('infinity', 'zero'):
        data = data_for(place, 4, full_core=True)
        for name, base_arguments, output_map, input_map, target_arguments in (
            ('common_flag', {}, post, pre,
             {'repair_spec': (post, wx, wy), 'input_spec': (pre, wx, wy)}),
            ('arbitrary_identity_transport', {'repair_kind': 'inverse'}, V, W,
             {'identity_spec': (C, W)})):
            base = coordinate_matrix(data, **base_arguments)
            transformed = transform_coordinate_matrix(base, data, output_map, input_map)
            native = coordinate_matrix(data, **target_arguments)
            assert base[2]['exact_polynomial_matrix'] and native[2]['exact_polynomial_matrix']
            count = compare_exact(transformed, native)
            shortened = transform_coordinate_matrix(truncate_fixture(base), data, output_map, input_map)
            assert not shortened[2]['exact_polynomial_matrix']
            cutoff = shortened[2]['known_below']
            truncated_count = compare_exact(shortened, transformed, below=cutoff)
            records.append({'place': place, 'M': 4, 'case': name, 'exact_rational_entries_checked': count,
                            'truncated_entries_checked': truncated_count, 'truncated_known_below': cutoff,
                            'transform_seconds': transformed[2]['coordinate_transform_seconds'],
                            'precision_certificate': shortened[2]['coordinate_transform_precision']})
            print(json.dumps({'place': place, 'case': name, 'exact_entries': count}), flush=True)
    result = {'status': 'exact cached coordinate updates agree with independently rebuilt CAR matrices',
              'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
