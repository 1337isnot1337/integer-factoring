#!/usr/bin/env python3
"""Independent direct-tensor and Gaussian controls for compact frozen rows."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

from check_rectangular_f4_repair import ZERO, binary_mul
from frozen_core_row_profiles import compile_interval, combine_components, reduce_components
from native_matrix_smith_f4 import MatrixOracle
from rectangular_frozen_core_rows import row_components, combine_row


ROWS = ([1]+[0]*11,
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0],
        [2, 3, 1, 2, 3, 2, 1, 1, 3, 2, 3, 1])


def direct_profile(seed, place, lam, mu):
    bulk, M = seed['bulk_size'], seed['ambient']
    counts = [1+a.bit_count() for a in range(bulk)]
    heights = seed['logical_heights']
    if place == 'infinity':
        wx = [(M+mu)*c for c in counts]
        wy = [(-M+lam)*c for c in counts]
        wu = [(1-M+mu)*c+2*h for c, h in zip(counts, heights)]
        wv = [(1-M+lam)*c+2*h for c, h in zip(counts, heights)]
        rw = [wv[v]+wy[y] for y in range(bulk) for v in range(bulk)]
        cw = [wu[u]+wx[x] for x in range(bulk) for u in range(bulk)]
    else:
        wx = [(1-2*M+lam)*c+2*h for c, h in zip(counts, heights)]
        wy = [(1-lam)*c+2*h for c, h in zip(counts, heights)]
        wu = [-mu*c for c in counts]
        wv = [mu*c for c in counts]
        rw = [wx[x]+wy[y] for y in range(bulk) for x in range(bulk)]
        cw = [wv[v]+wu[u] for u in range(bulk) for v in range(bulk)]
    n = bulk*bulk
    matrix = [[ZERO]*n for _ in range(n)]
    for x, y, u, v in product(range(bulk), repeat=4):
        row, column = ((v+bulk*y, u+bulk*x) if place == 'infinity'
                       else (x+bulk*y, v+bulk*u))
        matrix[row][column] = seed['values'][x+bulk*(y+bulk*(u+bulk*v))]
    return MatrixOracle(matrix, n).profile(rw, cw, infinity=place == 'infinity')['orders']


def check_unit_restoration(raw, reduced):
    original, shifts, _ = raw
    normalized, normalized_shifts, metadata = reduced
    assert shifts == normalized_shifts
    common = int(metadata['common_unit_removed'], 16)
    rows = [int(p, 16) for p in metadata['row_units_removed']]
    columns = [int(p, 16) for p in metadata['column_units_removed']]
    assert all(p & 1 for p in [common]+rows+columns)
    checks = 0
    for i in range(len(original[0])):
        for j in range(len(original[0])):
            factor = binary_mul(common, binary_mul(rows[i], columns[j]))
            for source in range(12):
                restored = tuple(binary_mul(p, factor) for p in normalized[source][i][j])
                assert restored == original[source][i][j]
                checks += 1
    return checks


def direct_controls():
    records, unit_checks = [], 0
    for beta in (0, 1):
        _, raw = compile_interval(4, beta, normalize_units=False)
        reduced = [reduce_components(*value) for value in raw]
        unit_checks += sum(check_unit_restoration(a, b) for a, b in zip(raw, reduced))
        components = row_components(beta=beta)
        for g in ROWS:
            seed = combine_row(components, g)
            for side, place in enumerate(('infinity', 'zero')):
                component_matrices, shifts, metadata = reduced[side]
                oracle = MatrixOracle(combine_components(component_matrices, g), 16, shifts,
                                      precision=metadata['profile_precision'])
                for lam, mu in ((0, 0), (0.25, -0.5), (-0.75, 1.25)):
                    actual = oracle.profile([lam*c for c in metadata['row_charges']],
                                            [mu*c for c in metadata['column_charges']])['orders']
                    expected = direct_profile(seed, place, lam, mu)
                    assert actual == expected
                    records.append({'M': 4, 'beta': beta, 'row': g, 'place': place,
                                    'lambda': lam, 'mu': mu, 'orders_checked': len(actual),
                                    'maximum_error': 0})
    return records, unit_checks


def gaussian_controls(M):
    data, raw = compile_interval(M, 1, normalize_units=False)
    reduced = [reduce_components(*value) for value in raw]
    unit_checks = sum(check_unit_restoration(a, b) for a, b in zip(raw, reduced))
    records = []
    for side, (d, (components, shifts, metadata)) in enumerate(zip(data, reduced)):
        total = (list(range(5-2*M, 2*M-4, 2)) if side == 0
                 else list(range(1, 2*M-4, 2))*2)
        total, spectators = Counter(total), Counter(d['spectator_orders'])
        assert not spectators-total
        active = list((total-spectators).elements())
        expected = [4-4*M if side == 0 else 0]
        for order in active:
            expected += [value+order for value in expected]
        n = len(components[0])
        oracle = MatrixOracle(combine_components(components, ROWS[0]), n, shifts,
                              precision=metadata['profile_precision'])
        for lam in (0, M/2, -M):
            actual = oracle.profile([lam*c for c in metadata['row_charges']],
                                    [-lam*c for c in metadata['column_charges']])['orders']
            assert actual == sorted(expected)
        records.append({'M': M, 'place': d['place'], 'matrix_order': n,
                        'active_one_particle_factors': active, 'expected_orders': sorted(expected),
                        'balanced_charges_checked': [0, M/2, -M],
                        'unit_reduced_maximum_bits': metadata['unit_reduced_maximum_bits']})
    return records, unit_checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--large-M', type=int, default=8)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    direct, unit_checks = direct_controls()
    gaussian, large_unit_checks = gaussian_controls(args.large_M)
    result = {'status': 'exact compact frozen-row profile and simultaneous unit-normalization controls',
              'direct_original_tensor_controls': direct,
              'gaussian_spectator_controls': gaussian,
              'restored_component_entries_checked': unit_checks+large_unit_checks,
              'scope': 'finite controls of the proved construction and inherited regular coordinate bases; no strict pressure claim'}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
