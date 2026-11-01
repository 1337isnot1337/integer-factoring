#!/usr/bin/env python3
"""Exact quadratic-minor obstruction and independent small Cayley orders."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import time

from high_axis_cayley import crossblock
from native_matrix_smith_f4 import MatrixOracle
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import pack


def dense_binary_rank(rows):
    rows = [row[:] for row in rows]
    rank = 0
    for j in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank+1, len(rows)):
            if rows[i][j]:
                rows[i] = [a ^ b for a,b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    begun, n, degree = time.monotonic(), 256, 2520
    e2 = [sum(1 << col for col in range(n)
              if not col & ~row and (row ^ col).bit_count() == 2)
          for row in range(n)]
    pivots, selected_rows, selected_cols = {}, [], []
    for i, original in enumerate(e2):
        value = original
        while value:
            j = value.bit_length()-1
            if j not in pivots:
                pivots[j] = value
                selected_rows.append(i)
                selected_cols.append(j)
                break
            value ^= pivots[j]
    rank = len(pivots)
    assert rank == 120 == (2**8-2**4)//2
    minor = [[(e2[row] >> col) & 1 for col in selected_cols] for row in selected_rows]
    assert dense_binary_rank(minor) == rank
    weight_difference = sum(i.bit_count() for i in selected_rows)-sum(j.bit_count() for j in selected_cols)
    assert weight_difference == 2*rank
    tangent_checks = []
    for r in range(9):
        gap = Fraction(degree, r+1)-Fraction(degree,3)+Fraction(degree*(r-2),9)
        assert gap == Fraction(degree*(r-2)**2,9*(r+1))
        assert (gap == 0) == (r == 2)
        tangent_checks.append({'augmentation_degree': r, 'valuation': degree//(r+1),
                               'strict_tangent_gap': [gap.numerator,gap.denominator]})
    assert 3*rank > n  # The Jensen multiplicity beats m*rho at rho=1/3.
    prefix = rank*(degree//3)
    assert prefix == 100800
    # All-one endpoint in the group basis, independently by subset zeta sums.
    endpoint = []
    for t in range(512):
        endpoint.append(sum(1 for s in range(512) if not t & ~s) % 2)
    assert endpoint == [0]*511+[1]
    # Compare complete, untruncated fraction-free profiles with the native engine.
    records, count = [], 0
    for dimension in (3,4,5):
        for values in ([1]*(dimension-1), [1+i % 3 for i in range(dimension-1)]):
            matrix, small_degree = crossblock(dimension, values)
            packed = [[pack(value) for value in row] for row in matrix]
            exact, _ = WeightedSmithOracle(packed).profile([0]*len(matrix), [0]*len(matrix), 1)
            got = MatrixOracle(matrix, len(matrix), precision=128).profile()['orders']
            assert list(exact) == list(map(int, got))
            assert all(a+b == 2*small_degree for a,b in zip(exact, reversed(exact)))
            records.append({'dimension': dimension, 'other_values': values, 'orders': list(exact)})
            count += len(exact)
    names = ['high_axis_cayley.py', 'check_high_axis_cayley.py', 'weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'dimension': 9, 'kernel_degree': degree,
              'one_particle_matrix_order': n, 'quadratic_multiplication_rank': rank,
              'minor_rows': selected_rows, 'minor_columns': selected_cols,
              'independent_dense_minor_rank': rank, 'total_augmentation_degree_difference': weight_difference,
              'specialized_minor_exact_valuation': prefix,
              'generic_prefix_upper_bound': prefix, 'tangent_checks': tangent_checks,
              'all_one_group_endpoint_coefficients_checked': len(endpoint),
              'untruncated_small_orders': count, 'small_controls': records,
              'pressure_floor': '1/3, for every positive scale and common physical frequency',
              'seconds': time.monotonic()-begun,
              'source_sha256': {name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'scope': 'Exact coefficient/minor proof for the explicit symmetric monomial Cayley '
                       'kernel. Not an obstruction to arbitrary nine-axis kernels, physical codes, '
                       'or stronger coefficient-rank methods. Complete factoring score unchanged.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status', 'quadratic_multiplication_rank',
                     'generic_prefix_upper_bound', 'untruncated_small_orders', 'seconds')}))


if __name__ == '__main__':
    main()
