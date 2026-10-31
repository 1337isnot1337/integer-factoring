#!/usr/bin/env python3
"""Independent exhaustive-minor checks for the weighted Smith oracle."""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import random

from check_schur_downset_crossing import sparse_kernel
from check_schur_gaussian_reduction import polynomial_multiply
from weighted_smith_oracle import WeightedSmithOracle


def enumerate_minors(matrix):
    n = len(matrix)
    @lru_cache(None)
    def determinant(rows, columns):
        if not rows:
            return 1
        row_bit = rows & -rows
        i = row_bit.bit_length() - 1
        result = 0
        remaining = columns
        while remaining:
            column_bit = remaining & -remaining
            remaining ^= column_bit
            j = column_bit.bit_length() - 1
            result ^= polynomial_multiply(matrix[i][j], determinant(rows ^ row_bit, columns ^ column_bit))
        return result
    result = []
    for k in range(1, n + 1):
        subsets = [mask for mask in range(1 << n) if mask.bit_count() == k]
        for rows in subsets:
            for columns in subsets:
                value = determinant(rows, columns)
                if value:
                    result.append((k, rows, columns, value))
    return result


def audit_matrix(matrix, rng, trials):
    n = len(matrix)
    all_minors = enumerate_minors(matrix)
    assert any(k == n for k, _, _, _ in all_minors)
    comparisons = 0
    for infinity in (False, True):
        oracle = WeightedSmithOracle(matrix, infinity)
        for trial in range(trials):
            denominator = 1 + trial % 3
            rows = [Fraction(rng.randrange(-13, 14), denominator) for _ in range(n)]
            columns = [Fraction(rng.randrange(-13, 14), denominator) for _ in range(n)]
            scale = Fraction(trial % 4, denominator)
            orders, prefixes = oracle.profile(rows, columns, scale)
            assert all(a <= b for a, b in zip(orders, orders[1:]))
            minima = {}
            polynomials = {}
            for k, row_mask, column_mask, value in all_minors:
                weight = scale * oracle.valuation(value)
                weight += sum(rows[i] for i in range(n) if row_mask >> i & 1)
                weight += sum(columns[j] for j in range(n) if column_mask >> j & 1)
                minima[k] = min(minima.get(k, weight), weight)
                polynomials[row_mask, column_mask] = value
            for k, prefix in enumerate(prefixes, 1):
                assert sum(orders[:k]) == minima[k]
                assert prefix.polynomial == polynomials[prefix.row_mask, prefix.column_mask]
                comparisons += 1
    return comparisons


def check():
    rng = random.Random(470172)
    comparisons = matrices = 0
    for n in range(1, 6):
        for sample in range(4):
            matrix = [[(rng.randrange(32) << 1) ^ int(i == j) for j in range(n)] for i in range(n)]
            comparisons += audit_matrix(matrix, rng, 12)
            matrices += 1
    kernel = sparse_kernel((0, 1, 2, 4), 8)
    for groups in (((1, 3), (2, 0)), ((1, 0), (2, 3))):
        for direction in (0, 1):
            rows = [4 * g + i for g in groups[direction] for i in range(4)]
            columns = [4 * g + i for g in groups[1 - direction] for i in range(4)]
            matrix = [[kernel[i][j] for j in columns] for i in rows]
            comparisons += audit_matrix(matrix, rng, 4)
            matrices += 1
    singular = WeightedSmithOracle([[1, 1], [1, 1]])
    try:
        singular.profile([0, 0], [0, 0])
    except ValueError as error:
        assert 'invertible' in str(error)
    else:
        raise AssertionError('singular negative control was not rejected')
    return {'status': 'pass', 'matrices': matrices,
            'independent_exhaustive_prefix_minima': comparisons,
            'weights': 'integer and rational row/column weights; zero and positive rational scales',
            'places': ['zero', 'infinity'], 'singular_negative_control': True,
            'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest()
                              for path in (Path(__file__), Path(__file__).with_name('weighted_smith_oracle.py'))},
            'scope': 'Exact weighted Smith oracle and determinantal separation; no improved factoring exponent.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rendered = json.dumps(check(), indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
