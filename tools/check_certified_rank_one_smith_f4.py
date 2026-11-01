#!/usr/bin/env python3
"""Compare structured profiles with full matrices and independent Toeplitz ranks."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random

from check_two_copy_crossing_couplings import binary_toeplitz_rank
from f4_recentered_lattices import pack, unpack
from native_matrix_smith_f4 import MatrixOracle
from native_certified_rank_one_smith_f4 import certified_rank_one_profile as rank_one_profile
from weighted_smith_oracle_f4 import polynomial_multiply


def literal_matrix(diagonal, left, right):
    return [[unpack(polynomial_multiply(pack(left[i]), pack(right[j]))
                    ^ (pack(diagonal[i]) if i == j else 0))
             for j in range(len(diagonal))] for i in range(len(diagonal))]


def verify():
    rng = random.Random(31102713)
    cases = []; counts = Counter(); pivots = Counter()
    for n in (1, 2, 3, 4, 5, 8, 12, 16):
        for trial in range(24):
            diagonal = [(1 << rng.randrange(13), 0) for _ in range(n)]
            left = []
            for _ in range(n):
                a, b = rng.getrandbits(7), rng.getrandbits(7)
                left.append((a ^ (a << 1), b ^ (b << 1)))
            right = [(rng.getrandbits(7), rng.getrandbits(7)) for _ in range(n)]
            matrix = literal_matrix(diagonal, left, right)
            # At one, the complete matrix is the identity, proving generic rank n.
            for i, row in enumerate(matrix):
                for j, (a, b) in enumerate(row):
                    assert a.bit_count() % 2 == int(i == j) and b.bit_count() % 2 == 0
                    counts['literal_full_rank_endpoints'] += 1
            ordinary = tuple(map(round, MatrixOracle(matrix, n, precision=128).profile()['orders']))
            first = rank_one_profile(diagonal, left, right, cutoff=64, work_precision=256)
            second = rank_one_profile(diagonal, left, right, cutoff=64, work_precision=512)
            assert tuple(first['orders']) == tuple(second['orders']) == ordinary
            pivots.update(first['pivot_counts'])
            counts['full_matrix_profile_comparisons'] += 1
            counts['independent_work_precision_comparisons'] += 1
            if n <= 5 and trial < 4:
                for precision in range(1, max(ordinary)+2):
                    assert binary_toeplitz_rank(matrix, precision) == 2*sum(max(0, precision-a) for a in ordinary)
                    counts['independent_binary_toeplitz_ranks'] += 1
            cases.append({'dimension': n, 'trial': trial, 'orders': ordinary,
                          'structured_pivot_counts': first['pivot_counts']})
    assert all(pivots[k] for k in ('diagonal_update', 'offdiagonal_to_arrow', 'arrow_center',
                                  'arrow_leaf_diagonal', 'arrow_leaf_column', 'arrow_leaf_row'))
    for diagonal, left, right, cutoff in (([(0, 0)], [(0, 0)], [(0, 0)], 64),
                                         ([(1 << 70, 0)], [(0, 0)], [(0, 0)], 64)):
        try:
            rank_one_profile(diagonal, left, right, cutoff=cutoff, work_precision=256)
        except ArithmeticError:
            counts['unresolved_negative_controls'] += 1
        else:
            raise AssertionError('A singular or out-of-window profile was accepted.')
    root = Path(__file__).resolve().parent
    return {'scope': 'Finite independent complete-matrix and Toeplitz validation of structured elimination with certified absolute error exponents.',
            'checks': dict(counts), 'all_pivot_branch_counts': dict(pivots), 'cases': cases,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in
                              ('native_certified_rank_one_smith_f4.py', 'native_certified_rank_one_smith_f4.cpp',
                               'native_rectangular_smith_f4.cpp', Path(__file__).name)}}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
