#!/usr/bin/env python3
"""Independent controls for sequential low-rank profiles and probe propagation."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random

from check_two_copy_crossing_couplings import binary_toeplitz_rank
from f4_recentered_lattices import pack, unpack
from native_certified_low_rank_smith_f4 import certified_low_rank_profile
from native_matrix_smith_f4 import MatrixOracle
from weighted_smith_oracle_f4 import polynomial_multiply


def literal_matrix(diagonal, updates):
    n = len(diagonal)
    result = [[pack(diagonal[i]) if i == j else 0 for j in range(n)] for i in range(n)]
    for left, right in updates:
        for i in range(n):
            for j in range(n):
                result[i][j] ^= polynomial_multiply(pack(left[i]), pack(right[j]))
    return [[unpack(value) for value in row] for row in result]


def verify():
    rng = random.Random(311026218); checks = Counter(); cases = []
    pivots = Counter(); probes = Counter(); margins = []
    for n in (1, 2, 3, 4, 8, 12, 16):
        for rank in (1, 2, 3, 6):
            for trial in range(24):
                diagonal = [(1 << rng.randrange(13), 0) for _ in range(n)]
                updates = []
                for _ in range(rank):
                    left = []
                    for _ in range(n):
                        a, b = rng.getrandbits(7), rng.getrandbits(7)
                        left.append((a ^ (a << 1), b ^ (b << 1)))
                    right = [(rng.getrandbits(7), rng.getrandbits(7)) for _ in range(n)]
                    updates.append((left, right))
                matrix = literal_matrix(diagonal, updates)
                for i, row in enumerate(matrix):
                    for j, (a, b) in enumerate(row):
                        assert a.bit_count() % 2 == int(i == j) and b.bit_count() % 2 == 0
                        checks['literal_full_rank_endpoints'] += 1
                original = tuple(map(round, MatrixOracle(matrix, n, precision=128).profile()['orders']))
                for precision in (256, 512):
                    result = certified_low_rank_profile(diagonal, updates, cutoff=64, work_precision=precision)
                    assert tuple(result['orders']) == original
                    checks['full_matrix_profile_comparisons'] += 1
                    if precision == 256:
                        pivots.update(result['pivot_counts']); probes.update(result['probe_pivot_counts'])
                        margins.append(result['minimum_certified_pivot_error_margin'])
                if trial < 2 and n <= 4:
                    for precision in range(1, max(original)+2):
                        assert binary_toeplitz_rank(matrix, precision) == 2*sum(max(0, precision-a) for a in original)
                        checks['independent_binary_toeplitz_ranks'] += 1
                if trial < 3 and rank > 1:
                    reverse = certified_low_rank_profile(diagonal, list(reversed(updates)), cutoff=64, work_precision=256)
                    assert tuple(reverse['orders']) == original
                    checks['update_order_invariance_controls'] += 1
                cases.append({'dimension': n, 'number_of_updates': rank, 'trial': trial, 'orders': original})
    branches = ('diagonal_update', 'offdiagonal_to_arrow', 'arrow_center',
                'arrow_leaf_diagonal', 'arrow_leaf_column', 'arrow_leaf_row')
    assert all(probes[k] > 0 for k in branches), probes
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'native_certified_low_rank_smith_f4.py',
             'native_certified_low_rank_smith_f4.cpp', 'native_certified_rank_one_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp')
    return {'scope': 'Exact finite original-matrix and independent Toeplitz controls for tracked-error sequential updates; all six probe-propagation branches exercised.',
            'checks': dict(checks), 'all_pivot_counts': dict(pivots),
            'probe_pivot_counts': dict(probes), 'minimum_certified_pivot_error_margin': min(margins),
            'cases': cases, 'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names}}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
