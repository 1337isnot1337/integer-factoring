#!/usr/bin/env python3
"""Exact all-scale floor for every polynomial boundary-only deformation."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from math import comb
import json
from pathlib import Path

from check_coherent_hot_boundary import logical
from check_hot_boundary_core_barrier import integer_at, integer_hash
from check_smith_profile_transfer import smith_orders
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import flatten
from restrict_implicit_reflection import restricted_tensor

ROOT = Path(__file__).resolve().parents[1]


def exact_core():
    words = [0, 1, 2, 3, 4, 5, 8]
    heights = [0, 1, 3, 4, 5, 6, 7]
    raw = restricted_tensor(16, list(map(logical, words)))
    profiles = []
    for cut in (0, 1):
        matrix, rows, columns = flatten(raw, cut)
        regular = []
        for i, (a, b) in enumerate(rows):
            row = []
            for j, (c, d) in enumerate(columns):
                low, high = matrix[i][j]
                assert high == 0
                if cut == 0:
                    shift = 2 * (heights[a] + heights[b]) - 128
                    if shift >= 0:
                        value = low << shift
                    else:
                        assert not low & ((1 << -shift) - 1)
                        value = low >> -shift
                else:
                    maximum = 128 + 2 * (heights[b] + heights[d])
                    assert low.bit_length() - 1 <= maximum
                    value = sum(1 << (maximum - k) for k in range(low.bit_length()) if low >> k & 1)
                row.append(value)
            regular.append(row)
        degree = max(value.bit_length() - 1 for row in regular for value in row)
        independent = smith_orders(regular, degree, precision=256)
        native = MatrixOracle([[(value, 0) for value in row] for row in regular], 49).profile()['orders']
        assert independent == list(map(int, native))
        assert len(independent) == 49 and all(not value % 2 for value in independent)
        profiles.append([value // 2 for value in independent])
    assert profiles[0] == profiles[1]
    return profiles[0], heights


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    orders, heights = exact_core()
    Z, Q = Counter(orders), Counter(heights + [8])
    a, b = 5, 7
    G = Counter()
    for i, ci in Z.items():
        for j, cj in Q.items():
            if i:
                G[a * (i - 1) + b * j] += i * ci * cj
            if j:
                G[a * i + b * (j - 1) + b - a] -= 2 * j * ci * cj
    G = {degree: value for degree, value in G.items() if value}
    degree = max(G)
    transformed = [sum(value * comb(degree - i, k - i)
                       for i, value in G.items() if i <= k) for k in range(degree + 1)]
    signs = [(value > 0) - (value < 0) for value in transformed]
    nonzero = [value for value in signs if value]
    variations = sum(x != y for x, y in zip(nonzero, nonzero[1:]))
    assert variations == 2 and G[0] == 1
    denominator = 10**8
    left, right = 90000111, 90000112
    gl, gr = (integer_at(G, numerator, denominator) for numerator in (left, right))
    assert gl < 0 < gr
    # Exact lower bound throughout the sole interior-minimum box.
    z_poly, q_poly = {a * i: ci for i, ci in Z.items()}, {b * j: cj for j, cj in Q.items()}
    zn, qn = integer_at(z_poly, left, denominator), integer_at(q_poly, right, denominator)
    zd, qd = denominator ** max(z_poly), denominator ** max(q_poly)
    lhs, rhs = zn**b * qd ** (2 * a), qn ** (2 * a) * zd**b
    assert lhs > rhs > 0
    assert sum(Z.values()) ** b > sum(Q.values()) ** (2 * a)
    assert 2 * a * a > b * b
    names = ['tools/check_moving_boundary_core_floor.py', 'tools/check_smith_profile_transfer.py',
             'tools/check_hot_boundary_core_barrier.py', 'tools/restrict_implicit_reflection.py',
             'tools/native_matrix_smith_f4.py', 'tools/native_matrix_smith_f4.cpp',
             'tools/native_rectangular_smith_f4.cpp']
    report = {'scope': 'All-scale natural common-type pressure floor for arbitrary polynomial tensors retaining the seven-letter coordinate subtensor and the stated eight physical heights; no factoring lower bound.',
              'status': 'passed', 'rho_floor': [a, b], 'above_post_half_integer_gap': 2 * a * a - b * b,
              'core_orders_t': orders, 'core_heights': heights, 'candidate_heights': heights + [8],
              'independent_binary_and_native_profiles': 2, 'profile_order_comparisons': 98,
              'derivative_coefficients': G, 'derivative_degree': degree,
              'mobius_signs': signs, 'mobius_sign_variations': variations,
              'minimum_box': [[left, denominator], [right, denominator]],
              'derivative_sign_hashes': [integer_hash(gl), integer_hash(gr)],
              'minimum_box_positive_gap_bits': (lhs - rhs).bit_length(),
              'minimum_box_gap_sha256': integer_hash(lhs - rhs),
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in ('status', 'rho_floor', 'derivative_degree', 'mobius_sign_variations', 'profile_order_comparisons')}))


if __name__ == '__main__':
    main()
