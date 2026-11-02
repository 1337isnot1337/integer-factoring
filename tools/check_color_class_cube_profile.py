#!/usr/bin/env python3
"""Specialized profiles of the one-hot color-class Gaussian cube tensor.

The physical alphabet has one vacuum letter and q-1 colors.  A vertex
coloring is a word in the q letters, and

    T_q = product_{c=1}^{q-1} G_{S_c},

where G is the homogeneous eight-mode Gaussian cube tensor and S_c is the
set of vertices carrying color c.  This is a homogeneous tensor with the
same antipodal endpoint as G, and q=2 is exactly G.

The script computes the block Toeplitz rank over F_256 of the specialized
cut matrix.  A specialization can only increase Smith orders, so the
resulting small-order counts are lower bounds for the generic tensor.
The q=2 output is checked against the known exact Gaussian profile.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import product
from math import comb
import argparse
import json
from pathlib import Path

import numpy as np

import cube_homogeneous_gaussian as cube

MODULUS = 0x11B


def field_multiply(a: int, b: int) -> int:
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a & 256:
            a ^= MODULUS
        b >>= 1
    return result


MULTIPLY = np.zeros((256, 256), dtype=np.uint8)
for _a in range(256):
    for _b in range(256):
        MULTIPLY[_a, _b] = field_multiply(_a, _b)
INVERSE = np.zeros(256, dtype=np.uint8)
for _a in range(1, 256):
    for _b in range(1, 256):
        if MULTIPLY[_a, _b] == 1:
            INVERSE[_a] = _b
            break


def field_power(value: int, exponent: int) -> int:
    result = 1
    for _ in range(exponent):
        result = field_multiply(result, value)
    return result


def evaluate(polynomial, y: int, z: int, length: int) -> np.ndarray:
    result = np.zeros(length, dtype=np.uint8)
    for a, b, c in polynomial:
        if a < length:
            result[a] ^= field_multiply(field_power(y, b), field_power(z, c))
    return result


KERNEL = cube.kernel()
GAUSSIAN = cube.gaussian_tensor(KERNEL)
LEFT, RIGHT = cube.axis_parts(0)


def color_matrix(q: int, y: int, z: int, length: int) -> list[np.ndarray]:
    columns = list(product(range(q), repeat=4))
    evaluated = [evaluate(GAUSSIAN[mask], y, z, length) for mask in range(256)]

    def color_mask(left, right, color):
        mask = 0
        for i, vertex in enumerate(LEFT):
            if left[i] == color:
                mask |= 1 << vertex
        for i, vertex in enumerate(RIGHT):
            if right[i] == color:
                mask |= 1 << vertex
        return mask

    def entry(left, right):
        result = np.zeros(length, dtype=np.uint8)
        result[0] = 1
        for color in range(1, q):
            factor = evaluated[color_mask(left, right, color)]
            updated = np.zeros(length, dtype=np.uint8)
            for i, value in enumerate(result):
                if value:
                    stop = min(length - i, len(factor))
                    updated[i:i + stop] ^= MULTIPLY[value, factor[:stop]]
            result = updated
        return result

    size = len(columns)
    matrix = [[entry(left, right) for right in columns] for left in columns]
    return [np.array([[matrix[i][j][k] for j in range(size)]
                      for i in range(size)], dtype=np.uint8)
            for k in range(length)]


def block_toeplitz(blocks: list[np.ndarray], length: int) -> np.ndarray:
    size = blocks[0].shape[0]
    total = size * length
    result = np.zeros((total, total), dtype=np.uint8)
    for row in range(length):
        for column in range(row + 1):
            result[row * size:(row + 1) * size,
                   column * size:(column + 1) * size] = blocks[row - column]
    return result


def rank_over_f256(matrix: np.ndarray) -> int:
    matrix = matrix.copy()
    size = matrix.shape[0]
    rank = 0
    for column in range(size):
        if rank >= size:
            break
        nonzero = np.nonzero(matrix[rank:, column])[0]
        if nonzero.size == 0:
            continue
        pivot = rank + int(nonzero[0])
        if pivot != rank:
            matrix[[rank, pivot]] = matrix[[pivot, rank]]
        value = matrix[rank, column]
        if value != 1:
            matrix[rank] = MULTIPLY[INVERSE[value], matrix[rank]]
        below = matrix[rank + 1:, column]
        nonzero_rows = np.nonzero(below)[0]
        if nonzero_rows.size:
            rows = rank + 1 + nonzero_rows
            matrix[rows] ^= MULTIPLY[below[nonzero_rows, None],
                                     matrix[rank][None, :]]
        rank += 1
    return rank


def specialized_profile(q: int, y: int, z: int, max_order: int):
    length = max_order + 1
    blocks = color_matrix(q, y, z, length)
    ranks = [0]
    for cut in range(1, length + 1):
        ranks.append(rank_over_f256(block_toeplitz(blocks, cut)))
    cumulative = {order: ranks[order + 1] - ranks[order]
                  for order in range(length)}
    profile = {}
    previous = 0
    for order in range(length):
        count = cumulative[order]
        profile[order] = count - previous
        previous = count
    return {order: count for order, count in profile.items() if count}, ranks


def expected_q2_profile():
    weights = (2, 4, 8, 10)
    orders = [sum(weights[i] for i in range(4) if mask >> i & 1)
              for mask in range(16)]
    return dict(Counter(orders))


def graded_comparison():
    from cube_graded_minor_budget import envelope

    records = {}
    for q in (3, 4):
        colors = q - 1
        energies = Counter({0: 1, 3: 4 * colors, 6: 6 * colors ** 2,
                            9: 4 * colors ** 3, 12: colors ** 4})
        result = envelope(energies)
        prefix = sum(count for count, order in result["blocks"] if order <= 6)
        records[str(q)] = {
            "alphabet": result["alphabet"],
            "letter_energies": sorted(energies.items()),
            "blocks": [[count, str(order)] for count, order in result["blocks"][:8]],
            "prefix_count_through_order6": prefix,
        }
    assert records["3"]["prefix_count_through_order6"] == 65
    assert records["4"]["prefix_count_through_order6"] == 175
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    q2, q2_ranks = specialized_profile(2, 2, 3, 24)
    assert q2 == expected_q2_profile()
    assert sum(q2.values()) == 16

    q3_expected = {0: 1, 2: 2, 4: 3, 6: 4, 8: 5, 10: 8, 12: 9, 14: 12,
                   16: 10, 18: 4, 20: 7, 22: 4, 24: 6, 28: 3, 32: 2, 36: 1}
    q3, q3_ranks = specialized_profile(3, 2, 3, 36)
    assert q3 == q3_expected
    assert sum(q3.values()) == 81
    assert sum(order * count for order, count in q3.items()) == 1224

    stability = []
    for y, z in ((3, 5), (5, 7), (7, 11)):
        low, _ = specialized_profile(3, y, z, 10)
        stability.append({str((y, z)): sorted(low.items())})
        assert low == {order: count for order, count in q3.items() if order <= 10}

    low_orders = {}
    for q in (4, 5, 6):
        profile, _ = specialized_profile(q, 2, 3, 6)
        assert profile == {0: 1, 2: q - 1, 4: comb(q, 2), 6: comb(q + 1, 3)}
        assert sum(profile.values()) == comb(q + 2, 3)
        low_orders[str(q)] = sorted(profile.items())

    comparison = graded_comparison()
    result = {
        "status": "passed",
        "q2_profile": sorted(q2.items()),
        "q2_ranks": q2_ranks,
        "q3_specialized_profile": sorted(q3.items()),
        "q3_ranks": q3_ranks,
        "q3_specialization_stability": stability,
        "low_order_profiles": low_orders,
        "graded_comparison": comparison,
        "source_sha256": {
            name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("check_color_class_cube_profile.py",
                         "cube_homogeneous_gaussian.py",
                         "cube_graded_minor_budget.py")
        },
        "scope": "One-hot color-class Gaussian cube tensor; specialized F256 "
                 "block Toeplitz ranks.  Specialization increases orders, so "
                 "the reported profiles are lower bounds for the generic tensor. "
                 "No factoring exponent changes.",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": "passed",
        "q2_profile": sorted(q2.items()),
        "q3_specialized_profile": sorted(q3.items()),
        "q4_prefix_through_order6": comb(6, 3),
        "q4_graded_prefix_budget": comparison["4"]["prefix_count_through_order6"],
        "q3_prefix_through_order6": comb(5, 3),
        "q3_graded_prefix_budget": comparison["3"]["prefix_count_through_order6"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
