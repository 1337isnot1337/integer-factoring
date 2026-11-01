#!/usr/bin/env python3
"""Exact charged Pfaffian gauges over F2[x,y,n]/(n^2)."""

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from pfaffian_particle_hole import local_exterior_table, particle_hole_mask


# The x,y exponents in this check stay below16. Truncation removes n^2.
X, Y, N = 1 << 1, 1 << 16, 1 << 256
MASK = (1 << 512) - 1
E = 1 ^ N


def mul(first, second):
    return multiply(first, second) & MASK


def matrix_product(first, second):
    result = [[0] * len(second[0]) for _ in first]
    for i, row in enumerate(first):
        for k, scalar in enumerate(row):
            if scalar:
                for j, other in enumerate(second[k]):
                    result[i][j] ^= mul(scalar, other)
    return result


def kronecker(first, second):
    return [[mul(first[i // len(second)][j // len(second[0])],
                 second[i % len(second)][j % len(second[0])])
             for j in range(len(first[0]) * len(second[0]))]
            for i in range(len(first) * len(second))]


def minor(matrix, rows, columns):
    @lru_cache(None)
    def recurse(position, remaining):
        if position == len(rows):
            return 1
        result = 0
        for j in range(len(columns)):
            if remaining >> j & 1:
                result ^= mul(matrix[rows[position]][columns[j]],
                              recurse(position + 1, remaining ^ (1 << j)))
        return result
    return recurse(0, (1 << len(columns)) - 1)


def exterior(matrix, scalar_parity=False):
    result = [[0] * 16 for _ in range(16)]
    subsets = [[i for i in range(4) if mask >> i & 1] for mask in range(16)]
    for i, rows in enumerate(subsets):
        for j, columns in enumerate(subsets):
            if len(rows) != len(columns):
                continue
            value = minor(matrix, rows, columns)
            if scalar_parity and (len(rows) + 1) % 2:
                value = mul(value, E)
            result[i][j] = value
    return result


def pfaffian_matrix(p, q):
    table = local_exterior_table(p, q, 1)
    result = [[0] * 16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        result[particle_hole_mask(x, y)][particle_hole_mask(v, u)] = table[((x * 4 + y) * 4 + u) * 4 + v] & MASK
    return result


def grades(mask):
    return (1 + ((mask >> 2) & 1) - ((mask >> 1) & 1),
            1 + ((mask >> 3) & 1) - (mask & 1))


def gauge_shifts(matrix):
    shifts = Counter()
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            value ^= int(i == j)
            assert value in (0, N)
            if value:
                a, b = grades(i)
                c, d = grades(j)
                shifts[a - c, b - d] += 1
    return [[list(shift), count] for shift, count in sorted(shifts.items())]


def verify_gauges():
    identity2 = [[1, 0], [0, 1]]
    identity16 = [[int(i == j) for j in range(16)] for i in range(16)]
    g = [[1, 1 ^ X], [1, 1]]
    h = [[1, 1], [1 ^ Y, 1 ^ X ^ Y]]
    ge = [[1, 1 ^ mul(X, E)], [1, 1]]
    he = [[1, 1], [1 ^ mul(Y, E), 1 ^ mul(X ^ Y, E)]]
    right_g = [[1, N], [0, E]]
    left_h = [[1, 0], [N, E]]
    left_g = [[E, N], [0, 1]]
    assert matrix_product(g, right_g) == matrix_product(left_g, g) == ge
    assert matrix_product(left_h, h) == he
    first = exterior(kronecker(left_h, identity2), scalar_parity=True)
    second = exterior(kronecker(identity2, right_g))
    left = exterior(kronecker(left_h, left_g), scalar_parity=True)
    uncolored = pfaffian_matrix(1 ^ X, 1 ^ Y)
    colored = pfaffian_matrix(1 ^ mul(X, E), 1 ^ mul(Y, E))
    assert matrix_product(matrix_product(first, uncolored), second) == colored
    assert matrix_product(left, uncolored) == colored
    assert matrix_product(first, first) == matrix_product(second, second) == matrix_product(left, left) == identity16
    first_shifts, second_shifts, left_shifts = map(gauge_shifts, (first, second, left))
    assert {tuple(row[0]) for row in first_shifts} <= {(0, 0), (1, 1)}
    assert {tuple(row[0]) for row in second_shifts} <= {(0, 0), (1, -1)}
    assert {tuple(row[0]) for row in left_shifts} <= {(0, 0), (1, 1), (1, -1)}
    return {"two_sided_polynomial_entries_checked": 256,
            "left_gauge_polynomial_entries_checked": 256,
            "two_sided_left_shift_counts": first_shifts,
            "two_sided_right_shift_counts": second_shifts,
            "one_sided_shift_counts": left_shifts,
            "constant_invertible_gauges_checked": 3}


def augmentation_pairing_rank(q, a, b, charge):
    first = [mask for mask in range(1 << q) if mask.bit_count() >= a]
    second = [mask for mask in range(1 << q) if mask.bit_count() >= b]
    rows = []
    for s in first:
        row = 0
        for j, t in enumerate(second):
            if not s & t and charge & ~(s | t) == 0:
                row |= 1 << j
        rows.append(row)
    return rank_two(rows)


def verify_augmentation_pairings():
    checked = 0
    for q in range(7):
        charges = range(1 << q) if q <= 4 else (0, (1 << q) - 1, (1 << q) // 3)
        for charge, a, b in product(charges, range(q + 2), range(q + 2)):
            expected = sum(comb(q, j) for j in range(a, q - b + 1))
            assert augmentation_pairing_rank(q, a, b, charge) == expected
            checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "gauges": verify_gauges(),
              "augmentation_pairing_rank_checks": verify_augmentation_pairings()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
