#!/usr/bin/env python3
"""Exact controls for the full diagonal-word rank-sum barrier."""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from random import Random


def binary_rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


@lru_cache(maxsize=None)
def subset_indices(width, subset):
    positions = [j for j in range(width) if subset >> j & 1]
    return tuple(sum(((a >> i) & 1) << j for i, j in enumerate(positions))
                 for a in range(1 << len(positions)))


@lru_cache(maxsize=None)
def square_zero_rank(width, polynomial):
    size = 1 << width
    support = [a for a in range(size) if polynomial >> a & 1]
    columns = [sum(1 << (a | b) for a in support if not a & b) for b in range(size)]
    return binary_rank(columns)


def word_row_rank(width, polynomial):
    result = 0
    for subset in range(1 << width):
        indices = subset_indices(width, subset)
        restricted = sum(((polynomial >> a) & 1) << i for i, a in enumerate(indices))
        result += (1 << (width - subset.bit_count())) * square_zero_rank(subset.bit_count(), restricted)
    return result


def complement_coefficients(width, polynomial):
    full = (1 << width) - 1
    return sum(((polynomial >> a) & 1) << (full ^ a) for a in range(1 << width))


def lower_bounds(width, support):
    d = 1 << width
    upper, lower = [], []
    for subset in range(d):
        upper.append(max((1 << (width - a.bit_count()) for a in support if a & ~subset == 0), default=0))
        lower.append(max((1 << a.bit_count() for a in support if subset & ~a == 0), default=0))
    r, c = sum(upper), sum(lower)
    assert all(upper[a] * lower[a] >= d for a in support)
    assert r * c >= d * sum(x * y for x, y in zip(upper, lower))
    assert r * c >= d * d * len(support)
    return r, c


def check_all_binary_words():
    total = 0
    by_width = []
    for width in range(1, 5):
        d = 1 << width
        ranks = [word_row_rank(width, polynomial) for polynomial in range(1 << d)]
        for polynomial in range(1, 1 << d):
            support = [a for a in range(d) if polynomial >> a & 1]
            r, c = ranks[polynomial], ranks[complement_coefficients(width, polynomial)]
            rr, cc = lower_bounds(width, support)
            assert r >= rr and c >= cc
            assert r * c >= d * d * len(support)
            total += 1
        by_width.append({"width": width, "nonzero_polynomials": (1 << d) - 1})
    return {"nonzero_polynomials": total, "by_width": by_width}


def full_tensor(width, polynomial):
    d = 1 << width
    support = [a for a in range(d) if polynomial >> a & 1]
    rows = [0] * (d * d)
    cols = [0] * (d * d)
    for x in range(d):
        for y in range(d):
            for u in range(d):
                for v in range(d):
                    entry = 0
                    for a in support:
                        value = 1
                        for j in range(width):
                            xx, yy, uu, vv = ((label >> j) & 1 for label in (x, y, u, v))
                            identity = int(xx == uu and yy == vv)
                            crossing = int(xx == vv and yy == uu)
                            value &= (identity ^ crossing) if a >> j & 1 else identity
                        entry ^= value
                    rows[x * d + y] |= entry << (u * d + v)
                    cols[x * d + u] |= entry << (y * d + v)
    return rows, cols


def check_direct_tensors():
    rng = Random(20261030061)
    count = 0
    for width in (1, 2, 3):
        d = 1 << width
        polynomials = range(1, 1 << d) if width <= 2 else [rng.randrange(1, 1 << d) for _ in range(16)]
        for polynomial in polynomials:
            rows, cols = full_tensor(width, polynomial)
            assert binary_rank(rows) == word_row_rank(width, polynomial)
            assert binary_rank(cols) == word_row_rank(width, complement_coefficients(width, polynomial))
            count += 1
    return count


def gf4_multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def gf4_rank(matrix):
    rows = [row[:] for row in matrix]
    height, width = len(rows), len(rows[0])
    pivot = 0
    for column in range(width):
        selected = next((i for i in range(pivot, height) if rows[i][column]), None)
        if selected is None:
            continue
        rows[pivot], rows[selected] = rows[selected], rows[pivot]
        inverse = gf4_multiply(rows[pivot][column], rows[pivot][column])
        rows[pivot] = [gf4_multiply(inverse, x) for x in rows[pivot]]
        for i in range(pivot + 1, height):
            scalar = rows[i][column]
            if scalar:
                rows[i] = [x ^ gf4_multiply(scalar, y) for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == height:
            break
    return pivot


def gf4_word_rank(width, coefficients):
    result = 0
    for subset in range(1 << width):
        indices = subset_indices(width, subset)
        values = [coefficients[a] for a in indices]
        size = len(values)
        columns = []
        for b in range(size):
            column = [0] * size
            for a, value in enumerate(values):
                if not a & b:
                    column[a | b] ^= value
            columns.append(column)
        result += (1 << (width - subset.bit_count())) * gf4_rank(columns)
    return result


def check_extension_field():
    rng = Random(20261030062)
    count = 0
    for width in range(1, 5):
        d = 1 << width
        for _ in range(20):
            coefficients = [rng.randrange(4) for _ in range(d)]
            if not any(coefficients):
                coefficients[0] = 1
            support = [a for a, x in enumerate(coefficients) if x]
            r = gf4_word_rank(width, coefficients)
            c = gf4_word_rank(width, [coefficients[(d - 1) ^ a] for a in range(d)])
            rr, cc = lower_bounds(width, support)
            assert r >= rr and c >= cc and r * c >= d * d * len(support)
            count += 1
    return count


def convolution_rank(kernel, width):
    d = 1 << width
    return binary_rank([sum(((kernel >> (x ^ y)) & 1) << y for y in range(d)) for x in range(d)])


def general_kernel_negative_control():
    width, d = 4, 16
    rows = [sum(((s & w).bit_count() % 2) << w for w in range(d)) for s in range(d)]
    r = sum(convolution_rank(row, width) for row in rows)
    c = r  # The displayed bilinear kernel is symmetric.
    q = binary_rank(rows)
    assert (r, c, q) == (30, 30, 4)
    assert r * c < d * d * q
    return {"width": width, "row_rank": r, "column_rank": c,
            "ordinary_kernel_rank": q, "rank_product": r * c,
            "false_unrestricted_lower_bound": d * d * q}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    binary = check_all_binary_words()
    direct = check_direct_tensors()
    extension = check_extension_field()
    negative = general_kernel_negative_control()
    result = {
        "status": "PASS", "scope": "full unsymmetrized diagonal I,N word span; characteristic two",
        "universal_claim": "r*c >= d^2*word_support; crossing rank-sum >= d^(3/2)",
        "new_factoring_exponent": False,
        "binary_exhaustion": binary, "direct_tensor_controls": direct,
        "extension_field": "F4 = F2[a]/(a^2+a+1)", "extension_field_controls": extension,
        "general_kernel_negative_control": negative,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
