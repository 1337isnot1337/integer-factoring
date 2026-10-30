#!/usr/bin/env python3
"""Complete scalar crossing profiles for small coded central groups.

The column cut is block diagonal by the horizontal abelianization invariant.
This computes actual piece ranks, beyond the one-cut bounds in the coded note.
"""

import argparse
import hashlib
import json
from pathlib import Path

from check_noncommutative_crossing_coded import (
    CodedGroup, constant_rank_polynomial, row_rank_upper_bounds,
)
from check_noncommutative_crossing_orbit import (
    binary_rank, group_word_polynomial, group_word_polynomial_independent,
)


def profile(columns, center, dense_check=False):
    group = CodedGroup(columns, center)
    d, v, z = group.d, group.v, 1 << center
    weights = (1,) * group.n + (2,) * center
    maximum = sum(weights)
    order = tuple(range(group.width))
    kernels = [group_word_polynomial(group, S, order, weights) for S in range(d)]
    if d <= 64:
        assert all(kernel == group_word_polynomial_independent(group, S, order, weights)
                   for S, kernel in enumerate(kernels))
    products = [[group.mul(a, b) for b in range(d)] for a in range(d)]
    row_ranks = [0] * (maximum + 1)
    by_horizontal = [[] for _ in range(v)]
    for S, kernel in enumerate(kernels):
        matrices = [[0] * d for _ in row_ranks]
        assert kernel.get(0, 0).bit_count() & 1
        for a, polynomial in kernel.items():
            assert (polynomial.bit_count() & 1) == (a == 0)
            by_horizontal[(S ^ a) & (v - 1)].append((S, a, polynomial))
            for y in range(d):
                u = products[y][a]
                bits = polynomial
                while bits:
                    bit = bits & -bits
                    matrices[bit.bit_length() - 1][y] ^= 1 << u
                    bits ^= bit
        for i, matrix in enumerate(matrices):
            row_ranks[i] += binary_rank(matrix)
    column_ranks = [0] * (maximum + 1)
    block_order = v * z * z
    dense = [[0] * (d * d) for _ in row_ranks] if dense_check else None
    for horizontal, entries in enumerate(by_horizontal):
        matrices = [[0] * block_order for _ in row_ranks]
        for S, a, polynomial in entries:
            for y in range(d):
                x, u = products[S][y], products[y][a]
                vv = products[x][a]
                xh, yh = x & (v - 1), y & (v - 1)
                assert xh ^ (u & (v - 1)) == yh ^ (vv & (v - 1)) == horizontal
                row = z * z * xh + z * (x // v) + u // v
                col = z * z * yh + z * (y // v) + vv // v
                bits = polynomial
                while bits:
                    bit = bits & -bits
                    i = bit.bit_length() - 1
                    matrices[i][row] ^= 1 << col
                    if dense is not None:
                        dense[i][x * d + u] ^= 1 << (y * d + vv)
                    bits ^= bit
        for row in range(block_order):
            xh, pair = divmod(row, z * z)
            xc, uc = divmod(pair, z)
            expected = z * z * (xh ^ horizontal) + z * uc + xc
            summed = 0
            for matrix in matrices:
                summed ^= matrix[row]
            assert summed == 1 << expected
        for i, matrix in enumerate(matrices):
            column_ranks[i] += binary_rank(matrix)
    if dense is not None:
        assert column_ranks == [binary_rank(matrix) for matrix in dense]
    assert column_ranks[0] == constant_rank_polynomial(group)["constant_column_rank"]
    upper = row_rank_upper_bounds(group)
    assert all(r <= bound for r, bound in zip(row_ranks, upper))
    return {"horizontal_pairs": group.m, "center_bits": center, "columns": list(columns),
            "labels": d, "label_bits": group.width, "weights": list(weights),
            "row_ranks": row_ranks, "column_ranks": column_ranks,
            "row_upper_bounds": upper,
            "maximum_rank_product": max(r * c for r, c in zip(row_ranks, column_ranks)),
            "column_block_count": v, "column_block_order": block_order,
            "dense_column_comparison": dense_check, "complete_crossing_entries": d ** 4}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--larger", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [profile((1,), 1, True), profile((1, 1), 1, True), profile((1, 2), 2)]
    if args.larger:
        cases += [profile((1, 1, 1), 1), profile((1, 1, 2), 2), profile((1, 2, 3), 2)]
    result = {"status": "pass", "scope": "exact finite complete coded-group crossings; no improved factoring exponent",
              "cases": cases, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
