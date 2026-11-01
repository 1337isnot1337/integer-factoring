#!/usr/bin/env python3
"""Exact conditioning, padding, finite-precision trace, and inverse controls."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import random

from check_pole_smith_stabilization import Quad, matrix_at, smith_orders
from pole_cauchy_conditioning import (balanced_padding, inverse_entry_valuation,
                                      normalized_inverse_entry, residue_count)


def expand(domain: tuple) -> tuple[int, ...]:
    return tuple(first + 2 * i for first, size in domain for i in range(size))


def inverse(matrix: list[list[Quad]]) -> list[list[Quad]]:
    size = len(matrix)
    a = [row[:] + [Quad(Fraction(int(i == j))) for j in range(size)]
         for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if a[i][column])
        a[column], a[pivot] = a[pivot], a[column]
        divisor = a[column][column]
        a[column] = [value / divisor for value in a[column]]
        for row in range(size):
            if row != column:
                multiple = a[row][column]
                a[row] = [x - multiple * y for x, y in zip(a[row], a[column])]
    return [row[size:] for row in a]


def mod_pair(value: Quad, precision: int) -> tuple[int, int]:
    modulus = 1 << precision
    return tuple(coordinate.numerator * pow(coordinate.denominator, -1, modulus) % modulus
                 for coordinate in (value.a, value.b))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(20261101)
    matrices = 0
    entry_checks = 0
    trace_checks = 0
    sharp = None
    for _ in range(200):
        n = rng.randrange(17, 150, 2)
        left = (rng.randrange(1, 8, 2), rng.randrange(1, 5))
        right = (rng.randrange(1, 8, 2), rng.randrange(1, 5))
        aa, bb, padding_modulus = balanced_padding(n, left, right)
        rows, cols = expand(aa), expand(bb)
        original_count = sum(a * b == n for a in expand((left,)) for b in expand((right,)))
        assert original_count == sum(a * b == n for a in rows for b in cols)
        h = max(1, max(abs(n - a * b) for a in rows for b in cols).bit_length())
        k = h + 4
        matrix = matrix_at(n, rows, cols, k)
        inv = inverse(matrix)
        profile = smith_orders(matrix)
        assert profile[-1] <= k - 1
        if profile[-1] == k - 1 and sharp is None:
            sharp = {"n": n, "rows": rows, "cols": cols, "k": k, "orders": profile}
        for row, a in enumerate(rows):
            for col, b in enumerate(cols):
                actual = inv[col][row]
                value = actual.valuation()
                assert value is not None and value >= 1 - k
                assert value == inverse_entry_valuation(n, aa, bb, a, b, k)
                normalized = normalized_inverse_entry(n, aa, bb, a, b, k, 8)
                assert normalized["valuation"] == value
                scaled = actual / (1 << value) if value >= 0 else actual * (1 << -value)
                assert tuple(normalized["unit"]) == mod_pair(scaled, 8)
                entry_checks += 1
        p = 4
        raw_precision = p + 2 * (k - 1)
        rounded = [[Quad(*(Fraction(v) for v in mod_pair(value, raw_precision)))
                    for value in row] for row in matrix]
        rounded_inverse = inverse(rounded)
        scalar = sum((rounded_inverse[j][i] * rounded[i][j] * rounded[i][j]
                      for i in range(len(rows)) for j in range(len(cols))), Quad())
        assert mod_pair(scalar, p) == (original_count % (1 << p), 0)
        trace_checks += 1
        for level in range(1, h + 1):
            modulus = 1 << level
            for residue in {a % modulus for a in rows} | {n * pow(b, -1, modulus) % modulus for b in cols}:
                delta = residue_count(aa, residue, level) - residue_count(bb, n * pow(residue, -1, modulus), level)
                assert abs(delta) <= 1
        matrices += 1
    n = (1 << 401) + 45
    aa, bb, modulus = balanced_padding(n, (1, (1 << 180) + 7), (3, (1 << 200) + 13))
    max_row = max(first + 2 * (size - 1) for first, size in aa)
    max_col = max(first + 2 * (size - 1) for first, size in bb)
    h = max(abs(n - aa[0][0] * bb[0][0]), abs(n - max_row * max_col)).bit_length()
    k = h + 12
    queries = ((aa[0][0], bb[0][0]), (aa[-1][0], bb[0][0] + 2),
               (1 + (1 << 91), 3 + (1 << 151)))
    long_controls = []
    for a, b in queries:
        result = normalized_inverse_entry(n, aa, bb, a, b, k, 12)
        assert result["valuation"] >= 1 - k
        long_controls.append({"row": a, "column": b, **result})
    result = {"claim": "Balanced Cauchy padding and dimension-independent inverse precision",
              "padded_matrix_checks": matrices, "exact_inverse_entry_checks": entry_checks,
              "finite_precision_trace_checks": trace_checks, "sharp_small_control": sharp,
              "long_implicit_control": {"input_bits": n.bit_length(), "domain_size": (1 << 200) + 13,
                                        "padding_modulus": modulus, "k": k, "left": aa,
                                        "right": bb, "queries": long_controls},
              "scope": "Polynomial-cost individual inverse entries; global trace summation remains open."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
