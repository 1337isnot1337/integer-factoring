#!/usr/bin/env python3
"""Exact controls for the all-size nonresonant Cauchy Smith reduction.

This is a finite checker, not a compressed evaluator or a factoring algorithm.
All field operations are over Q[w]/(w*w+w+1); no floating point is used.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path


def v2(q: Fraction) -> int | None:
    if not q:
        return None
    numerator, denominator = abs(q.numerator), q.denominator
    return (numerator & -numerator).bit_length() - (denominator & -denominator).bit_length()


@dataclass(frozen=True)
class Quad:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def of(value: Quad | int | Fraction) -> Quad:
        return value if isinstance(value, Quad) else Quad(Fraction(value))

    def __add__(self, other: Quad | int | Fraction) -> Quad:
        rhs = Quad.of(other)
        return Quad(self.a + rhs.a, self.b + rhs.b)

    __radd__ = __add__

    def __neg__(self) -> Quad:
        return Quad(-self.a, -self.b)

    def __sub__(self, other: Quad | int | Fraction) -> Quad:
        return self + -Quad.of(other)

    def __rsub__(self, other: Quad | int | Fraction) -> Quad:
        return Quad.of(other) + -self

    def __mul__(self, other: Quad | int | Fraction) -> Quad:
        rhs = Quad.of(other)
        return Quad(self.a * rhs.a - self.b * rhs.b,
                    self.a * rhs.b + self.b * rhs.a - self.b * rhs.b)

    __rmul__ = __mul__

    def inverse(self) -> Quad:
        norm = self.a * self.a - self.a * self.b + self.b * self.b
        if not norm:
            raise ZeroDivisionError
        return Quad((self.a - self.b) / norm, -self.b / norm)

    def __truediv__(self, other: Quad | int | Fraction) -> Quad:
        return self * Quad.of(other).inverse()

    def __rtruediv__(self, other: Quad | int | Fraction) -> Quad:
        return Quad.of(other) * self.inverse()

    def __bool__(self) -> bool:
        return bool(self.a or self.b)

    def valuation(self) -> int | None:
        values = [value for value in (v2(self.a), v2(self.b)) if value is not None]
        return min(values) if values else None


W = Quad(Fraction(0), Fraction(1))


def smith_orders(matrix: list[list[Quad]]) -> list[int]:
    """DVR pivot elimination; independent of the Cauchy minor formula."""
    if not matrix:
        return []
    a = [row.copy() for row in matrix]
    m, n = len(a), len(a[0])
    orders = []
    for start in range(min(m, n)):
        candidates = [(a[i][j].valuation(), i, j)
                      for i in range(start, m) for j in range(start, n) if a[i][j]]
        if not candidates:
            raise ValueError("Rank-deficient matrix")
        value, row, col = min(candidates)
        assert value is not None
        orders.append(value)
        a[start], a[row] = a[row], a[start]
        for current in a:
            current[start], current[col] = current[col], current[start]
        pivot = a[start][start]
        for i in range(start + 1, m):
            multiplier = a[i][start] / pivot
            for j in range(start + 1, n):
                a[i][j] = a[i][j] - multiplier * a[start][j]
            a[i][start] = Quad()
    assert orders == sorted(orders)
    return orders


def matrix_at(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int) -> list[list[Quad]]:
    scale = (1 << k) * W
    return [[scale / (n + scale - a * b) for b in cols] for a in rows]


def nonresonant_data(n: int, rows: tuple[int, ...], cols: tuple[int, ...]):
    matched_rows = {a for a in rows if n % a == 0 and n // a in cols}
    matched_cols = {n // a for a in matched_rows}
    unmatched_rows = tuple(a for a in rows if a not in matched_rows)
    unmatched_cols = tuple(b for b in cols if b not in matched_cols)
    kernel = [[Quad(Fraction(1, n - a * b)) for b in unmatched_cols] for a in unmatched_rows]
    return len(matched_rows), unmatched_rows, unmatched_cols, smith_orders(kernel)


def verify_case(n: int, rows: tuple[int, ...], cols: tuple[int, ...]) -> dict:
    assert n % 2 and len(set(rows)) == len(rows) and len(set(cols)) == len(cols)
    h = max((abs(n - a * b).bit_length() for a in rows for b in cols), default=1)
    h = max(1, h)
    count, unmatched_rows, unmatched_cols, base_orders = nonresonant_data(n, rows, cols)
    records = []
    for k in (h, h + 1, h + 4):
        actual = smith_orders(matrix_at(n, rows, cols, k))
        expected = [0] * count + [k + value for value in base_orders]
        assert actual == expected, (n, rows, cols, k, actual, expected)
        records.append({"k": k, "orders": actual})
    return {"n": n, "rows": rows, "cols": cols, "h": h,
            "factor_count": count, "unmatched_rows": unmatched_rows,
            "unmatched_cols": unmatched_cols, "nonresonant_orders": base_orders,
            "scaled_profiles": records}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = []
    for n in range(3, 32, 2):
        for size in range(1, 5):
            for first in range(1, 5):
                rows = tuple(range(first, first + size))
                for second in range(1, 5):
                    for other_size in range(1, 5):
                        cols = tuple(range(second, second + other_size))
                        cases.append(verify_case(n, rows, cols))
    controls = [
        (15, (3, 4, 5), (3, 4, 5)),
        (45, (1, 3, 5, 9, 15), (1, 3, 5, 9, 15)),
        (105, (3, 5, 7, 15, 21), (3, 5, 7, 15, 21)),
        (21, (1, 2, 3, 7, 22), (1, 2, 3, 7, 23)),
        (105, (3, 5, 7), (3, 5, 7, 15, 21, 35, 106)),
        (105, (3, 5, 7, 15, 21, 35, 106), (3, 5, 7)),
    ]
    witnesses = [verify_case(*case) for case in controls]
    result = {"theorem": "resonant unit-pivot Smith stabilization",
              "arithmetic": "exact Q[w]/(w^2+w+1)",
              "exhaustive_interval_cases": len(cases),
              "profile_checks": 3 * (len(cases) + len(witnesses)),
              "witnesses": witnesses,
              "scope": "Finite controls. No compressed evaluator and no improved factoring score."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
