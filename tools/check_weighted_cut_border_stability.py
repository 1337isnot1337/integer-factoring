#!/usr/bin/env python3
"""Exact Laurent-polynomial cancellation controls for the DVR proof."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from check_weighted_cut_rank_budget import (
    column_flattening, crossing, exact_above_sqrt2_power,
    harmonic_bound, matrix_rank,
)

Poly = dict[int, int]
Matrix = list[list[Poly]]


def add(a: Poly, b: Poly, p: int) -> Poly:
    result = dict(a)
    for degree, coefficient in b.items():
        result[degree] = (result.get(degree, 0) + coefficient) % p
        if not result[degree]:
            del result[degree]
    return result


def multiply(a: Poly, b: Poly, p: int) -> Poly:
    result: Poly = {}
    for i, x in a.items():
        for j, y in b.items():
            result[i + j] = (result.get(i + j, 0) + x * y) % p
    return {degree: value for degree, value in result.items() if value}


def constant(matrix: list[list[int]], p: int) -> Matrix:
    return [[{0: x % p} if x % p else {} for x in row] for row in matrix]


def identity(d: int, p: int) -> Matrix:
    return constant([[int(i == j) for j in range(d)] for i in range(d)], p)


def matrix_add(a: Matrix, b: Matrix, p: int) -> Matrix:
    return [[add(x, y, p) for x, y in zip(arow, brow)] for arow, brow in zip(a, b)]


def shift(a: Matrix, degree: int, scalar: int, p: int) -> Matrix:
    return [[{j + degree: scalar * x % p for j, x in entry.items() if scalar * x % p}
             for entry in row] for row in a]


def compose(a: Matrix, b: Matrix, p: int) -> Matrix:
    result: Matrix = [[{} for _ in b[0]] for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if not x:
                continue
            for j, y in enumerate(b[k]):
                if y:
                    result[i][j] = add(result[i][j], multiply(x, y, p), p)
    return result


def kron(a: Matrix, b: Matrix, p: int) -> Matrix:
    return [[multiply(x, y, p) for x in ar for y in br] for ar in a for br in b]


def reduction(a: Matrix) -> list[list[int]]:
    assert all(degree >= 0 for row in a for entry in row for degree in entry)
    return [[entry.get(0, 0) for entry in row] for row in a]


def at_one(a: Matrix, p: int) -> list[list[int]]:
    return [[sum(entry.values()) % p for entry in row] for row in a]


def verify_case(d: int, s: int, k: int, pole: int, p: int) -> None:
    ell = k * s
    assert 0 < ell < d
    unit = identity(d, p)
    g = [[add(unit[i][j], {1: 1} if i == j + 1 else {}, p)
          for j in range(d)] for i in range(d)]
    inverse = [[{i - j: (-1)**(i - j) % p} if i >= j else {}
                for j in range(d)] for i in range(d)]
    assert compose(g, inverse, p) == unit == compose(inverse, g, p)
    quotient = inverse[ell:]
    assert matrix_rank(reduction(quotient), p) == d - ell
    b: Matrix = [[{} for _ in range(d * d)] for _ in range(d * d)]
    for j in range(k):
        selector = constant([[int(v < s and u == j * s + v)
                              for v in range(d)] for u in range(d)], p)
        first = shift(compose(g, selector, p), -pole - j, 1, p)
        second = unit if j == 0 else constant(
            [[int(u == (v + 1) % d) for v in range(d)] for u in range(d)], p)
        b = matrix_add(b, kron(first, second, p), p)
    # The two constant second factors and disjoint first-factor blocks
    # are independent.  This finite specialization verifies c(B)=k.
    assert matrix_rank(column_flattening(at_one(b, p), d), p) == k
    n = d * d
    perturbation = constant([[(i + 2 * j + 1) % p for j in range(n)] for i in range(n)], p)
    total = matrix_add(constant(crossing(d), p), shift(perturbation, 1, 1, p), p)
    a = matrix_add(total, shift(b, 0, -1, p), p)
    assert matrix_add(a, b, p) == total
    assert any(degree < 0 for row in a for entry in row for degree in entry)
    assert any(degree < 0 for row in b for entry in row for degree in entry)
    projection = kron(quotient, unit, p)
    restricted_a = [row[:s * d] for row in a]
    restricted_b = [row[:s * d] for row in b]
    restricted_total = [row[:s * d] for row in total]
    assert all(not entry for row in compose(projection, restricted_b, p) for entry in row)
    projected = compose(projection, restricted_a, p)
    assert projected == compose(projection, restricted_total, p)
    assert matrix_rank(reduction(projected), p) == s * (d - ell)
    assert s * (d - ell) + s * s * k == s * d
    # Negative control for reducing an unsaturated set of generators.
    saturated_generators = [row[:ell] for row in g]
    unsaturated_generators = shift(saturated_generators, 1, 1, p)
    assert matrix_rank(reduction(saturated_generators), p) == ell
    assert matrix_rank(reduction(unsaturated_generators), p) == 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = 0
    for p in (2, 3, 5):
        for d in range(2, 7):
            for s in range(1, d):
                for k in range(1, min(2, (d - 1) // s) + 1):
                    for pole in (1, 3):
                        verify_case(d, s, k, pole, p)
                        cases += 1
    assert 577 * 577 - 2 * 408 * 408 == 1
    for d in range(2, 42):
        assert exact_above_sqrt2_power(harmonic_bound(d), d)
    result = {
        "status": "PASS", "scope": "supplied one-parameter Laurent-series degenerations",
        "finite_field_cancellation_cases": cases,
        "fields": [2, 3, 5], "dimensions": [2, 3, 4, 5, 6],
        "product_counts": [1, 2], "initial_pole_orders": [1, 3],
        "unsaturated_negative_controls": cases,
        "harmonic_dimensions_certified": list(range(2, 42)),
        "no_extension_of_refined_cutoff_asserted": True,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_checker_sha256": hashlib.sha256(Path(__file__).with_name(
            "check_weighted_cut_rank_budget.py").read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "harmonic_dimensions_certified"}, indent=2))


if __name__ == "__main__":
    main()
