#!/usr/bin/env python3
"""Exact weighted-cut certificates for finite crossing decompositions.

Only Python's standard library is used.  The universal proof is in
research/weighted-cut-rank-budget.md; finite field examples are controls,
not a substitute for that proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
from random import Random


def encoded_fraction(value: Fraction) -> dict[str, object]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": float(value),
    }


@lru_cache(maxsize=None)
def weights(k: int) -> tuple[Fraction, ...]:
    if k < 2:
        raise ValueError("the interpolant requires k >= 2")
    return (Fraction(1, 6),) + tuple(
        Fraction(2, (2 * s - 1) * (2 * s + 1)) for s in range(2, k)
    )


def lower_interpolant(k: int, r: int, c: int) -> Fraction:
    total = Fraction(r + c, 2 * k - 1)
    for s, weight in enumerate(weights(k), start=1):
        total += weight * (min(r, s * s * c) + min(s * s * r, c))
    return total


@lru_cache(maxsize=None)
def odd_harmonic(k: int) -> Fraction:
    return sum((Fraction(1, 2 * s - 1) for s in range(1, k + 1)), Fraction())


def harmonic_bound(d: int) -> Fraction:
    return d * (2 * odd_harmonic(d) - 1)


def truncated_bound(d: int, k: int) -> Fraction:
    if not 2 <= k <= d:
        raise ValueError("B(d,k) is used only for 2 <= k <= d")
    closed = d * (2 * odd_harmonic(k) - 1) + Fraction(2 * d * (d - k), 2 * k - 1)
    direct = Fraction(2 * d * d, 2 * k - 1) + 2 * d * sum(
        (s * weight for s, weight in enumerate(weights(k), start=1)), Fraction()
    )
    assert direct == closed
    return closed


def quotient_complement(d: int, b: int) -> int:
    """max_{integer 0 <= s <= d} s(d-bs), by concavity."""
    if b < 1:
        raise ValueError("positive opposite-cut rank required")
    s = d // (2 * b)
    return max(0, s * (d - b * s), (s + 1) * (d - b * (s + 1)))


def extra_rank_budget(d: int, k: int) -> tuple[int, int, int]:
    """Return (extra total, extremal small rank, checked scalar cases)."""
    n = d * d
    best_extra = (k * k - 3 * k + 4) // 2
    assert 2 * best_extra == k * k - 3 * k + 4
    best_b = 1
    cases = 1
    for b in range(2, (n - 1) // ((k - 1) ** 2) + 1):
        a_lower = (k - 1) ** 2 * b + 1
        assert a_lower <= n and a_lower <= k * k * b
        extra = max(0, a_lower + quotient_complement(d, b) - n)
        cases += 1
        if extra < best_extra:
            best_extra, best_b = extra, b
    return best_extra, best_b, cases


def exact_above_sqrt2_power(value: Fraction, d: int) -> bool:
    # 577/408 > sqrt(2), certified by the Pell equation below.
    return value.numerator**408 > d**577 * value.denominator**408


def exact_below_sqrt2_power(value: Fraction, d: int) -> bool:
    # 1393/985 < sqrt(2).
    return value.numerator**985 < d**1393 * value.denominator**985


def dimension_certificate(d: int) -> dict[str, object]:
    # The k=1 case is separate.  Substituting k=1 into B(d,k)
    # would give an invalid bound and is deliberately forbidden above.
    best = Fraction(d * d)
    selected = {"k": 1, "extra_rank_budget": 0, "extremal_small_rank": None}
    scalar_cases = 1
    candidates: list[dict[str, object]] = []
    for k in range(2, d + 1):
        extra, b, checked = extra_rank_budget(d, k)
        scalar_cases += checked
        value = truncated_bound(d, k) + Fraction(extra, 2 * k - 1)
        candidates.append({
            "k": k,
            "extra_rank_budget": extra,
            "extremal_small_rank": b,
            "bound": encoded_fraction(value),
        })
        if value < best:
            best = value
            selected = {"k": k, "extra_rank_budget": extra, "extremal_small_rank": b}
    return {
        "d": d,
        "harmonic_bound": encoded_fraction(harmonic_bound(d)),
        "harmonic_above_target": exact_above_sqrt2_power(harmonic_bound(d), d),
        "bound": encoded_fraction(best),
        "above_target": exact_above_sqrt2_power(best, d),
        "below_target": exact_below_sqrt2_power(best, d),
        "minimizing_case": selected,
        "scalar_cases": scalar_cases,
        "candidates": candidates,
    }


def matrix_rank(matrix: list[list[int]], p: int) -> int:
    rows = [[x % p for x in row] for row in matrix]
    height = len(rows)
    width = len(rows[0]) if rows else 0
    pivot = 0
    for col in range(width):
        found = next((i for i in range(pivot, height) if rows[i][col]), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        inverse = pow(rows[pivot][col], -1, p)
        rows[pivot] = [(v * inverse) % p for v in rows[pivot]]
        for i in range(pivot + 1, height):
            multiple = rows[i][col]
            if multiple:
                rows[i] = [(a - multiple * b) % p for a, b in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == height:
            break
    return pivot


def crossing(d: int) -> list[list[int]]:
    return [[int(x == v and y == u) for u in range(d) for v in range(d)]
            for x in range(d) for y in range(d)]


def column_flattening(row_matrix: list[list[int]], d: int) -> list[list[int]]:
    return [[row_matrix[x * d + y][u * d + v] for y in range(d) for v in range(d)]
            for x in range(d) for u in range(d)]


def product_operator(left: list[list[int]], right: list[list[int]], p: int) -> list[list[int]]:
    d = len(left)
    return [[left[x][u] * right[y][v] % p for u in range(d) for v in range(d)]
            for x in range(d) for y in range(d)]


def random_matrix(d: int, p: int, rng: Random) -> list[list[int]]:
    return [[rng.randrange(p) for _ in range(d)] for _ in range(d)]


def verify_interpolants() -> dict[str, int]:
    knots = 0
    strict_midpoints = 0
    integer_pairs = 0
    for k in range(2, 48):
        for j in range(1, k + 1):
            assert lower_interpolant(k, j * j, 1) == j
            assert lower_interpolant(k, 1, j * j) == j
            knots += 2
        for j in range(1, k):
            r, c = j * j + (j + 1) ** 2, 2
            lower = lower_interpolant(k, r, c)
            assert lower > 0 and lower * lower < r * c
            assert lower_interpolant(k, c, r) == lower
            strict_midpoints += 1
    for k in range(2, 13):
        for r in range(1, min(k * k, 15) + 1):
            for c in range(1, min(k * k, 15) + 1):
                lower = lower_interpolant(k, r, c)
                assert lower >= 0 and lower * lower <= r * c
                integer_pairs += 1
    return {"knots": knots, "strict_midpoints": strict_midpoints, "integer_pairs": integer_pairs}


def verify_scalar_formulas() -> dict[str, object]:
    quadratic_cases = 0
    rank_one_cases = 0
    for d in range(2, 21):
        for b in range(1, d * d + 1):
            assert quotient_complement(d, b) == max(s * (d - b * s) for s in range(d + 1))
            quadratic_cases += 1
    for k in range(2, 48):
        guaranteed = (k * k - 3 * k + 4) // 2
        for a in range((k - 1) ** 2 + 1, k * k + 1):
            h = isqrt(a)
            assert a - h * (h + 1) // 2 >= guaranteed
            rank_one_cases += 1
    fake_profile = [(18, 18), (2, 464), (464, 2)]
    assert sum(r for r, _ in fake_profile) == 22**2
    assert sum(c for _, c in fake_profile) == 22**2
    assert sum(min(r, c) for r, c in fake_profile) == 22
    weighted = sum(min(r, 4 * c) for r, c in fake_profile)
    assert weighted == 28 and weighted < 44
    try:
        truncated_bound(4, 1)
    except ValueError:
        pass
    else:
        raise AssertionError("k=1 must use its separate proof")
    return {
        "quadratic_maxima": quadratic_cases,
        "rank_one_intervals": rank_one_cases,
        "old_budget_negative_control": {"d": 22, "s": 2, "left": weighted, "required": 44},
        "k_one_guard": True,
    }


def verify_finite_fields() -> dict[str, object]:
    rng = Random(20261030)
    quotient_cases = 0
    rank_one_cases = 0
    decomposition_cases = 0
    weighted_inequalities = 0
    for p in (2, 3, 5):
        for d in range(2, 6):
            n = d * d
            swap = crossing(d)
            for terms in range(4):
                for _ in range(6):
                    products = [product_operator(random_matrix(d, p, rng), random_matrix(d, p, rng), p)
                                for _ in range(terms)]
                    b = [[sum(t[i][j] for t in products) % p for j in range(n)] for i in range(n)]
                    residual = [[(swap[i][j] - b[i][j]) % p for j in range(n)] for i in range(n)]
                    residual_rank = matrix_rank(residual, p)
                    opposite_rank = matrix_rank(column_flattening(b, d), p)
                    assert opposite_rank <= terms
                    for s in range(d + 1):
                        assert residual_rank + s * s * opposite_rank >= s * d
                        quotient_cases += 1
                    if opposite_rank == 1:
                        a = matrix_rank(b, p)
                        h = isqrt(a)
                        assert residual_rank >= n - h * (h + 1) // 2
                        rank_one_cases += 1
            for count in (2, 3, 4):
                for _ in range(10):
                    pieces = [random_matrix(n, p, rng) for _ in range(count - 1)]
                    pieces.append([[(swap[i][j] - sum(t[i][j] for t in pieces)) % p
                                    for j in range(n)] for i in range(n)])
                    profile = [(matrix_rank(t, p), matrix_rank(column_flattening(t, d), p)) for t in pieces]
                    assert all((r == 0) == (c == 0) for r, c in profile)
                    for s in range(1, d + 1):
                        assert sum(min(r, s * s * c) for r, c in profile) >= s * d
                        assert sum(min(s * s * r, c) for r, c in profile) >= s * d
                        weighted_inequalities += 2
                    decomposition_cases += 1
    return {
        "primes": [2, 3, 5],
        "dimensions": [2, 3, 4, 5],
        "quotient_inequalities": quotient_cases,
        "rank_one_complements": rank_one_cases,
        "exact_decompositions": decomposition_cases,
        "weighted_inequalities": weighted_inequalities,
        "random_seed": 20261030,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 577 * 577 - 2 * 408 * 408 == 1
    assert 1393 * 1393 - 2 * 985 * 985 == -1
    interpolants = verify_interpolants()
    scalar_formulas = verify_scalar_formulas()
    finite_fields = verify_finite_fields()
    dimensions = [dimension_certificate(d) for d in range(2, 48)]
    assert all(row["harmonic_above_target"] for row in dimensions if row["d"] <= 41)
    assert exact_below_sqrt2_power(harmonic_bound(42), 42)
    assert all(row["above_target"] for row in dimensions if row["d"] <= 46)
    assert dimensions[-1]["below_target"]
    assert not dimensions[-1]["above_target"]
    result = {
        "status": "PASS",
        "theorem": "weighted cut budgets and finite rank-sum exclusion through 46 labels",
        "scope": "all fields; finite rank-sum criterion only; no new factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "upper_sqrt2_pell_pair": [577, 408],
        "lower_sqrt2_pell_pair": [1393, 985],
        "interpolant_controls": interpolants,
        "scalar_controls": scalar_formulas,
        "finite_field_controls": finite_fields,
        "harmonic_last_certified_dimension": 41,
        "strengthened_last_certified_dimension": 46,
        "first_uncertified_dimension": 47,
        "scalar_cases": sum(row["scalar_cases"] for row in dimensions),
        "dimensions": dimensions,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "harmonic_last_certified_dimension": 41,
        "strengthened_last_certified_dimension": 46,
        "first_uncertified_dimension": 47,
        "scalar_cases": result["scalar_cases"],
        "bound_at_46": dimensions[-2]["bound"],
        "bound_at_47": dimensions[-1]["bound"],
        "interpolant_controls": interpolants,
        "finite_field_controls": finite_fields,
        "output": str(args.output) if args.output else None,
    }, indent=2))


if __name__ == "__main__":
    main()
