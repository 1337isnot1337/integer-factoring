#!/usr/bin/env python3
"""Exact residual weighted-cut certificate, extending the cutoff to 50."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
from random import Random

from check_weighted_cut_rank_budget import (
    column_flattening,
    crossing,
    encoded_fraction,
    exact_above_sqrt2_power,
    exact_below_sqrt2_power,
    lower_interpolant,
    matrix_rank,
    product_operator,
    quotient_complement,
    random_matrix,
    weights,
)


def chord_coefficients(k: int, s: int) -> tuple[Fraction, Fraction]:
    assert k >= 2 and 1 <= s < k
    denominator = k**4 - 1
    mu = Fraction(s * s * k * k - 1, denominator)
    nu = Fraction(k * k - s * s, denominator)
    assert mu >= 0 and nu >= 0
    return mu, nu


def case_parameters(d: int, k: int, b: int) -> dict[str, int]:
    n = d * d
    assert 2 <= k <= d and 1 <= b <= (n - 1) // ((k - 1) ** 2)
    a_lower = (k - 1) ** 2 * b + 1
    a_upper = min(k * k * b, n)
    f = quotient_complement(d, b)
    row_rest = max(n - a_upper, f)
    col_rest = n - b
    if b == 1:
        row_rest = max(row_rest, n - k * (k + 1) // 2)
        extra = (k * k - 3 * k + 4) // 2
    else:
        extra = max(0, a_lower + f - n)
    assert a_lower <= a_upper and min(row_rest, col_rest, extra) >= 0
    return {
        "a_lower": a_lower,
        "a_upper": a_upper,
        "row_rest": row_rest,
        "col_rest": col_rest,
        "extra_rank_budget": extra,
    }


def case_bound(d: int, k: int, b: int) -> Fraction:
    parameters = case_parameters(d, k, b)
    a = parameters["a_lower"]
    rr, cc = parameters["row_rest"], parameters["col_rest"]
    value = Fraction(2 * d * d + parameters["extra_rank_budget"], 2 * k - 1)
    for s, weight in enumerate(weights(k), start=1):
        mu, nu = chord_coefficients(k, s)
        first = max(Fraction(s * d), min(a, s * s * b) + mu * rr + nu * cc)
        second = max(Fraction(s * d), b + nu * rr + mu * cc)
        value += weight * (first + second)
    return value


def dimension_certificate(d: int) -> dict[str, object]:
    best = Fraction(d * d)  # Separate k=1 case.
    selected: dict[str, object] = {"k": 1, "b": None}
    case_count = 1
    by_k: list[dict[str, object]] = []
    for k in range(2, d + 1):
        local_best: Fraction | None = None
        local_b = 0
        upper_b = (d * d - 1) // ((k - 1) ** 2)
        for b in range(1, upper_b + 1):
            bound = case_bound(d, k, b)
            case_count += 1
            if local_best is None or bound < local_best:
                local_best, local_b = bound, b
        assert local_best is not None
        by_k.append({"k": k, "b": local_b, "bound": encoded_fraction(local_best), "b_cases": upper_b})
        if local_best < best:
            best = local_best
            selected = {"k": k, "b": local_b, **case_parameters(d, k, local_b)}
    return {
        "d": d,
        "bound": encoded_fraction(best),
        "above_target": exact_above_sqrt2_power(best, d),
        "below_target": exact_below_sqrt2_power(best, d),
        "minimizing_case": selected,
        "scalar_cases": case_count,
        "minimum_by_k": by_k,
    }


def verify_chords() -> dict[str, int]:
    inequalities = 0
    endpoint_equalities = 0
    for k in range(2, 21):
        for s in range(1, k):
            mu, nu = chord_coefficients(k, s)
            assert mu * k * k + nu == s * s
            assert mu + nu * k * k == 1
            endpoint_equalities += 2
            for r in range(1, min(k * k, 16) + 1):
                for c in range(1, min(k * k, 16) + 1):
                    assert min(r, s * s * c) >= mu * r + nu * c
                    assert min(s * s * r, c) >= nu * r + mu * c
                    inequalities += 2
    return {"endpoint_equalities": endpoint_equalities, "integer_pair_inequalities": inequalities}


def verify_profile(d: int, profile: list[tuple[int, int]]) -> int:
    nonzero = [(r, c) for r, c in profile if r or c]
    assert all(r > 0 and c > 0 for r, c in nonzero)
    ratio, index, exchange = max(
        [(Fraction(r, c), i, False) for i, (r, c) in enumerate(nonzero)]
        + [(Fraction(c, r), i, True) for i, (r, c) in enumerate(nonzero)]
    )
    if ratio == 1:
        assert all(r == c for r, c in nonzero)
        assert sum(r for r, _ in nonzero) >= d * d
        return 0
    if exchange:
        nonzero = [(c, r) for r, c in nonzero]
    a, b = nonzero[index]
    k = isqrt((a - 1) // b) + 1
    assert (k - 1) ** 2 * b < a <= k * k * b
    parameters = case_parameters(d, k, b)
    rest = nonzero[:index] + nonzero[index + 1:]
    assert sum(r for r, _ in rest) >= parameters["row_rest"]
    assert sum(c for _, c in rest) >= parameters["col_rest"]
    assert sum(r + c for r, c in nonzero) >= 2 * d * d + parameters["extra_rank_budget"]
    for s in range(1, k):
        mu, nu = chord_coefficients(k, s)
        row_lower = max(
            Fraction(s * d),
            min(parameters["a_lower"], s * s * b)
            + mu * parameters["row_rest"] + nu * parameters["col_rest"],
        )
        col_lower = max(
            Fraction(s * d),
            b + nu * parameters["row_rest"] + mu * parameters["col_rest"],
        )
        assert sum(min(r, s * s * c) for r, c in nonzero) >= row_lower
        assert sum(min(s * s * r, c) for r, c in nonzero) >= col_lower
    # This checks the complete rational certificate against the sum of
    # actual interpolant values, before the separate square-root bound.
    actual_dual = sum((lower_interpolant(k, r, c) for r, c in nonzero), Fraction())
    assert actual_dual >= case_bound(d, k, b)
    return 1


def verify_finite_decompositions() -> dict[str, object]:
    rng = Random(50120261030)
    profiles = 0
    nontrivial = 0
    for p in (2, 3, 5):
        for d in range(2, 6):
            n = d * d
            swap = crossing(d)
            for terms in (1, 2, 3):
                for _ in range(10):
                    products = [product_operator(random_matrix(d, p, rng), random_matrix(d, p, rng), p)
                                for _ in range(terms)]
                    first = [[sum(t[i][j] for t in products) % p for j in range(n)] for i in range(n)]
                    second = [[(swap[i][j] - first[i][j]) % p for j in range(n)] for i in range(n)]
                    profile = [(matrix_rank(t, p), matrix_rank(column_flattening(t, d), p))
                               for t in (first, second)]
                    nontrivial += verify_profile(d, profile)
                    profiles += 1
            # A three-piece decomposition with both a row-rank-one and
            # a column-rank-one piece exercises mixed extreme cuts.
            for _ in range(10):
                row_left = [rng.randrange(p) for _ in range(n)]
                row_right = [rng.randrange(p) for _ in range(n)]
                first = [[a * b % p for b in row_right] for a in row_left]
                second = product_operator(random_matrix(d, p, rng), random_matrix(d, p, rng), p)
                third = [[(swap[i][j] - first[i][j] - second[i][j]) % p for j in range(n)] for i in range(n)]
                profile = [(matrix_rank(t, p), matrix_rank(column_flattening(t, d), p))
                           for t in (first, second, third)]
                nontrivial += verify_profile(d, profile)
                profiles += 1
    return {"primes": [2, 3, 5], "dimensions": [2, 3, 4, 5],
            "decompositions": profiles, "nontrivial_ratio_classes": nontrivial,
            "random_seed": 50120261030}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 577 * 577 - 2 * 408 * 408 == 1
    assert 1393 * 1393 - 2 * 985 * 985 == -1
    chords = verify_chords()
    fields = verify_finite_decompositions()
    dimensions = [dimension_certificate(d) for d in range(47, 52)]
    assert all(row["above_target"] for row in dimensions[:-1])
    assert dimensions[-1]["below_target"]
    assert not dimensions[-1]["above_target"]
    here = Path(__file__)
    result = {
        "status": "PASS",
        "scope": "all-field finite rank-sum exclusion through 50; no new factoring exponent",
        "previous_cutoff": 46,
        "last_certified_dimension": 50,
        "first_uncertified_dimension": 51,
        "checker_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
        "base_checker_sha256": hashlib.sha256(here.with_name("check_weighted_cut_rank_budget.py").read_bytes()).hexdigest(),
        "chord_controls": chords,
        "finite_field_controls": fields,
        "scalar_cases": sum(row["scalar_cases"] for row in dimensions),
        "dimensions": dimensions,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "last_certified_dimension": 50,
        "first_uncertified_dimension": 51,
        "scalar_cases": result["scalar_cases"],
        "bounds": [{"d": row["d"], "bound": row["bound"], "class": row["minimizing_case"]}
                   for row in dimensions],
        "chord_controls": chords,
        "finite_field_controls": fields,
        "output": str(args.output) if args.output else None,
    }, indent=2))


if __name__ == "__main__":
    main()
