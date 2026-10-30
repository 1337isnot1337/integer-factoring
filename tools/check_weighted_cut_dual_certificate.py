#!/usr/bin/env python3
"""Replay positive rational cut duals excluding seeds through 52 labels.

The witnesses were discovered separately with a numerical optimizer.
This verifier uses only Python's standard library and exact arithmetic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import re

from check_weighted_cut_rank_budget import encoded_fraction, exact_above_sqrt2_power
from check_weighted_cut_residual_budget import case_bound, case_parameters


Key = tuple[int, int, int]


def parse_nonnegative_fraction(raw: object) -> Fraction:
    # Check the stored sign and denominator before Fraction normalizes them.
    assert isinstance(raw, str) and re.fullmatch(r"[0-9]+(?:/[1-9][0-9]*)?", raw)
    numerator, separator, denominator = raw.partition("/")
    assert int(numerator) >= 0
    assert not separator or int(denominator) > 0
    result = Fraction(raw)
    assert result >= 0 and result.denominator > 0
    return result


def dual_value(k: int, coefficients: tuple[Fraction, ...], r: int, c: int) -> Fraction:
    assert len(coefficients) == 2 * k and r > 0 and c > 0
    features = [r, c]
    features += [min(r, s * s * c) for s in range(1, k)]
    features += [min(s * s * r, c) for s in range(1, k)]
    return sum((a * x for a, x in zip(coefficients, features)), Fraction())


def verify_knots(k: int, coefficients: tuple[Fraction, ...]) -> dict[str, object]:
    assert k >= 2 and len(coefficients) == 2 * k
    assert all(a >= 0 and a.denominator > 0 for a in coefficients)
    slacks = [Fraction(j) - dual_value(k, coefficients, j * j, 1)
              for j in range(1, k + 1)]
    slacks += [Fraction(j) - dual_value(k, coefficients, 1, j * j)
               for j in range(2, k + 1)]
    assert len(slacks) == 2 * k - 1 and min(slacks) >= 0
    return {"knots": len(slacks), "minimum_slack": encoded_fraction(min(slacks)),
            "positive_coefficients": sum(a > 0 for a in coefficients)}


def checked_parameters(d: int, k: int, b: int) -> dict[str, int]:
    parameters = case_parameters(d, k, b)
    n = d * d
    assert 2 <= k <= d and 1 <= b <= (n - 1) // (k - 1) ** 2
    a0, a1 = parameters["a_lower"], parameters["a_upper"]
    assert a0 == (k - 1) ** 2 * b + 1
    assert a1 == min(k * k * b, n)
    assert b < a0 <= a1 <= n
    # Independently enumerate the quotient maximum used by the imported helper.
    f = max(s * (d - b * s) for s in range(d + 1))
    rr, cc = max(n - a1, f), n - b
    if b == 1:
        # The stable-rank swap/product complement theorem uses h <= k.
        rr = max(rr, n - k * (k + 1) // 2)
        extra = (k * k - 3 * k + 4) // 2
        assert 2 * extra == k * k - 3 * k + 4
    else:
        extra = max(0, a0 + f - n)
    assert parameters["row_rest"] == rr >= 0
    assert parameters["col_rest"] == cc >= 0
    assert parameters["extra_rank_budget"] == extra >= 0
    return parameters


def witness_bound(key: Key, coefficients: tuple[Fraction, ...],
                  sqrt_denominator: int) -> tuple[Fraction, dict[str, object]]:
    d, k, b = key
    knots = verify_knots(k, coefficients)
    parameters = checked_parameters(d, k, b)
    a0, a1 = parameters["a_lower"], parameters["a_upper"]
    left = [max(0, s * d - min(a1, s * s * b)) for s in range(1, k)]
    right = [max(0, s * d - b) for s in range(1, k)]
    budgets = [parameters["row_rest"], parameters["col_rest"], *left, *right]
    assert len(budgets) == len(coefficients) == 2 * k
    assert all(x >= 0 for x in budgets)
    residual = sum((a * x for a, x in zip(coefficients, budgets)), Fraction())
    square = sqrt_denominator * sqrt_denominator * a0 * b
    numerator = isqrt(square)
    assert numerator * numerator <= square < (numerator + 1) ** 2
    piece = Fraction(numerator, sqrt_denominator)
    bound = piece + residual
    assert exact_above_sqrt2_power(bound, d)
    return bound, {
        "d": d, "k": k, "b": b,
        "parameters": parameters,
        "residual_budgets": budgets,
        "piece_lower_bound": encoded_fraction(piece),
        "residual_lower_bound": encoded_fraction(residual),
        "bound": encoded_fraction(bound),
        "above_target": True,
        **knots,
    }


def load_witnesses(path: Path) -> tuple[int, dict[Key, tuple[Fraction, ...]]]:
    document = json.loads(path.read_text())
    denominator = document["sqrt_denominator"]
    assert type(denominator) is int and denominator == 10**12
    witnesses: dict[Key, tuple[Fraction, ...]] = {}
    for row in document["witnesses"]:
        assert set(row) == {"d", "k", "b", "coefficients"}
        key = (row["d"], row["k"], row["b"])
        assert all(type(x) is int for x in key)
        assert key[0] in (51, 52) and key not in witnesses
        coefficients = tuple(parse_nonnegative_fraction(raw) for raw in row["coefficients"])
        assert len(coefficients) == 2 * key[1]
        witnesses[key] = coefficients
    assert len(witnesses) == 6
    return denominator, witnesses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witnesses", type=Path,
                        default=Path(__file__).resolve().parent.parent / "evidence" /
                        "weighted-cut-dual-witnesses.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 577 * 577 - 2 * 408 * 408 == 1
    assert Fraction(577, 408) > 1
    denominator, witnesses = load_witnesses(args.witnesses)
    needed: set[Key] = set()
    dimensions: list[dict[str, object]] = []
    certificate_rows: list[dict[str, object]] = []
    for d in (51, 52):
        n = d * d
        best = Fraction(n)  # k=1: every nonzero piece has equal cut ranks.
        assert exact_above_sqrt2_power(best, d)
        selected: dict[str, object] = {"k": 1, "b": None, "method": "equal ranks"}
        cases, inherited, new_duals, rank_one = 1, 0, 0, 0
        for k in range(2, d + 1):
            for b in range(1, (n - 1) // (k - 1) ** 2 + 1):
                key = (d, k, b)
                checked_parameters(d, k, b)
                lower = case_bound(d, k, b)
                method = "residual budget"
                cases += 1
                rank_one += b == 1
                if exact_above_sqrt2_power(lower, d):
                    inherited += 1
                else:
                    needed.add(key)
                    assert key in witnesses, f"missing dual for {key}"
                    lower, record = witness_bound(key, witnesses[key], denominator)
                    certificate_rows.append(record)
                    method = "positive rational dual"
                    new_duals += 1
                assert exact_above_sqrt2_power(lower, d)
                if lower < best:
                    best = lower
                    selected = {"k": k, "b": b, "method": method}
        assert rank_one == d - 1
        assert cases == 1 + inherited + new_duals
        assert exact_above_sqrt2_power(best, d)
        dimensions.append({
            "d": d, "bound": encoded_fraction(best), "above_target": True,
            "minimizing_case": selected, "scalar_cases": cases,
            "inherited_residual_cases": inherited, "new_dual_cases": new_duals,
            "rank_one_classes_checked": rank_one,
        })
    assert needed == set(witnesses), "unused or extra witness keys"
    assert needed == {(51, 34, 2), (51, 35, 2), (52, 34, 2),
                      (52, 35, 2), (52, 36, 2), (52, 37, 2)}
    here = Path(__file__).resolve()
    hashes = {name: hashlib.sha256(here.with_name(name).read_bytes()).hexdigest()
              for name in (here.name, "check_weighted_cut_rank_budget.py",
                           "check_weighted_cut_residual_budget.py")}
    hashes[args.witnesses.name] = hashlib.sha256(args.witnesses.read_bytes()).hexdigest()
    result = {
        "status": "PASS",
        "scope": "all-field finite rank-sum exclusion through 52; no new factoring exponent",
        "previous_cutoff": 50, "last_certified_dimension": 52,
        "dimensions_tested_here": [51, 52],
        "no_claim_about_dimension_53": True,
        "arithmetic": "exact integers and fractions; no numerical optimizer",
        "source_sha256": hashes,
        "scalar_cases": sum(row["scalar_cases"] for row in dimensions),
        "ratio_knot_tests": sum(row["knots"] for row in certificate_rows),
        "witnesses": certificate_rows,
        "dimensions": dimensions,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"], "last_certified_dimension": 52,
        "scalar_cases": result["scalar_cases"], "ratio_knot_tests": result["ratio_knot_tests"],
        "dimensions": dimensions, "output": str(args.output) if args.output else None,
    }, indent=2))


if __name__ == "__main__":
    main()
