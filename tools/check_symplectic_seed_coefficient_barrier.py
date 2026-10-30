#!/usr/bin/env python3
"""Exact character-minor bounds for common-parameter radical-seed powers."""

import argparse
import hashlib
import json
from math import comb, exp, isqrt, log
from pathlib import Path


def centered_diagonal_lower(power):
    row = [1]
    for _ in range(power):
        successor = [0] * (len(row) + 3)
        for degree, value in enumerate(row):
            for shift, multiplicity in ((1, 40), (2, 436), (3, 200)):
                successor[degree + shift] += value * multiplicity
        row = successor
    column = [0] * (2 * power + 1)
    polynomial = 1
    for j in range(power + 1):
        multiplicity = comb(power, j) * 668 ** (power - j) * 8 ** j
        for degree in range(2 * j + 1):
            if (polynomial >> degree) & 1:
                column[degree] += multiplicity
        polynomial ^= (polynomial << 1) ^ (polynomial << 2)
    # Integer square roots round downward, so this is a rigorous lower bound.
    lower = sum(isqrt(row[k] * column[k]) for k in range(2 * power + 1))
    return {"power": power, "rank_cost_integer_lower": lower,
            "nth_root_for_display": exp(log(lower) / power),
            "exceeds_seed_threshold_by_integer_comparison": lower ** 70 > 26 ** (99 * power)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", type=Path,
                        default=Path("evidence/symplectic-radical-seed-algebra.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    algebra = json.loads(args.algebra.read_text())
    row, column = algebra["cuts"]
    assert row["cut"] == "row" and column["cut"] == "column"
    assert row["characters_on_E_A_B_C_D_F"] == [
        [1, 0, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0], [1, 1, 1, 1, 0, 0]]
    assert row["representation_radical_filtration"][0]["character_multiplicities"] == [436, 40, 200]
    assert column["characters_on_E_A_B_C_D_F"] == [
        [1, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0]]
    assert column["representation_radical_filtration"][0]["character_multiplicities"] == [668, 8]
    # The rational character type (6/10,4/10,0) matches column type 7/10.
    numerator = 10 ** 20 * 40 ** 6 * 436 ** 4 * 668 ** 3 * 8 ** 7
    denominator = 6 ** 6 * 4 ** 4 * 7 ** 7 * 3 ** 3
    assert 99 ** 2 > 2 * 70 ** 2
    assert numerator ** 70 > denominator ** 70 * 26 ** 1980
    # For every power divisible by ten and at least 80, the elementary
    # binomial-type bound B^n/(n+1) already exceeds 26^(n*99/70).
    assert numerator ** 280 > denominator ** 280 * 81 ** 70 * 26 ** 7920
    assert numerator ** 35 * 81 ** 70 > denominator ** 35 * 91 ** 70 * 26 ** 990
    result = {
        "scope": "literal polynomial coefficient decompositions; arbitrary recombinations are not excluded",
        "labels": 26,
        "row_character_factors": [[1, 40], [2, 436], [3, 200]],
        "column_character_multiplicities": [668, 8],
        "uniform_polynomial_curve_lower": {
            "family": "R(t(x),z(x)) with nonconstant polynomial t and arbitrary polynomial z",
            "twentieth_power_numerator": numerator,
            "twentieth_power_denominator": denominator,
            "base_for_display": exp((log(numerator) - log(denominator)) / 20),
            "exact_integer_threshold_comparison_passed": True,
            "finite_bound_for_powers_divisible_by_ten_at_least": 80,
        },
        "centered_power_checks": [centered_diagonal_lower(n) for n in (1, 2, 4, 8, 16, 32, 64, 128)],
        "algebra_sha256": hashlib.sha256(args.algebra.read_bytes()).hexdigest(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
