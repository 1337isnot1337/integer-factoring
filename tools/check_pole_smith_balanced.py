#!/usr/bin/env python3
"""Exact controls for the overlap-count formula for balanced pole profiles."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random

from check_pole_smith_stabilization import matrix_at, smith_orders
from pole_cauchy_conditioning import balanced_padding
from pole_smith_balanced import balanced_smith_histogram, dyadic_balanced_padding
from pole_smith_flow import smith_profile


def expanded_orders(result: dict) -> list[int]:
    return [item["order"] for item in result["histogram"] for _ in range(item["multiplicity"])]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(811)
    for index in range(1000):
        n = rng.randrange(51, 1000, 2)
        left = (rng.randrange(1, 14, 2), rng.randrange(1, 9))
        right = (rng.randrange(1, 14, 2), rng.randrange(1, 9))
        aa, bb, _ = balanced_padding(n, left, right)
        rows = tuple(a + 2 * j for a, count in aa for j in range(count))
        cols = tuple(a + 2 * j for a, count in bb for j in range(count))
        balanced = balanced_smith_histogram(n, rows, cols)
        orders = expanded_orders(balanced)
        assert orders == smith_profile(n, rows, cols)["orders"]
        overlap = balanced["overlap_counts"]
        assert orders == [sum(count < j for count in overlap) for j in range(1, len(rows) + 1)]
        if index < 100:
            assert orders == smith_orders(matrix_at(n, rows, cols, balanced["k"]))
        if index < 200:
            aa, bb, modulus = dyadic_balanced_padding(n, left, right)
            rows = tuple(a + 2 * j for a, count in aa for j in range(count))
            cols = tuple(a + 2 * j for a, count in bb for j in range(count))
            by_residue = {b % modulus: b for b in cols}
            assert len(by_residue) == len(cols)
            dyadic = balanced_smith_histogram(n, rows, cols)
            k = dyadic["k"]
            diagonal_orders = []
            for a in rows:
                b = by_residue[n * pow(a, -1, modulus) % modulus]
                difference = abs(n - a * b)
                valuation = k if difference == 0 else (difference & -difference).bit_length() - 1
                diagonal_orders.append(k - min(k, valuation))
            assert expanded_orders(dyadic) == sorted(diagonal_orders)
            assert max(diagonal_orders) <= k - (len(rows).bit_length() - 1) - 1
            if index < 100:
                assert sorted(diagonal_orders) == smith_orders(matrix_at(n, rows, cols, k))
    try:
        balanced_smith_histogram(101, (1, 3, 5, 7), (1, 9, 17, 25))
    except ValueError as error:
        assert "imbalance" in str(error)
    else:
        raise AssertionError("The imbalance hypothesis must be enforced")
    huge_k = 10 ** 80
    shifted = balanced_smith_histogram(77, (3, 5, 7, 9), (7, 9, 11, 13), huge_k)
    assert expanded_orders(shifted) == smith_profile(77, (3, 5, 7, 9), (7, 9, 11, 13), huge_k)["orders"]
    first, second = (1 << 100) + 3, (1 << 100) + 7
    size = 65536
    large = balanced_smith_histogram(first * second,
                                     tuple(first + 2 * i for i in range(size)),
                                     tuple(second + 2 * i for i in range(size)))
    digest = hashlib.sha256(json.dumps(expanded_orders(large)).encode()).hexdigest()
    assert digest == "5119a2e4c7efcd01771df1906c0608e5a40d88465147b05081747ee83d4a6b5f"
    result = {"claim": "Balanced pole Smith orders are the conjugate overlap-count partition",
              "full_profile_comparisons": 1000, "independent_dense_Qw_comparisons": 200,
              "dyadic_diagonal_comparisons": 200, "huge_displacement_control": shifted,
              "scope_failure_rejected": True, "large_control": {**large, "profile_sha256": digest},
              "scope": "One pass over explicit node-tree; implicit overlap counts remain unresolved."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
