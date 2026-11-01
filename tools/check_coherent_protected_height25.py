#!/usr/bin/env python3
"""Integer pressure certificate for all protected lower-28 spaces through height 25.

The same 784-order profile used by the height-16 boundary obstruction
dominates the physical entropy of the entire 104-letter height-25 ball.
This checker proves Z(s**17)**12 > A104(s**24)**17 on 0<s<=1.
The stronger written theorem protects literal polynomial submatrices and
does not require deletion closure of the larger physical map spaces.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from check_coherent_boundary_protected_floor import (
    LOWER, Z, height, homogeneous_numerator, profile_controls,
)


ROOT = Path(__file__).resolve().parents[1]
CAP = 25


def physical_height_counts() -> tuple[int, ...]:
    counts = [1] + [0] * CAP
    for odd in range(1, CAP + 1, 2):
        for degree in range(CAP, odd - 1, -1):
            counts[degree] += counts[degree - odd]
    direct = Counter(height(word) for word in range(1 << 13) if height(word) <= CAP)
    assert counts == [direct[i] for i in range(CAP + 1)]
    assert sum(counts) == 104
    return tuple(counts)


A104 = physical_height_counts()
DENOMINATOR_DIFFERENCE = 17 * (len(Z) - 1) * 12 - 24 * CAP * 17


@lru_cache(None)
def denominator_power(denominator: int) -> int:
    return denominator**DENOMINATOR_DIFFERENCE


def interval_positive(depth: int, index: int) -> bool:
    denominator = 5 << depth
    lower = (4 << depth) + index
    upper = lower + 1
    lhs = homogeneous_numerator(Z, lower, denominator, 17)**12
    rhs = homogeneous_numerator(A104, upper, denominator, 24)**17
    return lhs > rhs * denominator_power(denominator)


def exact_cover() -> list[list[int]]:
    leaves, stack = [], [(0, 0)]
    while stack:
        depth, index = stack.pop()
        if interval_positive(depth, index):
            leaves.append([depth, index])
        else:
            if depth >= 25:
                raise AssertionError("Height-25 protected cover did not close at depth 25")
            stack.extend(((depth + 1, 2 * index + 1), (depth + 1, 2 * index)))
    return leaves


def verify_cover(leaves: list[list[int]]) -> None:
    maximum = max(depth for depth, _ in leaves)
    previous = 0
    for depth, index in leaves:
        assert 0 <= index < 1 << depth
        assert index << (maximum - depth) == previous
        assert interval_positive(depth, index)
        previous = (index + 1) << (maximum - depth)
    assert previous == 1 << maximum


def certificate() -> dict:
    assert len(Z) == 83 and sum(Z) == 28**2 and Z[0] == Z[1] == 1
    assert A104[0] == A104[1] == 1 and A104[2] == 0
    assert DENOMINATOR_DIFFERENCE == 6528
    # For t=s**24 and 0<s<=4/5, A104(t)<=1+t+102*t**3<=1+2*t.
    assert 102 * 4**48 < 5**48
    assert 32 * 4**24 < 5**24
    assert all(comb(17, j) * 2**j <= 34 * 32**(j - 1) for j in range(1, 18))
    small_gap = 12 * (5**24 - 32 * 4**24) - 34 * 4**7 * 5**17
    assert small_gap == 182169346247795596
    assert 2 * 17**2 - 24**2 == 2
    leaves = exact_cover()
    verify_cover(leaves)
    return {
        "schema": 1,
        "status": "passed",
        "scope": "Natural common-scale two-cut pressure for height-homogeneous constant split maps retaining all 28 words below height 16 and supported through height 25; arbitrary characteristic-two field, separate paired alphabets and types, no deletion-closure hypothesis on the larger spaces",
        "factoring_exponent_improved": False,
        "pressure_infimum_lower_bound": [17, 24],
        "threshold_integer_gap": 2,
        "protected_words": list(LOWER),
        "protected_smith_counts": list(Z),
        "maximum_height": CAP,
        "largest_alphabet_dimension": sum(A104),
        "largest_physical_height_counts": list(A104),
        "polynomial_inequality": {"left": "Z(s^17)^12", "right": "A104(s^24)^17"},
        "small_interval": {"upper": [4, 5], "strict_gap_numerator": small_gap,
                           "tail_coefficient_bound": 102},
        "cover_interval": [[4, 5], [1, 1]],
        "cover_leaves": len(leaves),
        "maximum_depth": max(depth for depth, _ in leaves),
        "denominator_power_difference": DENOMINATOR_DIFFERENCE,
        "cover_sha256": sha256(json.dumps(leaves, separators=(",", ":")).encode()).hexdigest(),
        "dyadic_cover": leaves,
        "predecessor_checker_sha256": sha256(
            (ROOT / "tools" / "check_coherent_boundary_protected_floor.py").read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--recompute-profiles", action="store_true")
    args = parser.parse_args()
    result = certificate()
    if args.recompute_profiles:
        result["protected_profile_controls"] = profile_controls()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    summary = {key: value for key, value in result.items()
               if key not in ("dyadic_cover", "protected_smith_counts", "protected_words",
                              "protected_profile_controls")}
    if args.recompute_profiles:
        summary["recomputed_complete_orders"] = result["protected_profile_controls"]["checked_complete_orders"]
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
