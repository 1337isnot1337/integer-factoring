#!/usr/bin/env python3
"""Exact all-field pressure obstruction for the height-16 boundary family.

The integer certificate proves Z(s**17)**12 > A33(s**24)**17 on 0<s<=1.
Here Z is the complete natural Smith counting polynomial of the fixed
28-coordinate lower space.  Its submatrix protection gives rho >= 17/24
for every homogeneous retraction on the five-dimensional height-16 cell.

With --recompute-profiles, reconstruct both complete lower-space profiles
using the original binary constructor and the field-valued constructor.
The default polynomial certificate itself needs only the standard library.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Z = (
    1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6, 7,
    8, 8, 9, 11, 12, 12, 14, 16, 17, 18, 20, 21, 23, 25, 27, 26, 30,
    31, 34, 32, 34, 38, 38, 33, 29, 32, 21, 14, 18, 11, 12, 8, 8, 4,
    3, 7, 3, 0, 2, 10, 1, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 2,
)


def height(word: int) -> int:
    return sum(2 * bit + 1 for bit in range(word.bit_length()) if word >> bit & 1)


LOWER = tuple(sorted((word for word in range(256) if height(word) < 16),
                     key=lambda word: (height(word), word)))
BOUNDARY = tuple(word for word in range(256) if height(word) == 16)
HEIGHT_COUNTS = Counter(map(height, LOWER + BOUNDARY))
A33 = tuple(HEIGHT_COUNTS[i] for i in range(17))
DENOMINATOR_DIFFERENCE = 17 * (len(Z) - 1) * 12 - 24 * 16 * 17


def homogeneous_numerator(polynomial: tuple[int, ...], numerator: int,
                          denominator: int, power: int) -> int:
    """Return denominator**(power*degree) * P((numerator/denominator)**power)."""
    x, y = numerator**power, denominator**power
    result, multiple = polynomial[-1], 1
    for coefficient in reversed(polynomial[:-1]):
        multiple *= y
        result = result * x + coefficient * multiple
    return result


@lru_cache(None)
def denominator_power(denominator: int) -> int:
    return denominator**abs(DENOMINATOR_DIFFERENCE)


def interval_positive(depth: int, index: int) -> bool:
    # The interval is [(4*2**depth+index)/(5*2**depth),
    #                  (4*2**depth+index+1)/(5*2**depth)].
    denominator = 5 << depth
    lower = (4 << depth) + index
    upper = lower + 1
    lhs = homogeneous_numerator(Z, lower, denominator, 17)**12
    rhs = homogeneous_numerator(A33, upper, denominator, 24)**17
    if DENOMINATOR_DIFFERENCE >= 0:
        rhs *= denominator_power(denominator)
    else:
        lhs *= denominator_power(denominator)
    return lhs > rhs


def exact_cover() -> list[list[int]]:
    leaves, stack = [], [(0, 0)]
    while stack:
        depth, index = stack.pop()
        if interval_positive(depth, index):
            leaves.append([depth, index])
        else:
            if depth >= 24:
                raise AssertionError("Protected-floor cover did not close at depth 24")
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


def profile_controls() -> dict:
    from coherent_hot_compounds import profiles as original_profiles
    from coherent_f4_boundary import profiles as field_profiles
    binary_maps = [1 << word for word in LOWER]
    original = original_profiles(binary_maps, binary_maps, 8)
    field_maps = [{word: 1} for word in LOWER]
    extension = field_profiles(field_maps, field_maps, 8)
    expected = [i for i, count in enumerate(Z) for _ in range(count)]
    for result in (original, extension):
        assert result["heights"] == list(map(height, LOWER))
        assert result["zero_orders_t"] == result["infinity_orders_t"] == expected
    dependencies = (
        "coherent_hot_compounds.py", "coherent_hot_compounds.cpp",
        "coherent_f4_boundary.py", "coherent_f4_boundary.cpp",
        "coherent_hot_compounds_fast.cpp", "selected_gaussian_coefficients.cpp",
        "native_rectangular_smith_f4.cpp",
    )
    return {
        "status": "passed",
        "constructors": ["original binary compound", "field-valued normalized compound"],
        "checked_complete_orders": 4 * 28**2,
        "profile_sha256": sha256(json.dumps(expected, separators=(",", ":")).encode()).hexdigest(),
        "occupation_blocks": [
            {name: block[name] for name in ("charge", "size", "orders_t", "precision_t",
                                           "coefficient_degree_bound_t", "precision_attempts")}
            for block in original["cuts"][0]["blocks"]
        ],
        "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                           for name in dependencies},
    }


def certificate() -> dict:
    assert len(LOWER) == 28 and BOUNDARY == (15, 24, 36, 66, 129)
    assert len(Z) == 83 and sum(Z) == 28**2 and Z[0] == Z[1] == 1
    assert A33[0] == A33[1] == 1 and A33[2] == 0 and sum(A33) == 33
    assert DENOMINATOR_DIFFERENCE == 10200

    # For 0<s<=4/5, t=s**24 satisfies A33(t)<=1+t+31*t**3<=1+2*t.
    assert 31 * 4**48 < 5**48
    # Bound the binomial tail by a geometric series, then use Z(s**17)>=1+s**17.
    assert 32 * 4**24 < 5**24
    for j in range(1, 18):
        assert comb(17, j) * 2**j <= 34 * 32**(j - 1)
    small_gap = 12 * (5**24 - 32 * 4**24) - 34 * 4**7 * 5**17
    assert small_gap > 0

    leaves = exact_cover()
    verify_cover(leaves)
    # The required strict threshold separation is (17/24)**2 > 1/2.
    threshold_gap = 2 * 17**2 - 24**2
    assert threshold_gap == 2
    return {
        "schema": 1,
        "status": "passed",
        "scope": "Natural common-scale two-cut profile pressure for arbitrary homogeneous retractions on height 16 retaining all 28 lower coordinates, over every characteristic-two field",
        "factoring_exponent_improved": False,
        "pressure_infimum_lower_bound": [17, 24],
        "threshold_integer_gap": threshold_gap,
        "protected_words": list(LOWER),
        "boundary_words": list(BOUNDARY),
        "protected_smith_counts": list(Z),
        "largest_physical_height_counts": list(A33),
        "polynomial_inequality": {"left": "Z(s^17)^12", "right": "A33(s^24)^17"},
        "small_interval": {"upper": [4, 5], "strict_gap_numerator": small_gap},
        "cover_interval": [[4, 5], [1, 1]],
        "cover_leaves": len(leaves),
        "maximum_depth": max(depth for depth, _ in leaves),
        "denominator_power_difference": DENOMINATOR_DIFFERENCE,
        "cover_sha256": sha256(json.dumps(leaves, separators=(",", ":")).encode()).hexdigest(),
        "dyadic_cover": leaves,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--recompute-profiles", action="store_true")
    args = parser.parse_args()
    result = certificate()
    if args.recompute_profiles:
        result["profile_controls"] = profile_controls()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    summary = {key: value for key, value in result.items()
               if key not in ("dyadic_cover", "protected_smith_counts", "protected_words", "profile_controls")}
    if args.recompute_profiles:
        summary["recomputed_complete_orders"] = result["profile_controls"]["checked_complete_orders"]
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
