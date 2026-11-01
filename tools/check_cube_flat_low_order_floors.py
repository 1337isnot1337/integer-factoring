#!/usr/bin/env python3
"""Exact low-order entropy obstructions for the common-degree cube criterion."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / "evidence" / "cube-flat-low-order-floors.json"


def low_order_obstructions(profiles, degree: int) -> dict:
    """Classify already certified generic profiles; do not use sampled minima.

    This function verifies elementary consistency bounds, not genericity.
    Actual polynomial tensor certificates remain the caller's responsibility.
    """
    if (len(profiles) != 3 or any(len(p) != 16 for p in profiles)
            or any(list(p) != sorted(p) for p in profiles)
            or any(not isinstance(x, int) or x < 0 for p in profiles for x in p)
            or not isinstance(degree, int) or degree < 0):
        raise ValueError("Expected three sorted generic integer profiles and a degree.")
    common = [p[0] for p in profiles]
    residual_degree = degree - sum(common)
    normalized = [[x - p[0] for x in p] for p in profiles]
    if residual_degree < 0 or any(sum(p) > 16 * residual_degree for p in normalized):
        raise ValueError("The profiles violate the normalized polynomial degree bound.")
    if residual_degree == 0:
        assert all(not any(p) for p in normalized)
        return {"common_coordinate_factor": common, "residual_degree": 0,
                "constant_endpoint": True, "pressure": "1"}
    rank_axes = [a for a, p in enumerate(normalized) if p[1] == 0]
    second_sum = sum(p[1] for p in normalized)
    return {"common_coordinate_factor": common, "residual_degree": residual_degree,
            "face_rank_at_least_two_axes": rank_axes,
            "normalized_second_order_sum": second_sum,
            "small_average_second_order": second_sum <= residual_degree,
            "certified_floor": "rho > 29/50 > 1/sqrt(3)"
            if rank_axes or second_sum <= residual_degree else None}


def comparison(name: str, left: int, right: int) -> dict:
    assert left > right
    difference = left - right
    return {"name": name, "strict_integer_comparison": True,
            "left_bits": left.bit_length(), "right_bits": right.bit_length(),
            "positive_difference_bits": difference.bit_length(),
            "difference_sha256": hashlib.sha256(difference.to_bytes(
                (difference.bit_length() + 7) // 8, "big")).hexdigest()}


def run() -> dict:
    f = Fraction
    assert 50 + 15 * 2 == 80
    assert 2 * 33 + 14 == 80
    assert (2 * f(2, 80) + f(1, 80)) * 16 / 3 == f(1, 3)
    assert 33 + 12 + 14 == 59
    assert f(16, 59) + f(11, 59) / 3 == f(1, 3)
    assert 3 * 29**2 > 50**2
    certificates = [
        comparison("a_face_of_rank_at_least_two", 80**1200,
                   2**2784 * (50**50 * 2**30)**10 * (33**66)**5),
        comparison("small_average_second_order", 59**1475,
                   2**3422 * (33**33 * 12**12)**25),
    ]
    even_seed = [0, 2, 2, 4, 4, 4, 6, 6, 6, 6, 8, 8, 8, 10, 10, 12]
    odd_seed = [0, 2, 2, 2, 2, 2, 2, 6, 6, 10, 10, 10, 10, 10, 10, 12]
    controls = {
        "even_seed": low_order_obstructions([even_seed] * 3, 6),
        "odd_seed": low_order_obstructions([odd_seed] * 3, 6),
        "one_rank_two_face": low_order_obstructions(
            [[0, 4, 4, 4] + [7] * 12] * 2 + [[0, 0] + [6] * 14], 6),
        "common_scalar_factor": low_order_obstructions(
            [[x + 2 for x in even_seed], [x + 1 for x in even_seed], even_seed], 9),
        "unexcluded_relaxation": low_order_obstructions([[0, 4, 4, 4] + [7] * 12] * 3, 6),
        "constant_endpoint": low_order_obstructions([[2] * 16, [3] * 16, [1] * 16], 6),
    }
    assert controls["even_seed"]["small_average_second_order"]
    assert controls["odd_seed"]["small_average_second_order"]
    assert controls["one_rank_two_face"]["face_rank_at_least_two_axes"] == [2]
    assert controls["common_scalar_factor"]["residual_degree"] == 6
    assert controls["unexcluded_relaxation"]["certified_floor"] is None
    assert controls["constant_endpoint"]["constant_endpoint"]
    try:
        low_order_obstructions([[5] * 16] * 3, 6)
    except ValueError:
        pass
    else:
        raise AssertionError("Impossible common coordinate degree was accepted.")
    return {"scope": "Certified generic profiles of arbitrary homogeneous flat binary cube tensors; odd coefficients allowed.",
            "proved_floor_under_either_condition": "rho > 29/50 > 1/sqrt(3)",
            "face_rank_weights": {"denominator": 80, "A_numerators": [50] + [2] * 15,
                                  "B_numerators": [33, 33] + [1] * 14},
            "second_order_weights": {"denominator": 59, "numerators": [33, 12] + [1] * 14},
            "entropy_certificates": certificates, "interface_controls": controls,
            "floating_point_used_for_certificate": False,
            "generic_smith_profiles_computed_by_this_checker": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-evidence", action="store_true")
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    result = run()
    if args.write_evidence:
        args.evidence.parent.mkdir(parents=True, exist_ok=True)
        args.evidence.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif args.evidence.exists():
        assert json.loads(args.evidence.read_text()) == result
    print("PASS: two exact flat low-order obstructions >29/50>1/sqrt(3)")


if __name__ == "__main__":
    main()
