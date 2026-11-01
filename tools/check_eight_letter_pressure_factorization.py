#!/usr/bin/env python3
"""Certify the wider exact pressure-gain interval for the q=8 retraction."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from eight_letter_reflection_retraction import displayed_maps, homogeneous_restriction
from check_eight_letter_reflection_retraction import counting_polynomial, multiply


def add(*terms):
    output = [0] * max(map(len, terms))
    for term in terms:
        for j, value in enumerate(term):
            output[j] += value
    while len(output) > 1 and output[-1] == 0:
        output.pop()
    return output


def negative(values):
    return [-value for value in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    coordinate, mixed = [homogeneous_restriction(*displayed_maps(flag))
                         for flag in (False, True)]
    assert coordinate["zero_orders"] == coordinate["infinity_orders"]
    c = counting_polynomial(coordinate["zero_orders"])
    z = counting_polynomial(mixed["zero_orders"])
    i = counting_polynomial(mixed["infinity_orders"])
    K = [0] * 16 + multiply([1, 1], [1, 0, 0, 0, -1])
    A = multiply([1, 0, 0, 1], [1, 0, 0, 0, 1])
    B = [1, 1, 0, 1, 1, 0, 0, 1]
    r4B = [0] * 4 + B
    assert z == add(c, negative(multiply(K, A)))
    assert i == add(c, multiply(K, r4B))
    difference = add(multiply(c, c), negative(multiply(z, i)))
    bracket = add(A, negative(r4B))
    assert bracket == [1, 0, 0, 1, 0, -1, 0, 0, -1, 0, 0, -1]
    decomposed = add(multiply(multiply(K, c), bracket),
                     multiply(multiply(K, K), multiply(A, r4B)))
    assert difference == decomposed
    numerator = 10**11 - 9**8 * 10**3 - 9**11
    assert numerator > 0
    root = Path(__file__).resolve().parents[1]
    result = {
        "status": "pass",
        "scope": "Same characteristic-two q8 crossing and original physical type charges",
        "strict_pressure_gain_for_every_finite_scale_at_least": "10/9",
        "bracket_lower_bound_for_zero_to_nine_tenths": str(Fraction(numerator, 10**11)),
        "four_complete_profiles_recomputed": True,
        "profile_difference_identities_checked": 2,
        "pressure_identity_coefficient_count": len(difference),
        "K_coefficients": K, "A_coefficients": A, "B_coefficients": B,
        "bracket_coefficients": bracket,
        "comparison_numerator": numerator,
        "dependency_sha256": {name: sha256((root / name).read_bytes()).hexdigest() for name in (
            "tools/eight_letter_reflection_retraction.py",
            "tools/check_eight_letter_reflection_retraction.py")},
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
