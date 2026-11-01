#!/usr/bin/env python3
"""Certify the explicit eight-letter nonorthogonal reflection crossing."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from check_block_smith_crossings import toeplitz_rank
from eight_letter_reflection_retraction import (
    CODE, HEIGHTS, displayed_maps, homogeneous_restriction,
)


ROOT = Path(__file__).resolve().parents[1]


def counting_polynomial(orders):
    counts = Counter(orders)
    return [counts[index] for index in range(max(orders) + 1)]


def multiply(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def integer_hash(value):
    return sha256(value.to_bytes((value.bit_length() + 7) // 8, "big")).hexdigest()


def comparison(coordinate, mixed):
    assert coordinate["zero_orders"] == coordinate["infinity_orders"]
    c = counting_polynomial(coordinate["zero_orders"])
    z = counting_polynomial(mixed["zero_orders"])
    i = counting_polynomial(mixed["infinity_orders"])
    left, right = multiply(c, c), multiply(z, i)
    difference = [(left[j] if j < len(left) else 0)
                  - (right[j] if j < len(right) else 0)
                  for j in range(max(len(left), len(right)))]
    assert all(value == 0 for value in difference[:16])
    running, quotient = 0, []
    for coefficient in difference[16:-1]:
        running += coefficient
        quotient.append(running)
    assert running + difference[-1] == 0
    assert [0] * 16 + multiply(quotient, [1, -1]) == difference
    assert quotient[0] == 1
    first_negative = next(j for j, value in enumerate(quotient) if value < 0)
    negative_mass = -sum(value for value in quotient if value < 0)
    assert (first_negative, negative_mass) == (23, 749)
    assert negative_mass * 2**first_negative < 3**first_negative
    return {
        "coordinate_counting_polynomial": c,
        "mixed_zero_counting_polynomial": z,
        "mixed_infinity_counting_polynomial": i,
        "coordinate_pressure_minus_mixed_pressure": difference,
        "factorization": "r^16 * (1-r) * Q(r)",
        "Q_coefficients": quotient,
        "Q_constant_coefficient": 1,
        "Q_first_negative_order": first_negative,
        "Q_total_negative_mass": negative_mass,
        "Q_lower_bound_on_zero_to_two_thirds": str(Fraction(
            3**first_negative - negative_mass * 2**first_negative,
            3**first_negative)),
        "strict_pressure_gain_for_every_finite_scale_at_least": "3/2",
    }


def rate_certificate(coordinate, mixed):
    # s=(73/71)^18 and rho=18/25 make every type weight rational.
    a, b = 73, 71

    def numerator(orders, power):
        maximum = max(orders)
        return (sum(b**(power * order) * a**(power * (maximum - order))
                    for order in orders), power * maximum)

    zn, ze = numerator(mixed["zero_orders"], 18)
    inn, ie = numerator(mixed["infinity_orders"], 18)
    cn, ce = numerator(coordinate["zero_orders"], 18)
    an, ae = numerator(HEIGHTS, 25)
    left = (zn * inn)**25 * a**(72 * ae)
    right = an**72 * a**(25 * (ze + ie))
    coordinate_left = cn**50 * a**(72 * ae)
    coordinate_right = an**72 * a**(50 * ce)
    assert left < right
    assert coordinate_left > coordinate_right
    type_numerators = [b**(25 * height) * a**(25 * (max(HEIGHTS) - height))
                       for height in HEIGHTS]
    assert sum(type_numerators) == an
    return {
        "scale": "(73/71)^18",
        "rho_strict_upper_bound": "18/25",
        "counting_polynomial_variable": "r=(71/73)^18",
        "type_weight_formula": "p_j = (71/73)^(25*h_j) / sum_k (71/73)^(25*h_k)",
        "type_numerators": [str(value) for value in type_numerators],
        "type_denominator": str(an),
        "integer_inequality": "(Znum*Inum)^25 * 73^(72*Aexp) < Anum^72 * 73^(25*(Zexp+Iexp))",
        "numerators": {"Z": str(zn), "I": str(inn), "A": str(an)},
        "denominator_exponents": {"Z": ze, "I": ie, "A": ae},
        "strict_integer_comparison": True,
        "comparison_bit_lengths": [left.bit_length(), right.bit_length()],
        "comparison_integer_hashes": [integer_hash(left), integer_hash(right)],
        "coordinate_fails_same_trial_rate_at_same_scale": True,
        "coordinate_comparison_bit_lengths": [coordinate_left.bit_length(),
                                               coordinate_right.bit_length()],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records, controls = [], 0
    for nonorthogonal in (False, True):
        record = homogeneous_restriction(*displayed_maps(nonorthogonal), keep_matrices=True)
        matrices = record.pop("matrices")
        assert record["actual_polynomial_degree"] == 32
        assert record["heights"] == HEIGHTS
        for matrix, orders in zip(matrices, (record["zero_orders"], record["infinity_orders"])):
            for precision in (1, 2, 4, 8, 16, 24, 33):
                assert toeplitz_rank(matrix, precision) == sum(
                    max(0, precision - order) for order in orders)
                controls += 1
        record["nonorthogonal"] = nonorthogonal
        records.append(record)
    coordinate, mixed = records
    rejection_controls = 0
    for bad_maps in ((displayed_maps()[0], displayed_maps()[0][::-1], HEIGHTS),
                     (displayed_maps()[0], tuple(value ^ (1 if j == 4 else 0)
                       for j, value in enumerate(displayed_maps()[1])), HEIGHTS)):
        try:
            homogeneous_restriction(*bad_maps)
        except ValueError:
            rejection_controls += 1
        else:
            raise AssertionError("A malformed physical retraction was accepted.")
    result = {
        "status": "pass",
        "scope": "One explicit eight-letter homogeneous nonorthogonal M=8 reflection restriction, over every field of characteristic two",
        "retained_coordinate_labels": CODE,
        "ambient_endpoint_entries": 65536,
        "kernel_endpoint_duality_entries": 256,
        "compressed_endpoint_entries": sum(record["literal_endpoint_entries"] for record in records),
        "compressed_infinity_restoration_entries": sum(record["literal_infinity_restoration_entries"] for record in records),
        "complete_smith_orders": 256,
        "independent_toeplitz_rank_controls": controls,
        "malformed_retraction_rejections": rejection_controls,
        "records": records,
        "pointwise_pressure_comparison": comparison(coordinate, mixed),
        "rational_type_certificate": rate_certificate(coordinate, mixed),
        "engine_sha256": sha256((ROOT / "tools/eight_letter_reflection_retraction.py").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "floating_arithmetic_used_for_certification": False,
        "global_optimality_claimed": False,
        "factoring_exponent_improved": False,
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("records", "pointwise_pressure_comparison", "rational_type_certificate")}, indent=2))


if __name__ == "__main__":
    main()
