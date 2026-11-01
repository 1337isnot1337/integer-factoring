#!/usr/bin/env python3
"""A protected rectangular subtensor blocks one-letter auxiliary deformations.

The exact all-scale bound is rho >= 8/11 for the specified reflection
gauges and the full-word regular-channel pressure interface.  It does not
bound actual ranks after direct type restriction or exclude other gauges.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path

import adaptive_smith_lattices as adaptive
import check_cancellation_module_automata as arithmetic
import check_triangular_auxiliary_row_fiber as reflection


ROOT = Path(__file__).resolve().parent.parent
BASE_PROFILE = (0, 1, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11)
LETTER_PROFILES = {
    1: (0, 1, 3, 4, 5, 6, 7, 8, 8, 9, 12, 13),
    2: (0, 1, 3, 4, 5, 7, 8, 8, 10, 11, 12, 15),
    3: BASE_PROFILE,
}


def polynomial_product(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    result[i + j] += a * b
    return result


def polynomial_power(value, exponent):
    result = [1]
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, value)
        exponent //= 2
        if exponent:
            value = polynomial_product(value, value)
    return result


def profile_polynomial():
    result = [0] * 12
    for order in BASE_PROFILE:
        result[order] += 1
    factors = ([1, 1], [1, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0, 1])
    expected = polynomial_product(polynomial_product(factors[0], factors[1]), factors[2])
    assert result == expected
    return result


def bernstein_numerators(polynomial):
    """Common positive denominator n! for degree-n Bernstein coefficients.

    If p(t)=sum_j c_j t^j=sum_k b_k C(n,k)t^k(1-t)^(n-k),
    then n! b_k=sum_{j<=k} c_j C(k,j) j! (n-j)!.
    """
    degree = len(polynomial) - 1
    factorials = [factorial(i) for i in range(degree + 1)]
    terms = [(j, value * factorials[j] * factorials[degree - j])
             for j, value in enumerate(polynomial) if value]
    numerators = [sum(value * comb(k, j) for j, value in terms if j <= k)
                  for k in range(degree + 1)]
    return numerators, factorials[-1]


def all_scale_certificate():
    # g(t^8)^11 - (1+t^11)^16 (1+t^33)^16 = t^8 p(t).
    g = profile_polynomial()
    substituted = [0] * (8 * (len(g) - 1) + 1)
    for degree, coefficient in enumerate(g):
        substituted[8 * degree] = coefficient
    left = polynomial_power(substituted, 11)
    first = [int(i in (0, 11)) for i in range(12)]
    second = [int(i in (0, 33)) for i in range(34)]
    right = polynomial_product(polynomial_power(first, 16), polynomial_power(second, 16))
    difference = [value - (right[i] if i < len(right) else 0)
                  for i, value in enumerate(left)]
    assert difference[:8] == [0] * 8
    polynomial = difference[8:]
    assert len(polynomial) - 1 == 960
    numerators, denominator = bernstein_numerators(polynomial)
    assert all(value > 8 * denominator for value in numerators)
    minimum_index = min(range(len(numerators)), key=lambda i: numerators[i])
    assert minimum_index == 628
    assert 8 * denominator < numerators[minimum_index] < 9 * denominator
    # The rational lower bound is strictly beyond the target 1/sqrt(2).
    assert 2 * 8 * 8 > 11 * 11
    canonical = json.dumps(polynomial, separators=(",", ":"))
    return {
        "rho_lower_bound": "8/11",
        "strictly_above_one_over_sqrt_two": "128 > 121",
        "substitution": "r = t^8",
        "polynomial_identity": "g(t^8)^11 - (1+t^11)^16(1+t^33)^16 = t^8 p(t)",
        "degree_of_p": 960,
        "bernstein_coefficients_checked": len(numerators),
        "common_denominator": "960!",
        "every_bernstein_coefficient_strictly_greater_than": 8,
        "minimum_coefficient_index": minimum_index,
        "minimum_coefficient_strictly_less_than": 9,
        "polynomial_coefficients_sha256": sha256(canonical.encode()).hexdigest(),
        "floating_arithmetic_used": False,
    }


def additional_letter_certificate(letter):
    orders = LETTER_PROFILES[letter]
    if letter == 1:
        substitution, left_power, right_power, first_power, second_power = 8, 11, 16, 11, 33
        rate, floor, expected_degree = "8/11", 8, 1136
    else:
        assert letter == 2
        substitution, left_power, right_power, first_power, second_power = 5, 7, 10, 7, 21
        rate, floor, expected_degree = "5/7", 4, 520
    g = [0] * (substitution * max(orders) + 1)
    for order in orders:
        g[substitution * order] += 1
    left = polynomial_power(g, left_power)
    right = polynomial_product(
        polynomial_power([int(i in (0, first_power)) for i in range(first_power + 1)], right_power),
        polynomial_power([int(i in (0, second_power)) for i in range(second_power + 1)], right_power),
    )
    difference = [value - (right[i] if i < len(right) else 0)
                  for i, value in enumerate(left)]
    assert difference[:substitution] == [0] * substitution
    polynomial = difference[substitution:]
    assert len(polynomial) - 1 == expected_degree
    numerators, denominator = bernstein_numerators(polynomial)
    assert all(value > floor * denominator for value in numerators)
    numerator, divisor = map(int, rate.split("/"))
    assert 2 * numerator * numerator > divisor * divisor
    canonical = json.dumps(polynomial, separators=(",", ":"))
    return {
        "deformed_letter": letter,
        "rho_lower_bound": rate,
        "substitution_power": substitution,
        "left_power": left_power,
        "right_power": right_power,
        "right_monomial_powers": [first_power, second_power],
        "degree_of_divided_difference": expected_degree,
        "bernstein_coefficients_checked": len(numerators),
        "every_bernstein_coefficient_strictly_greater_than": floor,
        "polynomial_coefficients_sha256": sha256(canonical.encode()).hexdigest(),
        "floating_arithmetic_used": False,
    }


def protected_profiles(deformed_letter=3):
    values = reflection.reflection_seed()
    result = []
    for heights, cut, infinity in (
            (reflection.ZERO_HEIGHTS, 0, False),
            (reflection.INFINITY_HEIGHTS, 1, True)):
        matrix = reflection.interpolation.gauged_matrix(values, heights, cut, infinity)
        # Row number is 4*x+y at zero and 4*x+u at infinity.
        protected = tuple(tuple(arithmetic.rat(value) for value in row)
                          for i, row in enumerate(matrix) if i // 4 != deformed_letter)
        assert len(protected) == 12 and all(len(row) == 16 for row in protected)
        factors = adaptive.rational_schur_factors(protected)
        orders = tuple(order for order, left, right in factors)
        base = LETTER_PROFILES[deformed_letter]
        expected = base if infinity else tuple(order + 4 for order in base)
        assert orders == expected
        result.append({
            "endpoint": "infinity" if infinity else "zero",
            "shape": [12, 16],
            "orders": orders,
            "exact_reconstruction_entries": 12 * 16,
        })
    return result


def run():
    return {
        "scope": "arbitrary auxiliary deformation supported on physical x=3; fixed natural gauges and full-word regular-channel pressure",
        "protected_letters": [0, 1, 2],
        "profiles": protected_profiles(),
        "profile_factor": "g(r) = (1+r)(1+r^3)(1+r^5+r^7)",
        "raw_two_endpoint_lower_bound": "B(s) >= g(1/s)^2",
        "common_type_letter_energies": [0, 12, 4, 16],
        "all_scale_certificate": all_scale_certificate(),
        "other_nonzero_letters": [
            {"deformed_letter": letter, "profiles": protected_profiles(letter),
             "all_scale_certificate": additional_letter_certificate(letter)}
            for letter in (1, 2)
        ],
        "auxiliary_dimension_bounded": False,
        "deformation_degree_bounded": False,
        "blocking_can_escape_this_pressure_bound": False,
        "direct_type_restricted_rank_lower_bound_claimed": False,
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence" / "protected-subtensor-pressure.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
