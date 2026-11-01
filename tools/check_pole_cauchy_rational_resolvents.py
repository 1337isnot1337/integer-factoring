#!/usr/bin/env python3
"""Exact field comparisons for rational inverse contractions."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from random import Random

from interval_unit_jets import QuadraticJets
from pole_cauchy_rational_resolvents import (inverse_resolvent_bilinear_form,
                                            normalized_shifted_affine_product)
from check_pole_cauchy_circuit_moments import (node_list, normalized_inverse,
                                             qmultiply)


def valuation(value):
    if value == 0:
        return 10**9
    value = abs(value)
    return (value & -value).bit_length() - 1


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def exact_resolvent(n, left, right, xi, eta, k):
    inverse = normalized_inverse(n, left, right, k)
    rows, columns = node_list(left), node_list(right)
    result = [Fraction(0), Fraction(0)]
    for j, b in enumerate(columns):
        for i, a in enumerate(rows):
            denominator = (1 << k) * (xi - a) * (eta - b)
            for coordinate in (0, 1):
                result[coordinate] += inverse[j][i][coordinate] / denominator
    return result


def check_residue(result, exact):
    precision, scale = result["precision"], result["scale_exponent"]
    assert precision + scale == result["absolute_precision"]
    modulus = 1 << precision
    expected = []
    for coordinate in exact:
        scaled = coordinate * (1 << -scale)
        assert scaled.denominator % 2
        expected.append(scaled.numerator * pow(scaled.denominator, -1, modulus) % modulus)
    assert result["residue"] == expected, (result, expected)


def check():
    rng = Random(20261101)
    affine_comparisons = 0
    truncated_cases = zero_parts = 0
    for index in range(320):
        first = rng.randrange(-100, 101)
        step = rng.choice([x for x in range(-20, 21) if x])
        length = rng.randrange(0, 16)
        shift = rng.randrange(0, 9)
        precision = rng.randrange(1, 18)
        ring = QuadraticJets(precision, 0)
        actual_v, actual_u = normalized_shifted_affine_product(ring, first, step, length, shift)
        exact = (1, 0)
        for i in range(length):
            real = first + step * i
            truncated_cases += bool(real and valuation(real) >= shift)
            zero_parts += real == 0
            exact = qmultiply(exact, (real, 1 << shift))
        expected_v = min(valuation(exact[0]), valuation(exact[1]))
        assert actual_v == expected_v
        expected_u = tuple((coordinate >> expected_v) % ring.modulus for coordinate in exact)
        assert actual_u[0] == expected_u
        affine_comparisons += 1

    exact_comparisons = 0
    parity_cases = {"even_even": 0, "even_odd": 0, "odd_even": 0, "odd_odd": 0}
    negative_ratio_valuations = equal_query_products = nonseparated_cases = 0
    min_scale, max_precision = 0, 0
    for index in range(72):
        size = rng.randrange(1, 7)
        if index % 3:
            left = ((2 * rng.randrange(1, 7) + 1, size),)
            right = ((2 * rng.randrange(1, 7) + 1, size),)
        else:
            left = tuple((a, 1) for a in sorted(rng.sample(range(1, 40, 2), size)))
            right = tuple((b, 1) for b in sorted(rng.sample(range(1, 40, 2), size)))
        rows, columns = node_list(left), node_list(right)
        xi = rng.choice([x for x in range(1, 50) if x not in rows and x % 2 == index % 2])
        eta = rng.choice([x for x in range(1, 50) if x not in columns and x % 2 == (index // 2) % 2])
        n = xi * eta if index % 8 == 3 else 2 * rng.randrange(1, 2048) + 1
        assert n % 2
        horizon = max([1] + ([abs(xi - a).bit_length() for a in rows] if xi % 2 else [])
                      + ([abs(eta - b).bit_length() for b in columns] if eta % 2 else []))
        k = horizon + index % 3
        absolute_precision = rng.randrange(1, 10)
        result = inverse_resolvent_bilinear_form(n, left, right, xi, eta, k, absolute_precision)
        exact = exact_resolvent(n, left, right, xi, eta, k)
        check_residue(result, exact)
        exact_comparisons += 1
        parity_cases[("odd" if xi % 2 else "even") + "_" + ("odd" if eta % 2 else "even")] += 1
        negative_ratio_valuations += result["ratio_valuation"] < 0
        equal_query_products += xi * eta == n
        nonseparated_cases += max(abs(n - a * b) for a in rows for b in columns) >= 1 << k
        min_scale = min(min_scale, result["scale_exponent"])
        max_precision = max(max_precision, result["precision"])

    # Every unit appearing in an even/even query is congruent to 1 mod 2.
    # Its order modulo 2^P divides 2^(P-1). A complete odd residue cycle
    # has length 2^(P-1). Thus 2^100 leading nodes contribute identity
    # when 100 >= 2P-2; the remaining three nodes admit a literal inverse.
    n, large_size, xi, eta, k, absolute_precision = (1 << 201) + 123, (1 << 100) + 3, 2, 4, 1, 16
    left, right = ((1, large_size),), ((5, large_size),)
    large = inverse_resolvent_bilinear_form(n, left, right, xi, eta, k, absolute_precision)
    assert 100 >= 2 * large["precision"] - 2
    check_residue(large, exact_resolvent(n, ((1, 3),), ((5, 3),), xi, eta, k))

    # A singular rational weight with a binary-encoded 61-bit dimension.
    # This is a non-enumerating execution check; unlike the preceding
    # control, it is not an independent exact large-dimension reference.
    odd_size = (1 << 60) + 3
    odd_xi, odd_eta = 2 * odd_size + 5, 2 * odd_size + 7
    odd_n = odd_xi * odd_eta - 2
    large_odd = inverse_resolvent_bilinear_form(odd_n, ((1, odd_size),), ((5, odd_size),),
                                              odd_xi, odd_eta, 64, 8)
    assert large_odd["absolute_precision"] == 8

    # The distance horizon is needed for a polynomial absolute-precision
    # promise. With m=2^s, both external poles 2m+1, and k=1, the ratio
    # has valuation 2-2m and the final contraction has valuation -2m.
    missing_horizon_controls = []
    for size in (2, 4, 8):
        pole = 2 * size + 1
        exact = exact_resolvent(1, ((1, size),), ((1, size),), pole, pole, 1)
        value_v = min(valuation(x.numerator) - valuation(x.denominator) for x in exact)
        assert value_v == -2 * size
        missing_horizon_controls.append({"dimension": size, "query": pole,
                                         "exact_scalar_valuation": value_v})

    rejected = 0
    invalid = [
        (15, ((1, 3),), ((1, 3),), 3, 8, 8, 4),
        (15, ((1, 3),), ((1, 3),), 8, 3, 8, 4),
        (15, ((1, 3),), ((1, 2),), 8, 8, 8, 4),
        (15, ((1, 20),), ((1, 20),), 41, 43, 2, 4),
        (14, ((1, 3),), ((1, 3),), 8, 8, 8, 4),
    ]
    for arguments in invalid:
        try:
            inverse_resolvent_bilinear_form(*arguments)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Invalid input was accepted")
    return {"status": "passed", "seed": 20261101,
            "affine_exact_comparisons": affine_comparisons,
            "affine_rational_parts_at_or_above_shift": truncated_cases,
            "affine_zero_rational_parts": zero_parts,
            "inverse_exact_comparisons": exact_comparisons,
            "query_parity_cases": parity_cases,
            "negative_ratio_valuation_cases": negative_ratio_valuations,
            "equal_query_product_cases": equal_query_products,
            "original_rectangle_not_separated_cases": nonseparated_cases,
            "smallest_output_scale": min_scale,
            "largest_small_test_precision": max_precision,
            "large_independent_control": {"dimension": large_size, "dimension_bits": large_size.bit_length(),
                                          "number_bits": n.bit_length(), "queries": [xi, eta],
                                          "displacement": k, "result": large,
                                          "reference": "complete unit cycles and a literal 3x3 inverse"},
            "large_odd_query_execution": {"dimension": odd_size,
                                           "dimension_bits": odd_size.bit_length(),
                                           "number_bits": odd_n.bit_length(),
                                           "queries": [odd_xi, odd_eta], "displacement": 64,
                                           "result": large_odd,
                                           "reference": "no independent large-dimension inverse"},
            "missing_horizon_counterexamples": missing_horizon_controls,
            "rejected_invalid_inputs": rejected,
            "large_controls_sha256": digest([large, large_odd])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    evidence = check()
    rendered = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
