#!/usr/bin/env python3
"""Finite exact providers and independent controls for the supplied-jet decoder."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb, factorial, prod
from pathlib import Path
import random

from check_pole_smith_stabilization import Quad, W
from check_pole_cauchy_conditioning import mod_pair
from check_interval_bezout_precision import multiply
from check_balanced_bezout_denominator import order
from displaced_bezout_extraction import decode_factor_moments, extraction_precision
from regularized_bezout_factor_moments import (
    factor_moments_from_bezout_jet, regularized_bezout_parameters,
)


def product_polynomial(nodes):
    result = [1]
    for node in nodes:
        result = multiply(result, [-node, 1])
    return result


def unit_residue(value, precision):
    value = Fraction(value)
    modulus = 1 << precision
    assert value.denominator % 2
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def literal_inverse_polynomial(n, aa, bb, padding):
    t = (1 << padding) * W
    values = []
    for a in aa:
        value = Quad.of(1)
        for b in bb:
            value = value * (n + t - a * b)
        values.append(1 / value)
    coefficients = [values[0]]
    for degree in range(1, len(aa)):
        values = [(values[i + 1] - values[i]) / (aa[i + degree] - aa[i])
                  for i in range(len(values) - 1)]
        coefficients.append(values[0])
    polynomial = [Quad()] * len(aa)
    basis = [1]
    for degree, coefficient in enumerate(coefficients):
        for j, value in enumerate(basis):
            polynomial[j] = polynomial[j] + coefficient * value
        if degree + 1 < len(aa):
            basis = multiply(basis, [-aa[degree], 1])
    return polynomial


def fixed_residue_polynomial(n, aa, bb):
    ra, rb = len(aa), len(bb)
    ea, eb = 2 * ra - 2 - (ra.bit_length() - 1), 2 * rb - 2 - (rb.bit_length() - 1)
    positions = {b: j for j, b in enumerate(bb)}
    result = [Fraction(0)] * ra
    weights = []
    for i, a in enumerate(aa):
        if n % a or n // a not in positions:
            continue
        j = positions[n // a]
        da = (1 << (ra - 1)) * (-1) ** (ra - 1 - i) * factorial(i) * factorial(ra - 1 - i)
        db = (1 << (rb - 1)) * (-1) ** (rb - 1 - j) * factorial(j) * factorial(rb - 1 - j)
        alpha, beta = da // (1 << ea), db // (1 << eb)
        assert alpha % 2 and beta % 2
        weight = Fraction(1, a ** (rb - 1) * alpha * beta)
        basis = product_polynomial(aa[:i] + aa[i + 1:])
        for k, value in enumerate(basis):
            result[k] += weight * value
        weights.append((a, weight))
    return result, weights


def shifted(polynomial):
    result = [Quad()] * len(polynomial)
    for degree, value in enumerate(polynomial):
        for j in range(degree + 1):
            result[j] = result[j] + comb(degree, j) * value
    return result


def literal_case(n, first_a, rows, first_b, columns, bound):
    aa = tuple(first_a + 2 * i for i in range(rows))
    bb = tuple(first_b + 2 * j for j in range(columns))
    ell = max(1, (aa[-1] - aa[0]).bit_length())
    budget = extraction_precision(bound, aa[-1], ell)
    precision = budget["precision"]
    contract = regularized_bezout_parameters(n, first_a, rows, first_b, columns, precision)
    padding, height = contract["padding_exponent"], contract["height_bound"]
    total = contract["row_derivative_order"] + contract["column_derivative_order"]
    polynomial = literal_inverse_polynomial(n, aa, bb, padding)
    primitive_scale = -min(value.valuation() for value in polynomial if value)
    residue_polynomial, weights = fixed_residue_polynomial(n, aa, bb)
    degree = len(weights)
    assert degree <= bound
    if weights:
        mu = min(order(value) for value in residue_polynomial if value)
        roots = [a for a, _ in weights]
        derivative_orders = [sum(order(a - b) for b in roots if b != a) for a in roots]
        assert 0 <= mu <= min(derivative_orders) <= (degree - 1) * ell
        assert primitive_scale == padding + total - mu
        boundary = height + mu
        low = literal_inverse_polynomial(n, aa, bb, boundary)
        high = literal_inverse_polynomial(n, aa, bb, boundary + 1)
        assert -min(value.valuation() for value in low if value) == boundary + total - mu
        assert -min(value.valuation() for value in high if value) == boundary + 1 + total - mu
        assert primitive_scale > total + height - 1
    else:
        roots, mu, derivative_orders = [], None, []
        assert primitive_scale <= total + height - 1
    cleared = [W * (1 << (padding + total)) * value for value in polynomial]
    for actual, expected in zip(cleared, residue_polynomial):
        error = actual - expected
        assert not error or error.valuation() >= precision
    shifted_polynomial = shifted(cleared)
    reverse = []
    for depth, value in enumerate(reversed(shifted_polynomial)):
        actual = mod_pair(value, precision)
        assert actual[1] == 0
        assert actual[0] % (1 << min(depth, precision)) == 0
        reverse.append(actual[0])
    supplied = reverse[:min(rows, precision, 2 * bound)]
    converted = factor_moments_from_bezout_jet(first_a, rows, precision, supplied, 2 * bound)
    expected_moments = [unit_residue(sum((weight * a ** j for a, weight in weights), Fraction(0)), precision)
                        for j in range(2 * bound)]
    assert converted["ordinary_moments"] == expected_moments
    decoded = decode_factor_moments(converted["ordinary_moments"], factor_bound=bound,
                                    max_factor=aa[-1], separation_bits=ell, precision=precision)
    assert decoded["degree"] == degree
    assert decoded["factor_polynomial"] == product_polynomial(roots)
    return {"n": n, "first_a": first_a, "rows": rows, "first_b": first_b,
            "columns": columns, "factor_bound": bound, "factor_rows": roots,
            "primitive_scale": primitive_scale, "residue_content": mu,
            "factor_derivative_orders": derivative_orders,
            "regularization_contract": contract,
            "decoded_factor_polynomial": decoded["factor_polynomial"],
            "decoded_hankel_order": decoded["hankel_determinant_valuation"],
            "moment_count": 2 * bound,
            "supplied_jet_coefficients": len(supplied),
            "scope": "Finite provider explicitly interpolates the complete small interval."}


def generic_control(nodes, weights, bound, max_factor):
    assert len(nodes) == len(weights) and len(nodes) <= bound
    ell = max(1, max_factor.bit_length())
    budget = extraction_precision(bound, max_factor, ell)
    precision = budget["precision"]
    moments = [unit_residue(sum((weight * a ** j for a, weight in zip(nodes, weights)), Fraction(0)), precision)
               for j in range(2 * bound)]
    decoded = decode_factor_moments(moments, factor_bound=bound, max_factor=max_factor,
                                    separation_bits=ell, precision=precision)
    assert decoded["factor_polynomial"] == product_polynomial(nodes)
    return decoded


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [
        (3, 1, 1, 1, 1, 1), (3, 3, 1, 1, 2, 1),
        (3, 1, 2, 1, 2, 2), (5, 1, 4, 1, 2, 2),
        (15, 1, 8, 1, 8, 4), (45, 1, 8, 1, 8, 4),
        (81, 1, 8, 1, 8, 2), (17, 1, 8, 1, 8, 3),
        (35, 3, 4, 1, 8, 3), (77, 3, 8, 3, 4, 2),
        (101, 1, 4, 1, 8, 2), (105, 3, 8, 1, 8, 5),
    ]
    exact_cases = [literal_case(*case) for case in cases]
    rng = random.Random(2026110103)
    controls = []
    for degree in range(7):
        for trial in range(5):
            nodes = sorted(rng.sample(range(1, 90), degree))
            weights = [Fraction(rng.choice((-7, -5, -3, -1, 1, 3, 5, 7)),
                                rng.choice((1, 3, 5, 7, 9))) for _ in nodes]
            controls.append(generic_control(nodes, weights, max(1, degree + trial % 2), 90))
    vanishing_prefix = generic_control([3, 7], [Fraction(1), Fraction(-1)], 3, 12)
    assert vanishing_prefix["nonzero_leading_hankel_sizes"] == [2]
    large_nodes = [(1 << 99) + 3, (1 << 100) + 5, (1 << 101) + 7]
    large = generic_control(large_nodes, [Fraction(1, 3), Fraction(-5, 7), Fraction(9, 11)], 4, 1 << 102)
    result = {"theorem": "Primitive displaced Bezout sensitivity and supplied-moment factor recovery",
              "finite_exact_interval_providers": len(exact_cases),
              "generic_unit_moment_decodings": len(controls),
              "primitive_scale_boundary_pairs": sum(bool(case["factor_rows"]) for case in exact_cases),
              "exact_interval_cases": exact_cases,
              "vanishing_leading_hankel_control": vanishing_prefix,
              "large_integer_node_control": {"nodes": large_nodes, **large,
                                             "scope": "Supplied unit-weight moments, not an implicit interval jet provider"},
              "scope": "The decoder is polynomial in the supplied factor bound and bit parameters. Computing a global regularized Bezout jet remains unimplemented."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
