#!/usr/bin/env python3
"""Exact controls for rectangular isometric-prefix Bezout denominators."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb, factorial, prod
from pathlib import Path
import random

from check_balanced_bezout_denominator import (
    divided_differences, evaluate, isometry, newton_polynomial, order,
)
from check_interval_bezout_precision import multiply, remainder
from rectangular_bezout_scale import rectangular_bezout_scale


def binomial(value, degree):
    return prod((Fraction(value) - i for i in range(degree)), start=Fraction(1)) / factorial(degree)


def random_isometry(length, rng):
    bits = max(1, (length - 1).bit_length()) + 2
    swaps = {}
    values = []
    for index in range(length):
        value = 0
        for depth in range(bits):
            key = (depth, index % (1 << depth))
            if key not in swaps:
                swaps[key] = rng.randrange(2)
            value |= (((index >> depth) & 1) ^ swaps[key]) << depth
        values.append(value + (rng.randrange(4) << bits))
    assert isometry(values)
    return values


def c_functional(nodes, surplus):
    a = len(nodes) - 1
    values = [binomial(node, a + surplus) for node in nodes]
    result = factorial(a) * divided_differences(nodes, values)[-1]
    assert order(result) is None or order(result) >= 0
    return result


def elementary(nodes, degree):
    result = [Fraction(1)] + [Fraction(0)] * degree
    for node in nodes:
        for i in range(degree, 0, -1):
            result[i] += node * result[i - 1]
    return result


def complete(nodes, degree):
    result = [Fraction(1)] + [Fraction(0)] * degree
    for node in nodes:
        for i in range(1, degree + 1):
            result[i] += node * result[i - 1]
    return result


def c_from_moments(nodes, surplus):
    a = len(nodes) - 1
    ee = elementary(range(a + surplus), surplus)
    hh = complete(nodes, surplus)
    numerator = sum(((-1) ** (surplus - j) * ee[surplus - j] * hh[j]
                     for j in range(surplus + 1)), Fraction(0))
    return numerator / prod(range(a + 1, a + surplus + 1))


def mahler_coefficient(r, t):
    return Fraction((-1) ** t * (1 << (r + t)) * factorial(r + t),
                    prod(1 + 2 * v for v in range(-r, t + 1)))


def literal_difference(r, t):
    return sum((Fraction((-1) ** (r + t - i - j) * comb(r, i) * comb(t, j),
                         1 + 2 * (j - i))
                for i in range(r + 1) for j in range(t + 1)), Fraction(0))


def normalized_rectangle_unit(xx, yy, precision):
    a, b = len(xx) - 1, len(yy) - 1
    cx = [c_functional(xx, i) for i in range(precision)]
    cy = [c_functional(yy, j) for j in range(precision)]
    total = Fraction(0)
    for i in range(precision):
        for j in range(precision - i):
            denominator = (prod(1 - 2 * a - 2 * ell for ell in range(1, i + 1))
                           * prod(1 + 2 * b + 2 * ell for ell in range(1, j + 1)))
            numerator = ((-1) ** j * (1 << (i + j))
                         * prod(range(a + b + 1, a + b + i + j + 1)))
            total += Fraction(numerator, denominator) * cx[i] * cy[j]
    return total


def polynomial_product(factors):
    result = [1]
    for factor in factors:
        result = multiply(result, factor)
    return result


def quotient_monic(numerator, denominator):
    work = list(map(Fraction, numerator))
    degree = len(denominator) - 1
    result = [Fraction(0)] * max(1, len(work) - degree)
    for i in range(len(work) - 1, degree - 1, -1):
        value = work[i]
        result[i - degree] = value
        for j, coefficient in enumerate(denominator):
            work[i - degree + j] -= value * coefficient
    assert all(value == 0 for value in work)
    return result


def hensel_scale(n, aa, bb, height):
    m, s = len(aa), len(bb)
    unit = 1 << height
    pp = polynomial_product([[-(a - 1) // unit, 1] for a in aa])
    qq = polynomial_product([[(n - b) // unit, -b] for b in bb])
    precision = height * m + 2
    inverse = [comb(s + i - 1, i) % 2 for i in range(m)]
    bits = 1
    while bits < precision:
        bits = min(2 * bits, precision)
        modulus = 1 << bits
        error = [-value for value in remainder(multiply(qq, inverse, modulus), pp, modulus)]
        error += [0] * (m - len(error))
        error[0] += 2
        inverse = remainder(multiply(inverse, error, modulus), pp, modulus)
    inverse += [0] * (m - len(inverse))
    check = remainder(multiply(qq, inverse, 1 << precision), pp, 1 << precision)
    assert check[0] == 1 and not any(check[1:])
    values = [(value & -value).bit_length() - 1 if value else precision for value in inverse]
    return max(height * (s + i) - value for i, value in enumerate(values))


def shifted_coefficients(polynomial):
    result = [Fraction(0)] * len(polynomial)
    for degree, coefficient in enumerate(polynomial):
        for i in range(degree + 1):
            result[i] += coefficient * comb(degree, i)
    return result


def check_bezout(n, aa, bb, height):
    m, s = len(aa), len(bb)
    pp = polynomial_product([[-a, 1] for a in aa])
    qq = polynomial_product([[n, -b] for b in bb])
    values = [Fraction(1, evaluate(qq, a)) for a in aa]
    coefficients = divided_differences(aa, values)
    expected = [-height * (s + k) + order(comb(s + k - 1, k)) for k in range(m)]
    assert [order(value) for value in coefficients] == expected
    scale = -min(expected)
    vv = newton_polynomial(aa, coefficients)
    assert min(order(value) for value in vv if value) == -scale
    numerator = [-value for value in multiply(vv, qq)]
    numerator[0] += 1
    uu = quotient_monic(numerator, pp)
    assert len(uu) <= s
    assert all((value * (1 << scale)).denominator % 2 for value in uu + vv)
    assert hensel_scale(n, aa, bb, height) == scale
    tail_checks = 0
    if height == 1:
        assert rectangular_bezout_scale(m, s)["exponent"] == scale
        defect = m + s - 1 - scale
        for polynomial, top in ((uu, s - 1), (vv, m - 1)):
            shifted = shifted_coefficients(polynomial)
            shifted += [Fraction(0)] * (top + 1 - len(shifted))
            for degree, coefficient in enumerate(shifted):
                depth = top - degree
                if coefficient:
                    assert scale + order(coefficient) >= depth - defect
                tail_checks += 1
    return {"n": n, "row_count": m, "column_count": s, "height": height,
            "dyadic_denominator": scale, "newton_valuations": expected}, tail_checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(2026110102)
    counts = {"double_forward_difference_identities": 0,
              "finite_double_mahler_reconstructions": 0,
              "parity_divided_difference_valuations": 0,
              "integral_normalized_functionals": 0,
              "low_moment_identities": 0,
              "normalized_unit_truncations": 0,
              "exact_inverse_newton_coefficients": 0,
              "independent_hensel_certificates": 0,
              "shifted_polynomial_tail_bounds": 0}
    for r in range(9):
        for t in range(9):
            actual = literal_difference(r, t)
            assert actual == mahler_coefficient(r, t)
            assert order(actual) == r + t + (r + t) - (r + t).bit_count()
            counts["double_forward_difference_identities"] += 1
    for x in range(8):
        for y in range(8):
            actual = sum((mahler_coefficient(r, t) * comb(x, r) * comb(y, t)
                          for r in range(x + 1) for t in range(y + 1)), Fraction(0))
            assert actual == Fraction(1, 1 + 2 * (y - x))
            counts["finite_double_mahler_reconstructions"] += 1
    parity_examples = []
    for m in range(1, 9):
        for s in range(1, 9):
            for variant in range(2):
                xx = list(map(Fraction, random_isometry(m, rng)))
                yy = list(map(Fraction, random_isometry(s, rng)))
                if variant:
                    xx = [(3 * x + 2) / (2 * x + 1) for x in xx]
                    yy = [(5 * y + 1) / (2 * y + 1) for y in yy]
                assert isometry(xx) and isometry(yy)
                nodes = [2 * x for x in xx] + [1 + 2 * y for y in yy]
                value = divided_differences(nodes, [0] * m + [1] * s)[-1]
                a, b = m - 1, s - 1
                leading = Fraction((-1) ** b * comb(a + b, a),
                                   prod(1 + 2 * v for v in range(-a, b + 1)))
                assert order(value) == order(comb(a + b, a))
                assert order(value / leading) == 0
                counts["parity_divided_difference_valuations"] += 1
                for group in (xx, yy):
                    for surplus in range(6):
                        actual = c_functional(group, surplus)
                        assert actual == c_from_moments(group, surplus)
                        if surplus == 0:
                            assert actual == 1
                        counts["integral_normalized_functionals"] += 1
                        counts["low_moment_identities"] += 1
                for precision in (1, 3, 6):
                    approximation = normalized_rectangle_unit(xx, yy, precision)
                    error = value / leading - approximation
                    assert not error or order(error) >= precision
                    counts["normalized_unit_truncations"] += 1
                if m in (1, 4, 8) and s in (1, 3, 8) and not variant:
                    parity_examples.append({"even_count": m, "odd_count": s,
                                            "valuation": order(value)})
    certificates = []
    for m in range(1, 9):
        for s in range(1, 9):
            for height in (1, 2, 4):
                xx, yy = random_isometry(m, rng), random_isometry(s, rng)
                aa = [1 + (1 << (height + 1)) * x for x in xx]
                bb = [1 + (1 << (height + 1)) * y for y in yy]
                n = 1 + (1 << height) + (1 << (height + 1)) * rng.randrange(50)
                certificate, tail_checks = check_bezout(n, aa, bb, height)
                counts["exact_inverse_newton_coefficients"] += m
                counts["independent_hensel_certificates"] += 1
                counts["shifted_polynomial_tail_bounds"] += tail_checks
                if m in (1, 4, 8) and s in (1, 3, 8):
                    certificates.append(certificate)
    large_scales = []
    for m, s in (((1 << 200) + 13, (1 << 137) + 7),
                 (1 << 200, 1 << 200), (1, 1 << 400), ((1 << 200) + 1, 3)):
        scale = rectangular_bezout_scale(m, s)
        assert scale["exponent"] == rectangular_bezout_scale(s, m)["exponent"]
        large_scales.append({"rows": m, "columns": s, **scale,
                             "scope": "Binary scale computation only; no nodes enumerated"})
    result = {"theorem": "Exact rectangular isometric-prefix Bezout denominators",
              "arithmetic": "Exact rational divided differences and independent modular Hensel inverses",
              **counts, "parity_examples": parity_examples,
              "certificate_examples": certificates, "large_binary_scales": large_scales,
              "scope": "Separated residue classes, with no factor incidence. No general compressed GCD or improved factoring exponent."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
