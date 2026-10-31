#!/usr/bin/env python3
"""Exact controls for the all-digit interval pressure floor.

The all-size analytic theorem is in the accompanying research note.
These finite controls use rational arithmetic and integer polynomials.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path


def digit_sum(value: Fraction, digits: int) -> Fraction:
    result = Fraction(0)
    term = Fraction(1)
    for _ in range(digits):
        result += term
        term *= value
    return result


def original_ratio(n: int, q: int, a: int, b: int,
                   root: Fraction) -> Fraction:
    """B(s)^b/Z(s)^a, where s=root^a and R=a/b."""
    bulk = root ** (a * (q - 1) * n * n)
    for j in range(1, n + 1):
        bulk *= digit_sum(root ** (-a * (2 * j - 1)), q) ** 4
    physical = Fraction(1)
    for j in range(n):
        physical *= digit_sum(root ** (-b * (8 * j - 4 * (n - 1))), q)
    return bulk ** b / physical ** a


def cancelled_ratio(n: int, q: int, a: int, b: int,
                    root: Fraction) -> Fraction:
    numerator = root ** (a * b * (q - 1)) if n % 2 else Fraction(1)
    for j in range(1, n + 1):
        numerator *= digit_sum(root ** (-a * (2 * j - 1)), q) ** (4 * b)
    if n % 2:
        denominator = Fraction(q) ** a
        for j in range(1, (n - 1) // 2 + 1):
            denominator *= digit_sum(root ** (-8 * b * j), q) ** (2 * a)
    else:
        denominator = Fraction(1)
        for j in range(1, n // 2 + 1):
            denominator *= digit_sum(root ** (-4 * b * (2 * j - 1)), q) ** (2 * a)
    return numerator / denominator


def polynomial_product(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def polynomial_power(value, exponent):
    result = [1]
    for _ in range(exponent):
        result = polynomial_product(result, value)
    return result


def shift_polynomial(value, shift):
    return [sum(value[j] * comb(j, i) * shift ** (j - i)
                for j in range(i, len(value))) for i in range(len(value))]


def derivative_coefficients(maximum):
    """Build 2F from independent sinh/cosh Taylor coefficients."""
    sinh = [Fraction(int(i % 2), factorial(i)) for i in range(maximum + 1)]
    cosh = [Fraction(int(not i % 2), factorial(i)) for i in range(maximum + 1)]
    square = polynomial_product(sinh, sinh)
    cross = polynomial_product(sinh, cosh)
    result = []
    for degree in range(maximum + 1):
        value = -square[degree]
        if degree >= 2:
            value += square[degree - 2]
        if degree >= 1:
            value -= cross[degree - 1]
        if degree == 2:
            value += 2
        result.append(2 * value)
    for degree, value in enumerate(result):
        if degree < 4 or degree % 2:
            assert value == 0
        else:
            r = degree // 2
            expected = Fraction(2 ** (2 * r - 2) * (4 * r * r - 6 * r - 4),
                                factorial(2 * r))
            assert value == expected
            assert value >= 0
            assert value > 0 or r == 2
    return result


def small_interval_certificates():
    # Exact strict comparison used in the decreasing exponential kernel.
    assert 79 ** 41 * 100 ** 200 > 7 ** 41 * 159 ** 200
    assert Fraction(141, 100) ** 2 < 2 < Fraction(142, 100) ** 2
    entropy = []
    for q in range(2, 6):
        weights = [comb(q - j + 1, 2) for j in range(q)]
        total = sum(weights)
        assert total == comb(q + 2, 3)
        assert 4 * sum(j * weights[j] for j in range(q)) == (q - 1) * total
        left = total ** (4 * total)
        right = q ** (3 * total)
        for weight in weights:
            right *= weight ** (4 * weight)
        assert left > right
        entropy.append({'digits': q, 'weights': weights,
                        'strict_entropy_comparison': True})
    left_poly = [625 * x for x in polynomial_product(
        [0, 1], polynomial_product(polynomial_power([1, 1], 4),
                                   polynomial_power([2, 1], 4)))]
    right_poly = [81 * x for x in polynomial_power([1, 6, 3], 4)]
    difference = [value - (right_poly[i] if i < len(right_poly) else 0)
                  for i, value in enumerate(left_poly)]
    shifted = shift_polynomial(difference, 6)
    assert all(value > 0 for value in shifted)
    for q in range(2, 129):
        weights = [comb(q - j + 1, 2) for j in range(q)]
        assert sum(weights) == comb(q + 2, 3)
        assert 4 * sum(j * weights[j] for j in range(q)) == (q - 1) * sum(weights)
        assert 60 * sum(w * w for w in weights) == q * (q + 1) * (q + 2) * (3 * q * q + 6 * q + 1)
    return {'two_mode_rational_kernel_certificate': True,
            'one_mode_small_digits': entropy,
            'renyi_polynomial_at_q_equal_6_plus_x': shifted}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    products = 0
    for n in range(1, 7):
        for q in range(2, 8):
            for a, b in ((3, 1), (5, 2), (17, 6)):
                for root in (Fraction(3, 2), Fraction(5, 4)):
                    assert original_ratio(n, q, a, b, root) == cancelled_ratio(n, q, a, b, root)
                    products += 1
    periodic = 0
    for q in range(2, 129):
        partial = 0
        for r in range(1, 8 * q + 1):
            partial += 1 - q * int(r % q == 0)
            assert partial == r % q
            periodic += 1
    derivative = derivative_coefficients(48)
    certificates = small_interval_certificates()
    result = {
        'status': 'Exact finite controls; all-size floor and fixed-digit asymptotic are proved in the note.',
        'source_predecessor': '18cec06deac47c8f142ea67a60e2278ac2234d67',
        'rational_product_equalities': products,
        'periodic_partial_sum_checks': periodic,
        'derivative_coefficients_checked_through_degree': len(derivative) - 1,
        **certificates,
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('checker_sha256', 'renyi_polynomial_at_q_equal_6_plus_x')}))


if __name__ == '__main__':
    main()
