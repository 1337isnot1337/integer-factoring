#!/usr/bin/env python3
"""Exact large interval pressure certificates with bounded mantissas.

All decisions use outward integer intervals. The exponent is stored
separately, so exponentially large products do not require huge integers.
The witnesses approach the half-rate threshold from above.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


class BinaryIntervals:
    """(lower, upper, exponent) encloses [lower, upper] * 2**exponent."""

    def __init__(self, precision=256):
        if precision < 8:
            raise ValueError("at least eight mantissa bits are required")
        self.precision = precision
        self.one = self.normalize(1, 1, 0)

    @staticmethod
    def ceil_div(a, b):
        return (a + b - 1) // b

    @staticmethod
    def shift_upper(value, shift):
        if shift <= 0:
            return value << -shift
        # Avoid allocating an integer with 'shift' bits for a tiny term.
        if shift >= value.bit_length():
            return int(value != 0)
        return (value >> shift) + int(bool(value & ((1 << shift) - 1)))

    def normalize(self, lower, upper, exponent):
        assert 0 <= lower <= upper and upper > 0
        shift = upper.bit_length() - self.precision
        if shift > 0:
            lower >>= shift
            upper = self.shift_upper(upper, shift)
        elif shift < 0:
            lower <<= -shift
            upper <<= -shift
        # Rounding the upper endpoint can add one leading bit. Retaining
        # that bit is harmless; every operation still uses O(precision) bits.
        return lower, upper, exponent + shift

    def rational(self, numerator, denominator):
        if numerator <= 0 or denominator <= 0:
            raise ValueError("a positive rational is required")
        shift = self.precision + denominator.bit_length() - numerator.bit_length()
        if shift >= 0:
            scaled, divisor = numerator << shift, denominator
        else:
            scaled, divisor = numerator, denominator << -shift
        return self.normalize(scaled // divisor, self.ceil_div(scaled, divisor), -shift)

    def multiply(self, first, second):
        a, b, e = first
        c, d, f = second
        return self.normalize(a * c, b * d, e + f)

    def add(self, first, second):
        a, b, e = first
        c, d, f = second
        exponent = max(e, f)
        lower = (a >> (exponent - e)) + (c >> (exponent - f))
        upper = self.shift_upper(b, exponent - e) + self.shift_upper(d, exponent - f)
        return self.normalize(lower, upper, exponent)

    def inverse(self, value):
        lower, upper, exponent = value
        assert lower > 0
        square = 1 << (2 * self.precision)
        return self.normalize(square // upper, self.ceil_div(square, lower),
                              -exponent - 2 * self.precision)

    def power(self, value, exponent):
        if exponent < 0:
            value, exponent = self.inverse(value), -exponent
        result = self.one
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            exponent >>= 1
            if exponent:
                value = self.multiply(value, value)
        return result

    def geometric_product(self, first, ratio, count):
        result = self.one
        term = first
        for j in range(count):
            result = self.multiply(result, self.add(self.one, term))
            if j + 1 < count:
                term = self.multiply(term, ratio)
        return result

    @staticmethod
    def fraction_bounds(value):
        lower, upper, exponent = value
        scale = Fraction(2) ** exponent
        return lower * scale, upper * scale


def arithmetic_controls():
    checked = 0
    for precision in (48, 127):
        arithmetic = BinaryIntervals(precision)
        for value in (Fraction(3, 2), Fraction(7, 9), Fraction(33, 32),
                      Fraction(101, 100), Fraction(2) ** -1000, Fraction(2) ** 1000):
            encoded = arithmetic.rational(value.numerator, value.denominator)
            for exponent in (-31, -4, -1, 0, 1, 7, 29):
                lower, upper = arithmetic.fraction_bounds(arithmetic.power(encoded, exponent))
                exact = value ** exponent
                assert lower <= exact <= upper
                checked += 1
            for other in (Fraction(2) ** -10000, Fraction(2) ** 10000, Fraction(13, 17)):
                encoded_other = arithmetic.rational(other.numerator, other.denominator)
                for operation, exact in ((arithmetic.add, value + other),
                                         (arithmetic.multiply, value * other)):
                    lower, upper = arithmetic.fraction_bounds(operation(encoded, encoded_other))
                    assert lower <= exact <= upper
                    checked += 1
        root = Fraction(33, 32)
        encoded = arithmetic.rational(33, 32)
        actual = arithmetic.geometric_product(arithmetic.power(encoded, -3),
                                               arithmetic.power(encoded, -6), 13)
        exact = Fraction(1)
        for j in range(13):
            exact *= 1 + root ** (-3 - 6 * j)
        lower, upper = arithmetic.fraction_bounds(actual)
        assert lower <= exact <= upper
        checked += 1
    return checked


def pressure_ratio(n, a, b, p, q, precision=256):
    """Enclose B_n(r**a)**b / Z_(a/b)(r**a)**a exactly."""
    assert n >= 1 and a > 0 and b > 0 and p > q > 0
    arithmetic = BinaryIntervals(precision)
    root = arithmetic.rational(p, q)
    numerator = arithmetic.geometric_product(arithmetic.power(root, -a),
                                              arithmetic.power(root, -2 * a), n)
    numerator = arithmetic.power(numerator, 4 * b)
    if n % 2:
        numerator = arithmetic.multiply(numerator, arithmetic.power(root, a * b))
        denominator = arithmetic.geometric_product(arithmetic.power(root, -8 * b),
                                                    arithmetic.power(root, -8 * b), (n - 1) // 2)
        denominator = arithmetic.power(denominator, 2 * a)
        denominator = arithmetic.multiply(denominator, arithmetic.normalize(1, 1, a))
    else:
        denominator = arithmetic.geometric_product(arithmetic.power(root, -4 * b),
                                                    arithmetic.power(root, -8 * b), n // 2)
        denominator = arithmetic.power(denominator, 2 * a)
    ratio = arithmetic.multiply(numerator, arithmetic.inverse(denominator))
    return ratio


def direct_fraction_controls():
    checked = 0
    for n in (1, 2, 3, 4, 7, 8):
        a, b, p, q = 7, 3, 33, 32
        root = Fraction(p, q)
        s = root ** a
        left = s ** (n * n)
        for j in range(1, n + 1):
            left *= (1 + s ** (-(2 * j - 1))) ** 4
        right = Fraction(1)
        for mode in range(n):
            right *= 1 + root ** (-b * (8 * mode - 4 * (n - 1)))
        exact = left ** b / right ** a
        lower, upper = BinaryIntervals.fraction_bounds(pressure_ratio(n, a, b, p, q))
        assert lower <= exact <= upper
        checked += 1
    return checked


def certify(n, a, b, p, q):
    interval = pressure_ratio(n, a, b, p, q)
    lower, upper = BinaryIntervals.fraction_bounds(interval)
    assert 0 < lower <= upper < 1
    rho = Fraction(a, 4 * b)
    gap = 1 - upper
    return {
        "modes": n, "R": [a, b], "rho_upper": [rho.numerator, rho.denominator],
        "root": [p, q], "s_root_power": a, "mantissa_precision": 256,
        "pressure_ratio_interval": {"lower_mantissa": str(interval[0]),
                                    "upper_mantissa": str(interval[1]), "binary_exponent": interval[2]},
        "strict_gap_lower_bound": [str(gap.numerator), str(gap.denominator)],
        "strict_pressure_inequality": True,
        "post_half_target_reached": rho * rho < Fraction(1, 2),
        "rational_mode_probabilities": {
            "formula": "r^(-b*(8*i-4*(n-1))) / (1+r^(-b*(8*i-4*(n-1))))",
            "mode_range": [0, n - 1], "r": [p, q], "b": b},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "status": "exact finite pressure witnesses approaching one over sqrt(2) from above",
        "scope": "displayed one-species interval gauges; no improved factoring exponent",
        "independent_fraction_arithmetic_controls": arithmetic_controls(),
        "independent_unfactored_pressure_controls": direct_fraction_controls(),
        "certificates": [certify(*values) for values in (
            (63, 17, 6, 177, 176), (511, 99, 35, 5756, 5755),
            (4095, 577, 204, 204681, 204680),
            (16383, 3363, 1189, 4111597, 4111596))],
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in (
            "tools/check_large_interval_pressure.py", "tools/check_schur_pascal_complement.py",
            "tools/check_schur_one_species_interval.py")},
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
