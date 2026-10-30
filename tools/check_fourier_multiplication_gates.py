#!/usr/bin/env python3
"""Exact small-field checks of the Fourier arithmetic-gate compiler.

This verifies local transformed tensors and a small product network.
It does not claim a faster contraction algorithm for the full network.
"""

from itertools import product
import json
from math import gcd


def polynomial_remainder(value: int, divisor: int) -> int:
    while value and value.bit_length() >= divisor.bit_length():
        value ^= divisor << (value.bit_length() - divisor.bit_length())
    return value


def polynomial_gcd(first: int, second: int) -> int:
    while second:
        first, second = second, polynomial_remainder(first, second)
    return first


class BinaryField:
    def __init__(self, polynomial: int):
        self.polynomial = polynomial
        self.degree = polynomial.bit_length() - 1
        self.size = 1 << self.degree
        self.order = self.size - 1
        value = 2
        for step in range(1, self.degree + 1):
            value = self.multiply(value, value)
            if step <= self.degree // 2:
                assert polynomial_gcd(value ^ 2, polynomial) == 1
        assert value == 2

        factors, remainder = [], self.order
        candidate = 2
        while candidate * candidate <= remainder:
            if remainder % candidate == 0:
                factors.append(candidate)
                while remainder % candidate == 0:
                    remainder //= candidate
            candidate += 1
        if remainder > 1:
            factors.append(remainder)
        self.generator = next(
            value for value in range(2, self.size)
            if all(self.power(value, self.order // prime) != 1 for prime in factors))
        self.powers = [1]
        for _ in range(1, self.order):
            self.powers.append(self.multiply(self.powers[-1], self.generator))
        assert len(set(self.powers)) == self.order
        assert self.multiply(self.powers[-1], self.generator) == 1

    def multiply(self, first: int, second: int) -> int:
        value = 0
        while second:
            if second & 1:
                value ^= first
            second >>= 1
            first <<= 1
            if first & self.size:
                first ^= self.polynomial
        return value

    def power(self, value: int, exponent: int) -> int:
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            value = self.multiply(value, value)
            exponent >>= 1
        return result

    def phase(self, exponent: int) -> int:
        return self.powers[exponent % self.order]

    def root_phase(self, root_order: int, exponent: int) -> int:
        assert self.order % root_order == 0
        return self.phase((self.order // root_order) * exponent)


def geometric(value: int, length: int, multiply) -> tuple[int, int]:
    """Return (value**length, sum_{j<length} value**j), without division."""
    if length == 0:
        return 1, 0
    power, total = geometric(value, length // 2, multiply)
    doubled_power = multiply(power, power)
    doubled_total = total ^ multiply(power, total)
    if length % 2:
        return multiply(doubled_power, value), doubled_total ^ doubled_power
    return doubled_power, doubled_total


def bilinear_closed(field: BinaryField, modulus: int, a: int, b: int, c: int) -> int:
    common = gcd(c, modulus)
    if a % common or b % common:
        return 0
    reduced = modulus // common
    solution = 0 if reduced == 1 else (
        -(b // common) * pow(c // common, -1, reduced)) % reduced
    return field.root_phase(modulus, a * solution)


def bilinear_brute(field: BinaryField, modulus: int, a: int, b: int, c: int) -> int:
    value = 0
    roots = [field.root_phase(modulus, exponent) for exponent in range(modulus)]
    for first in range(modulus):
        for second in range(modulus):
            value ^= roots[(a * first + b * second + c * first * second) % modulus]
    return value


def kloosterman_direct(field: BinaryField, modulus: int, a: int, b: int) -> int:
    value = 0
    for variable in range(modulus):
        if gcd(variable, modulus) == 1:
            value ^= field.root_phase(
                modulus, a * variable + b * pow(variable, -1, modulus))
    return value


def kloosterman_fixed_points(field: BinaryField, modulus: int, a: int, b: int) -> int:
    assert gcd(a, modulus) == gcd(b, modulus) == 1
    value = 0
    for variable in range(modulus):
        if (a * variable * variable - b) % modulus == 0:
            assert gcd(variable, modulus) == 1
            value ^= field.root_phase(modulus, 2 * a * variable)
    return value


def carry_closed(field: BinaryField, base: int, digit: int, alpha: int,
                 beta: int, gamma: int, initial: bool, final: bool) -> int:
    modulus, length = base ** 4 + 1, base ** 3
    omega = lambda exponent: field.root_phase(modulus, exponent)
    phase = field.phase
    mult = field.multiply
    a = mult(omega(alpha * base), phase(gamma))
    q = phase(beta * base + gamma)
    z = mult(phase(beta), omega(-alpha))
    if initial and final:
        return omega(alpha * digit)
    if initial:
        result = mult(omega(alpha * digit), geometric(a, length, mult)[1])
        if digit == 0:
            result ^= field.power(a, length)
        return result
    low = mult(omega(alpha * digit), geometric(z, digit + 1, mult)[1])
    if final:
        return low
    high = mult(omega(alpha * (digit + base)), phase(gamma))
    high = mult(high, field.power(z, digit + 1))
    high = mult(high, geometric(z, base - digit - 1, mult)[1])
    middle = mult(geometric(a, length, mult)[1], low ^ high)
    middle ^= mult(field.power(a, length), phase(beta * digit))
    return mult(geometric(q, 2 * length, mult)[1], middle)


def carry_points(base: int, digit: int, initial: bool, final: bool):
    modulus, carry_limit = base ** 4 + 1, 2 * base ** 4
    for incoming in range(1 if initial else carry_limit):
        for coefficient in range(modulus):
            total = incoming + coefficient
            if total % base != digit:
                continue
            outgoing = total // base
            if final and outgoing:
                continue
            yield coefficient, incoming, outgoing


def cyclic_multiply(first: int, second: int, dimension: int) -> int:
    mask, value = (1 << dimension) - 1, 0
    for shift in range(dimension):
        if (second >> shift) & 1:
            value ^= ((first << shift) | (first >> (dimension - shift))) & mask
    return value


def verify() -> dict:
    fields = {
        2: BinaryField(0b111),
        4: BinaryField(0b10011),
        6: BinaryField(0b1000011),
        8: BinaryField(0x11B),
        16: BinaryField(0x1100B),
    }
    geometric_cases = 0
    for dimension in (3, 5):
        multiply = lambda a, b: cyclic_multiply(a, b, dimension)
        for value in range(1 << dimension):
            power, total = 1, 0
            for length in range(33):
                assert geometric(value, length, multiply) == (power, total)
                total ^= power
                power = multiply(power, value)
                geometric_cases += 1

    bilinear_cases = nonunit_cases = 0
    for modulus, degree in ((3, 2), (5, 4), (9, 6), (17, 8)):
        field = fields[degree]
        for a, b, c in product(range(modulus), repeat=3):
            assert bilinear_closed(field, modulus, a, b, c) == (
                bilinear_brute(field, modulus, a, b, c))
            bilinear_cases += 1
            nonunit_cases += gcd(c, modulus) > 1
    field = fields[16]
    for index in range(32):
        modulus = 257
        a, b, c = ((index * 23) % modulus, (index * 79 + 1) % modulus,
                   (index * 41) % modulus)
        assert bilinear_closed(field, modulus, a, b, c) == (
            bilinear_brute(field, modulus, a, b, c))
        bilinear_cases += 1
        nonunit_cases += gcd(c, modulus) > 1

    interval_cases = 0
    field = fields[8]
    for alpha in range(17):
        for free_bits in range(4):
            width = 1 << free_bits
            for start in range(0, 16, width):
                z = field.root_phase(17, alpha)
                direct = 0
                for value in range(start, start + width):
                    direct ^= field.root_phase(17, alpha * value)
                compiled = field.multiply(
                    field.power(z, start), geometric(z, width, field.multiply)[1])
                assert direct == compiled
                interval_cases += 1

    carry_cases = carry_settings = 0
    for base, degree in ((2, 8), (4, 16)):
        field = fields[degree]
        modulus = base ** 4 + 1
        triples = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                   (modulus - 1, field.order - 1, field.order - 1),
                   (3, 17, 19), (7, 23, 31), (0, modulus, modulus)]
        for digit in range(base):
            for initial, final in product((False, True), repeat=2):
                settings = tuple(carry_points(base, digit, initial, final))
                for alpha, beta, gamma in triples:
                    direct = 0
                    for coefficient, incoming, outgoing in settings:
                        exponent = ((field.order // modulus) * alpha * coefficient
                                    + beta * incoming + gamma * outgoing)
                        direct ^= field.phase(exponent)
                    compiled = carry_closed(
                        field, base, digit, alpha, beta, gamma, initial, final)
                    assert direct == compiled
                    carry_cases += 1
                    carry_settings += len(settings)

    # Contract a whole one-digit product network in the transformed basis.
    product_network_cases = 0
    for base, modulus, degree in ((2, 5, 4), (4, 17, 8)):
        field = fields[degree]
        inputs = [
            geometric(field.root_phase(modulus, alpha), base, field.multiply)[1]
            for alpha in range(modulus)
        ]
        for target in range((base - 1) ** 2 + 2):
            transformed = 0
            for alpha, beta, gamma in product(range(modulus), repeat=3):
                value = field.multiply(inputs[alpha], inputs[beta])
                value = field.multiply(
                    value, bilinear_closed(field, modulus, -alpha % modulus,
                                            -beta % modulus, gamma))
                value = field.multiply(value, field.root_phase(modulus, -gamma * target))
                transformed ^= value
            expected = sum(first * second == target
                           for first in range(base) for second in range(base)) % 2
            assert transformed == expected
            product_network_cases += 1

    kloosterman_cases = 0
    for modulus, degree in ((3, 2), (5, 4), (9, 6), (15, 4),
                            (17, 8), (21, 6), (51, 8), (85, 8)):
        field = fields[degree]
        units = [value for value in range(modulus) if gcd(value, modulus) == 1]
        for a, b in product(units, repeat=2):
            assert kloosterman_direct(field, modulus, a, b) == (
                kloosterman_fixed_points(field, modulus, a, b))
            kloosterman_cases += 1

    # A sum over a unit frequency already leaves the simple M-root phase class.
    field, modulus = fields[4], 5
    kloosterman = 0
    for value in range(1, modulus):
        kloosterman ^= field.root_phase(modulus, value - pow(value, -1, modulus))
    root_values = {field.root_phase(modulus, exponent) for exponent in range(modulus)}
    assert kloosterman != 0 and kloosterman not in root_values
    assert field.multiply(kloosterman, kloosterman) ^ kloosterman ^ 1 == 0

    # Interval factors need not be invariant under c -> b/(a*c).
    weighted_direct = weighted_fixed = 0
    for variable in range(1, modulus):
        weight = 1 ^ field.root_phase(modulus, variable)
        phase = field.root_phase(modulus, variable + pow(variable, -1, modulus))
        summand = field.multiply(weight, phase)
        weighted_direct ^= summand
        if variable * variable % modulus == 1:
            weighted_fixed ^= summand
    assert weighted_direct != weighted_fixed
    # For a=0, b=1, no quadratic fixed point exists, but the sum is one.
    assert kloosterman_direct(field, modulus, 0, 1) == 1

    return {
        "scope": "exact local Fourier-gate and small contraction identities; no faster global contraction theorem",
        "verified_binary_field_degrees": list(fields),
        "zero_divisor_geometric_sum_cases": geometric_cases,
        "bilinear_gauss_sum_cases": bilinear_cases,
        "nonunit_bilinear_cases": nonunit_cases,
        "prefix_interval_cases": interval_cases,
        "carry_transform_cases": carry_cases,
        "carry_settings_summed": carry_settings,
        "whole_product_network_cases": product_network_cases,
        "unit_kloosterman_fixed_point_cases": kloosterman_cases,
        "weighted_fixed_point_negative_control": {
            "modulus": modulus,
            "interval_length": 2,
            "direct_value": weighted_direct,
            "fixed_point_value": weighted_fixed,
            "different": True,
        },
        "nonunit_kloosterman_negative_control": {
            "modulus": modulus,
            "a": 0,
            "b": 1,
            "direct_value": 1,
            "quadratic_fixed_points": 0,
        },
        "phase_closure_negative_control": {
            "modulus": modulus,
            "value": kloosterman,
            "outside_zero_or_Mth_roots": True,
            "multiplicative_order": 3,
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
