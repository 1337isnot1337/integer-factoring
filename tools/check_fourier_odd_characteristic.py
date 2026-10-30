#!/usr/bin/env python3
"""Exact finite checks for the odd-characteristic Fourier gate compiler.

The checker caches the small field's primitive-root cycle to accelerate
direct sums. The polynomial-space production argument uses exponentiation
instead and does not require that table.
"""

from itertools import product
from math import gcd
import json


def coefficients(value, prime):
    result = []
    while value:
        result.append(value % prime)
        value //= prime
    return result


def encode(values, prime):
    result = 0
    for value in reversed(values):
        result = result * prime + value % prime
    return result


def polynomial_remainder(value, divisor, prime):
    current, base = coefficients(value, prime), coefficients(divisor, prime)
    inverse = pow(base[-1], -1, prime)
    while len(current) >= len(base):
        shift, scalar = len(current) - len(base), current[-1] * inverse % prime
        for index, value in enumerate(base):
            current[shift + index] = (current[shift + index] - scalar * value) % prime
        while current and current[-1] == 0:
            current.pop()
    return encode(current, prime)


class PrimePowerField:
    def __init__(self, prime, degree):
        self.prime, self.degree = prime, degree
        self.size = prime ** degree
        self.order = self.size - 1
        self.polynomial = next(
            candidate for candidate in range(self.size + 1, 2 * self.size)
            if candidate % prime
            and all(polynomial_remainder(candidate, divisor, prime)
                    for divisor_degree in range(1, degree // 2 + 1)
                    for divisor in range(prime ** divisor_degree,
                                         2 * prime ** divisor_degree))
        )
        self.modulus = coefficients(self.polynomial, prime)[:-1]
        remainder, candidate, factors = self.order, 2, []
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
            if all(self.power(value, self.order // divisor) != 1 for divisor in factors)
        )
        self.powers = [1]
        for _ in range(1, self.order):
            self.powers.append(self.multiply(self.powers[-1], self.generator))
        assert len(set(self.powers)) == self.order
        assert self.multiply(self.powers[-1], self.generator) == 1

    def vector(self, value):
        result = coefficients(value, self.prime)
        return result + [0] * (self.degree - len(result))

    def add(self, first, second):
        result, place = 0, 1
        for _ in range(self.degree):
            result += ((first + second) % self.prime) * place
            first //= self.prime
            second //= self.prime
            place *= self.prime
        return result

    def scale(self, value, scalar):
        return encode([scalar * entry % self.prime for entry in self.vector(value)], self.prime)

    def multiply(self, first, second):
        first, second = self.vector(first), self.vector(second)
        values = [0] * (2 * self.degree - 1)
        for i, a in enumerate(first):
            if a:
                for j, b in enumerate(second):
                    values[i + j] += a * b
        for position in range(len(values) - 1, self.degree - 1, -1):
            scalar = values[position] % self.prime
            if scalar:
                for index, coefficient in enumerate(self.modulus):
                    values[position - self.degree + index] -= scalar * coefficient
        return encode(values[:self.degree], self.prime)

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            value = self.multiply(value, value)
            exponent >>= 1
        return result

    def phase(self, exponent):
        return self.powers[exponent % self.order]

    def root_phase(self, order, exponent):
        assert self.order % order == 0
        return self.phase((self.order // order) * exponent)


def geometric(field, value, length):
    if length == 0:
        return 1, 0
    power, total = geometric(field, value, length // 2)
    doubled_power = field.multiply(power, power)
    doubled_sum = field.add(total, field.multiply(power, total))
    if length % 2:
        return field.multiply(doubled_power, value), field.add(doubled_sum, doubled_power)
    return doubled_power, doubled_sum


def carry_closed(field, base, modulus, carry_limit, digit,
                 alpha, beta, gamma, initial, final):
    omega = lambda exponent: field.root_phase(modulus, exponent)
    phase, mult = field.phase, field.multiply
    high, remainder = divmod(modulus, base)
    carry_high, carry_remainder = divmod(carry_limit, base)
    u = mult(omega(alpha * base), phase(gamma))
    v = phase(beta * base + gamma)
    z = mult(omega(-alpha), phase(beta))
    if initial and final:
        return omega(alpha * digit)
    if initial:
        return mult(omega(alpha * digit),
                    geometric(field, u, high + (digit < remainder))[1])
    if final:
        return mult(omega(alpha * digit), geometric(field, z, digit + 1)[1])
    cuts = sorted({max(0, min(base, value)) for value in (
        0, base, digit + 1, carry_remainder,
        digit - remainder + 1, digit + base - remainder + 1)})
    assert len(cuts) <= 6
    result = 0
    for lower, upper in zip(cuts, cuts[1:]):
        wrap = lower > digit
        residue = digit - lower + base * wrap
        coefficient_count = high + (residue < remainder)
        carry_count = carry_high + (lower < carry_remainder)
        value = mult(omega(alpha * (digit + base * wrap)), phase(gamma * wrap))
        value = mult(value, field.power(z, lower))
        value = mult(value, geometric(field, z, upper - lower)[1])
        value = mult(value, geometric(field, u, coefficient_count)[1])
        value = mult(value, geometric(field, v, carry_count)[1])
        result = field.add(result, value)
    return result


def carry_points(base, modulus, carry_limit, digit, initial, final):
    for incoming in range(1 if initial else carry_limit):
        for coefficient in range(modulus):
            total = incoming + coefficient
            if total % base != digit:
                continue
            outgoing = total // base
            if final and outgoing:
                continue
            yield coefficient, incoming, outgoing


def transform(values, root, modulus):
    length = len(values)
    powers = [pow(root, exponent, modulus) for exponent in range(length)]
    return [sum(value * powers[index * frequency % length]
                for index, value in enumerate(values)) % modulus
            for frequency in range(length)]


def check_integer_transform(prime, digit_bits, first, second):
    base = 1 << digit_bits
    modulus = (prime ** (4 * digit_bits) + 1) // 2
    length = 2 * digit_bits
    root = prime ** 4 % modulus
    assert modulus % 2 == 1 and gcd(length, modulus) == 1
    assert pow(root, digit_bits, modulus) == modulus - 1
    assert modulus > 1 << (4 * digit_bits)
    fa = transform(first + [0] * digit_bits, root, modulus)
    fb = transform(second + [0] * digit_bits, root, modulus)
    transformed = transform([a * b % modulus for a, b in zip(fa, fb)],
                            pow(root, -1, modulus), modulus)
    inverse = pow(length, -1, modulus)
    actual = [value * inverse % modulus for value in transformed]
    expected = [sum(first[i] * second[j]
                    for i in range(digit_bits) for j in range(digit_bits)
                    if i + j == position) for position in range(length)]
    assert actual == expected
    carry_limit = 1 << modulus.bit_length()
    carry, reconstructed, place = 0, 0, 1
    for coefficient in actual:
        assert coefficient < modulus and carry < modulus < carry_limit
        total = coefficient + carry
        reconstructed += total % base * place
        carry = total // base
        place *= base
    assert carry == 0
    first_integer = sum(value * base ** index for index, value in enumerate(first))
    second_integer = sum(value * base ** index for index, value in enumerate(second))
    assert reconstructed == first_integer * second_integer


def verify():
    field = PrimePowerField(3, 8)
    assert field.size == 6561 and field.order == 6560
    orthogonality_cases = 0
    for order in (5, 8, 41, field.order):
        for first, second in ((0, 0), (1, 1), (0, 1), (3, 7)):
            raw = 0
            for frequency in range(order):
                raw = field.add(raw, field.root_phase(order, frequency * (first - second)))
            normalized = field.scale(raw, pow(order % field.prime, -1, field.prime))
            assert normalized == int((first - second) % order == 0)
            orthogonality_cases += 1
    assert 41 % 3 != 1  # Normalization cannot silently be omitted.

    inverse_gate_cases = 0
    modulus, dimension = 5, 3
    linear_map = ((1, 2, 0), (0, 1, 1), (1, 0, 1))
    for sample in range(12):
        beta = tuple((sample * (index + 1) + index) % modulus
                     for index in range(dimension))
        aligned = tuple(-sum(beta[row] * linear_map[row][column]
                             for row in range(dimension)) % modulus
                        for column in range(dimension))
        for changed in (False, True):
            alpha = ((aligned[0] + changed) % modulus,) + aligned[1:]
            raw = 0
            for setting in product(range(modulus), repeat=dimension):
                output = tuple(sum(row[index] * setting[index]
                                   for index in range(dimension)) % modulus
                               for row in linear_map)
                exponent = sum(a * u + b * v
                               for a, u, b, v in zip(alpha, setting, beta, output))
                raw = field.add(raw, field.root_phase(modulus, exponent))
            assert raw == (0 if changed else pow(modulus, dimension, field.prime))
            inverse_gate_cases += 1

    bilinear_cases = nonunit_cases = 0
    for modulus in (5, 8, 41):
        triples = (product(range(modulus), repeat=3) if modulus < 10 else
                   ((index * 7 % modulus, (index * 13 + 1) % modulus,
                     index * 19 % modulus) for index in range(32)))
        for a, b, c in triples:
            raw = 0
            for u in range(modulus):
                for v in range(modulus):
                    raw = field.add(raw, field.root_phase(modulus, a * u + b * v + c * u * v))
            common = gcd(c, modulus)
            if a % common or b % common:
                expected = 0
            else:
                reduced = modulus // common
                solution = (0 if reduced == 1 else
                            -(b // common) * pow(c // common, -1, reduced) % reduced)
                expected = field.scale(field.root_phase(modulus, a * solution),
                                       modulus * common)
            assert raw == expected
            bilinear_cases += 1
            nonunit_cases += common > 1

    carry_cases = carry_settings = 0
    modulus = 41
    frequencies = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                   (modulus - 1, field.order - 1, field.order - 1), (7, 23, 31)]
    for base, carry_limit in product((2, 4, 16), (57, 64)):
        for digit in range(base):
            for initial, final in product((False, True), repeat=2):
                settings = tuple(carry_points(base, modulus, carry_limit, digit, initial, final))
                for alpha, beta, gamma in frequencies:
                    raw = 0
                    for coefficient, incoming, outgoing in settings:
                        raw = field.add(raw, field.phase(
                            (field.order // modulus) * alpha * coefficient
                            + beta * incoming + gamma * outgoing))
                    expected = carry_closed(field, base, modulus, carry_limit,
                                            digit, alpha, beta, gamma, initial, final)
                    assert raw == expected
                    carry_cases += 1
                    carry_settings += len(settings)

    digit_cases = zero_factor_cases = 0
    width = 4
    for alpha, signs_code, free in product((0, 1, 7, 16), range(1 << width), range(width + 1)):
        signs = [-1 if (signs_code >> bit) & 1 else 1 for bit in range(width)]
        z = field.root_phase(modulus, alpha)
        for start in range(0, 1 << width, 1 << free):
            raw = [0] * (width + 1)
            for value in range(start, start + (1 << free)):
                sign = 1
                for bit in range(width):
                    if (value >> bit) & 1:
                        sign *= signs[bit]
                term = field.scale(field.power(z, value), sign)
                raw[0] = field.add(raw[0], term)
                for bit in range(width):
                    if (value >> bit) & 1:
                        raw[bit + 1] = field.add(raw[bit + 1], term)
            fixed_sign = 1
            for bit in range(free, width):
                if (start >> bit) & 1:
                    fixed_sign *= signs[bit]
            common = field.scale(field.power(z, start), fixed_sign)
            factors = [field.add(1, field.scale(field.power(z, 1 << bit), signs[bit]))
                       for bit in range(free)]
            scalar = common
            for value in factors:
                scalar = field.multiply(scalar, value)
            expected = [scalar]
            for bit in range(width):
                if bit >= free:
                    expected.append(scalar if (start >> bit) & 1 else 0)
                else:
                    value = field.multiply(common, field.scale(field.power(z, 1 << bit), signs[bit]))
                    for other, factor in enumerate(factors):
                        if other != bit:
                            value = field.multiply(value, factor)
                    expected.append(value)
            assert raw == expected
            digit_cases += 1
            zero_factor_cases += 0 in factors

    transform_cases = 0
    for prime in (3, 5):
        for digit_bits in (2, 8, 32):
            base = 1 << digit_bits
            if digit_bits == 2:
                cases = ((list(first), list(second))
                         for first in product(range(base), repeat=digit_bits)
                         for second in product(range(base), repeat=digit_bits))
            else:
                cases = (([(17 * index + 23 * sample + index * sample) % base
                           for index in range(digit_bits)],
                          [(31 * index + 7 * sample + sample * sample) % base
                           for index in range(digit_bits)])
                         for sample in range(8 if digit_bits == 8 else 4))
            for first, second in cases:
                check_integer_transform(prime, digit_bits, first, second)
                transform_cases += 1

    return {
        "scope": "exact odd-characteristic gate and modified integer-transform identities; no faster block contraction",
        "coefficient_field": {
            "prime": field.prime,
            "degree": field.degree,
            "modulus_polynomial_encoding": field.polynomial,
            "primitive_element": field.generator,
            "verified_cycle_length": field.order,
        },
        "normalized_character_orthogonality_cases": orthogonality_cases,
        "inverse_gate_cases_with_nontrivial_scalar": inverse_gate_cases,
        "bilinear_gauss_sum_cases": bilinear_cases,
        "nonunit_bilinear_cases": nonunit_cases,
        "general_carry_cases": carry_cases,
        "carry_settings_summed": carry_settings,
        "walsh_tagged_digit_cases": digit_cases,
        "zero_factor_digit_cases": zero_factor_cases,
        "modified_integer_transform_cases": transform_cases,
        "integer_transform_digit_bits": [2, 8, 32],
        "missing_normalization_negative_control": {"modulus": 41, "raw_scalar": 2, "required_scalar": 1},
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
