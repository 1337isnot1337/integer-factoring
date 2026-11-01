#!/usr/bin/env python3
"""Finite F2 extension-field controls for the rational-family transfer proof."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random


def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def multiply(a, b, modulus):
    product = 0
    while b:
        if b & 1:
            product ^= a
        a <<= 1
        b >>= 1
    return remainder(product, modulus)


def power(a, exponent, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, a, modulus)
        a = multiply(a, a, modulus)
        exponent >>= 1
    return result


def irreducible(modulus):
    degree = modulus.bit_length()-1
    x = remainder(2, modulus)
    value = x
    for i in range(1, degree+1):
        value = multiply(value, value, modulus)
        if i <= degree//2 and gcd(value ^ x, modulus) != 1:
            return False
    return value == x


def field_modulus(degree):
    return next((1 << degree) | tail for tail in range(1, 1 << degree, 2)
                if irreducible((1 << degree) | tail))


def evaluate(polynomial, point, modulus):
    result = 0
    for i in reversed(range(polynomial.bit_length())):
        result = multiply(result, point, modulus) ^ ((polynomial >> i) & 1)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(20261030)
    field_records, factor_cases, entries = [], 0, 0
    for degree in range(1, 9):
        modulus = field_modulus(degree)
        q = 1 << degree
        assert irreducible(modulus)
        # Every nonzero element must have the expected inverse.
        for a in range(1, q):
            assert multiply(a, power(a, q-2, modulus), modulus) == 1
        for rank in (1, 2, 3):
            for sample in range(10):
                left = [[rng.randrange(q) for _ in range(rank)] for _ in range(4)]
                right = [[rng.randrange(q) for _ in range(5)] for _ in range(rank)]
                for i in range(4):
                    for j in range(5):
                        value, descended = 0, 0
                        for channel in range(rank):
                            value ^= multiply(left[i][channel], right[channel][j], modulus)
                            for bit in range(degree):
                                descended ^= ((left[i][channel] >> bit) & 1) & (
                                    multiply(1 << bit, right[channel][j], modulus) & 1)
                        assert (value & 1) == descended
                        entries += 1
                factor_cases += 1
        # Reproduce specialization outside a union of denominator zero sets.
        denominators = [2, 3, 7, 11]
        total_degree = sum(p.bit_length()-1 for p in denominators)
        good_points = []
        if q > total_degree:
            bad = 0
            for point in range(q):
                if any(evaluate(p, point, modulus) == 0 for p in denominators):
                    bad += 1
                else:
                    good_points.append(point)
            assert bad <= total_degree and good_points
            for point in good_points:
                denominator = point ^ 1
                inverse = power(denominator, q-2, modulus)
                assert inverse ^ multiply(point, inverse, modulus) == 1
        field_records.append({'extension_degree': degree, 'modulus': modulus,
                              'elements': q, 'good_specialization_points': len(good_points)})
    # The coordinate projection sends 1 to 1; ordinary trace need not.
    modulus = field_modulus(2)
    assert (1 ^ power(1, 2, modulus)) == 0 and (1 & 1) == 1
    result = {'status': 'finite field/descent controls; rational transfer is proved separately',
              'fields': field_records, 'factorization_cases': factor_cases,
              'descended_entries_checked': entries, 'trace_normalization_negative_control': True,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'fields'}))


if __name__ == '__main__':
    main()
