#!/usr/bin/env python3
"""Independent exact controls for the dyadic arithmetic representation."""

import argparse
import json
from math import isqrt
from pathlib import Path
import random

from dyadic_hyperbola import (affine_interval_count, divisor_count, floor_sum,
                              inverse_fibre, inverse_polynomial, mobius_section,
                              rectangle_count, recursive_rectangle_count,
                              smallest_factor)


def check():
    counts = dict(floor_sums=0, affine_intervals=0, inverse_fibre_entries=0,
                  rectangles=0, recursive_rectangles=0, mobius_section_entries=0,
                  polynomial_fibre_entries=0, divisor_prefixes=0, factorizations=0)
    rng = random.Random(20261031)
    for modulus in range(1, 33):
        for length in range(17):
            for _ in range(32):
                slope, offset = rng.randrange(-64, 65), rng.randrange(-64, 65)
                expected = sum((slope * i + offset) // modulus for i in range(length))
                assert floor_sum(length, modulus, slope, offset) == expected
                counts['floor_sums'] += 1
    for modulus in range(1, 18):
        for slope in range(modulus):
            for offset in range(modulus):
                for _ in range(8):
                    start, stop = sorted(rng.choices(range(modulus + 1), k=2))
                    lower, upper = sorted(rng.choices(range(modulus + 1), k=2))
                    expected = sum(lower <= (slope * t + offset) % modulus < upper
                                   for t in range(start, stop))
                    assert affine_interval_count(modulus, slope, offset, start,
                                                 stop, lower, upper) == expected
                    counts['affine_intervals'] += 1
    for bits in range(1, 10):
        modulus = 1 << bits
        for number in range(1, modulus, 2):
            for low_bits in range((bits + 1) // 2, bits + 1):
                base = 1 << low_bits
                high_modulus = modulus // base
                for root in range(1, base, 2):
                    b, c, slope = inverse_fibre(number, bits, low_bits, root)
                    for t in range(high_modulus):
                        x = root + base * t
                        y = b + base * ((c + slope * t) % high_modulus)
                        assert x * y % modulus == number
                        counts['inverse_fibre_entries'] += 1
            for _ in range(12):
                x0, x1 = sorted(rng.choices(range(modulus + 1), k=2))
                y0, y1 = sorted(rng.choices(range(modulus + 1), k=2))
                expected = sum(y0 <= number * pow(x, -1, modulus) % modulus < y1
                               for x in range(x0 + (1 - x0 % 2), x1, 2))
                for low_bits in {(bits + 1) // 2, bits}:
                    actual, _ = rectangle_count(number, bits, (x0, x1), (y0, y1),
                                                 low_bits)
                    assert actual == expected
                    counts['rectangles'] += 1
                assert recursive_rectangle_count(number, bits, (x0, x1),
                                                   (y0, y1)) == expected
                counts['recursive_rectangles'] += 1
    for bits in range(1, 8):
        modulus = 1 << bits
        for _ in range(32):
            matrix = (rng.randrange(1, modulus, 2), rng.randrange(modulus),
                      rng.randrange(0, modulus, 2))
            a, b, c = matrix
            for digit in (0, 1):
                out, (a1, b1, c1) = mobius_section(matrix, digit, bits)
                for t in range(modulus // 2):
                    x = digit + 2 * t
                    y = (a * x + b) * pow(c * x + 1, -1, modulus) % modulus
                    lifted = out + 2 * ((a1 * t + b1)
                                        * pow(c1 * t + 1, -1, modulus // 2)
                                        % (modulus // 2))
                    assert y == lifted
                    counts['mobius_section_entries'] += 1
        for number in range(1, modulus, 2):
            for low_bits in range(1, bits + 1):
                base, high_modulus = 1 << low_bits, modulus >> low_bits
                for root in range(1, base, 2):
                    b, coefficients = inverse_polynomial(number, bits, low_bits, root)
                    for t in range(high_modulus):
                        value = 0
                        for coefficient in reversed(coefficients):
                            value = (value * t + coefficient) % high_modulus
                        assert (root + base * t) * (b + base * value) % modulus == number
                        counts['polynomial_fibre_entries'] += 1
    for number in range(3, 1024, 2):
        for upper in range(2, isqrt(number) + 2):
            expected = sum(number % d == 0 for d in range(2, upper))
            assert divisor_count(number, 2, upper) == expected
            counts['divisor_prefixes'] += 1
    for number in range(2, 4097):
        expected = next((d for d in range(2, isqrt(number) + 1)
                         if number % d == 0), number)
        assert smallest_factor(number) == expected
        counts['factorizations'] += 1
    # A whole 2^64-element affine fibre is counted without expanding it.
    for number in (1, 3, 17, (1 << 127) - 1):
        for root in (1, 3, (1 << 64) - 1):
            _, c, slope = inverse_fibre(number, 128, 64, root)
            assert slope % 2
            assert affine_interval_count(1 << 64, slope, c, 0, 1 << 64,
                                         17, (1 << 63) + 19) == (1 << 63) + 2
    return {
        'scope': 'Exact identities and finite controls; no improved factoring exponent.',
        'seed': 20261031,
        'checks': counts,
        'large_full_fibres': 12,
        'complete_factorization_range': [2, 4096],
        'remaining_cost': 'Up to 2^(ceil(n/2)-1) low roots; polynomial bit cost per root.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
