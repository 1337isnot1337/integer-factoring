#!/usr/bin/env python3
"""Dense-content and long-interval controls for normalized integer jets."""

import argparse
import json
from pathlib import Path
import random

from interval_orbit_gcd import multiply_linear
from interval_unit_jets import QuadraticJets
from normalized_interval_jets import normalized_product_jet, valuation_two


def dense_check(interval, center, precision, degree, number):
    receipt = normalized_product_jet(interval, center, precision, degree, number)
    scale = 1 << receipt['scale_bits']
    polynomial = (1,)
    for a in range(*interval):
        constant, linear = ((center - a, scale) if number is None
                             else (number - a * center, -a * scale))
        polynomial = multiply_linear(polynomial, constant, linear)
    order = next(i for i, a in enumerate(polynomial) if a)
    valuation = min(valuation_two(a) for a in polynomial if a)
    assert (order, valuation) == (receipt['zero_order'], receipt['valuation'])
    expected = [((polynomial[i + order] >> valuation) % (1 << precision)
                 if i + order < len(polynomial) else 0) for i in range(degree + 1)]
    assert receipt['unit'] == expected


def odd_factorial_part(number, precision):
    """Independent period-table method, used only at small fixed precision."""
    modulus = 1 << precision
    block = 1
    for a in range(1, modulus, 2):
        block = block * a % modulus
    result = 1
    while number:
        quotient, remainder = divmod(number, modulus)
        result = result * pow(block, quotient, modulus) % modulus
        for a in range(1, remainder + 1, 2):
            result = result * a % modulus
        number //= 2
    return result


def split_check(interval, center, precision, degree, number):
    whole = normalized_product_jet(interval, center, precision, degree, number)
    middle = sum(interval) // 2
    parts = [normalized_product_jet(part, center, precision, degree, number)
             for part in ((interval[0], middle), (middle, interval[1]))]
    ring = QuadraticJets(precision, degree)
    value, valuation, order = ring.constant(1), 0, 0
    for part in parts:
        delta = whole['scale_bits'] - part['scale_bits']
        assert delta >= 0
        valuation += part['valuation'] + delta * part['zero_order']
        order += part['zero_order']
        adjusted = tuple((a * pow(2, delta * i, ring.modulus) % ring.modulus, 0)
                         for i, a in enumerate(part['unit']))
        value = ring.multiply(value, adjusted)
    assert (valuation, order) == (whole['valuation'], whole['zero_order'])
    assert [a for a, _ in value] == whole['unit']
    return whole


def check():
    rng = random.Random(20261036)
    dense_checks = 0
    for precision in range(1, 13):
        for degree in range(6):
            for _ in range(8):
                lower = rng.randrange(1, 17)
                interval = (lower, lower + rng.randrange(25))
                center = rng.randrange(-8, 41)
                for number in (None, rng.randrange(1, 257)):
                    dense_check(interval, center, precision, degree, number)
                    dense_checks += 1
    dense_check((7, 20), 11, 16, 8, 143)
    large = []
    for endpoint_bits, length_bits, precision, degree in ((120, 100, 8, 4),
                                                          (240, 200, 12, 5)):
        lower = (1 << endpoint_bits) + 7
        length = (1 << length_bits) + 12
        interval = (lower, lower + length)
        for number, center in ((None, lower - 1), (None, lower),
                               ((1 << (endpoint_bits + 7)) - 1, lower)):
            receipt = split_check(interval, center, precision, degree, number)
            if number is None and center == lower - 1:
                assert receipt['valuation'] == length - length.bit_count()
                expected = odd_factorial_part(length, precision)
                if length % 2:
                    expected = -expected % (1 << precision)
                assert receipt['unit'][0] == expected
            if number is not None and not receipt['zero_order']:
                valuation = 0
                bound = max(abs(number - a * center) for a in (lower, interval[1] - 1)).bit_length()
                for j in range(1, bound):
                    modulus = 1 << j
                    root = number * pow(center, -1, modulus) % modulus
                    valuation += -((root - interval[1]) // modulus) + ((root - lower) // modulus)
                assert valuation == receipt['valuation']
            saved = dict(receipt)
            saved['interval'] = [str(a) for a in saved['interval']]
            for key in ('center', 'number', 'valuation'):
                if saved[key] is not None:
                    saved[key] = str(saved[key])
            large.append(saved)
    return {'scope': 'Exact binary valuation and normalized scaled jets; no global root search.',
            'seed': 20261036, 'dense_content_and_jet_controls': dense_checks + 1,
            'large_split_controls': 6, 'large_factorial_unit_controls': 2,
            'large_cumulative_valuation_controls': 2,
            'large_cases': large}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'large_cases'}, indent=2))
