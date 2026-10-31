#!/usr/bin/env python3
"""Independent controls for the exact interval-gcd factoring interface."""

import argparse
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path
import random

from interval_orbit_gcd import interval_factor_gcd, interval_product, multiply_linear


def verify_rectangle(number, x_interval, y_interval):
    actual = interval_factor_gcd(number, x_interval, y_interval)
    expected = (1,)
    roots = []
    for d in range(*x_interval):
        if number % d == 0 and y_interval[0] <= number // d < y_interval[1]:
            roots.append(d)
            expected = multiply_linear(expected, -d, 1)
    assert actual == expected
    assert max(abs(c) for c in actual) <= (number + 1) ** len(roots)
    return len(roots)


def check():
    rng = random.Random(20261033)
    random_checks = box_checks = recurrence_checks = 0
    for _ in range(1024):
        number = rng.randrange(1, 513)
        x_interval = tuple(sorted(rng.choices(range(1, 18), k=2)))
        y_interval = tuple(sorted(rng.choices(range(1, 18), k=2)))
        verify_rectangle(number, x_interval, y_interval)
        random_checks += 1
        if random_checks <= 128:
            def evaluate(polynomial, value):
                result = 0
                for coefficient in reversed(polynomial):
                    result = result * value + coefficient
                return result
            x = Fraction(3, 2)
            p = interval_product(x_interval)
            assert ((x - x_interval[1] + 1) * evaluate(p, x + 1)
                    == (x - x_interval[0] + 1) * evaluate(p, x))
            q = interval_product(y_interval, number)
            degree = y_interval[1] - y_interval[0]
            moved = number * x / (number + x)
            assert ((number + x) ** degree * (number - (y_interval[1] - 1) * x)
                    * evaluate(q, moved)
                    == number ** degree * (number - (y_interval[0] - 1) * x)
                    * evaluate(q, x))
            recurrence_checks += 2
    for number in range(2, 257):
        count, power = 0, 2
        while power <= isqrt(number):
            x_interval = (power, min(2 * power, isqrt(number) + 1))
            y_interval = (number // (2 * power) + 1, number // power + 1)
            count += verify_rectangle(number, x_interval, y_interval)
            power *= 2
            box_checks += 1
        assert count == sum(number % d == 0 for d in range(2, isqrt(number) + 1))
    # A small-prime reduction can introduce a common root with no integer factor.
    assert interval_factor_gcd(35, (2, 5), (2, 9)) == (1,)
    assert (35 - 2 * 4) % 3 == 0
    return {
        'scope': 'Exact rational gcd identity and dense finite reference; no compressed gcd algorithm.',
        'seed': 20261033,
        'random_rectangle_checks': random_checks,
        'complete_input_range': [2, 256],
        'dyadic_factor_boxes': box_checks,
        'orbit_boundary_recurrence_controls': recurrence_checks,
        'small_prime_false_common_root_control': {'number': 35, 'prime': 3, 'x': 2, 'y': 4},
        'unresolved': 'Uniform deterministic elimination from binary interval endpoints without expanding input products.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
