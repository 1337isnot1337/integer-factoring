#!/usr/bin/env python3
"""Literal controls for complete unipotent-family interval aggregation."""

import argparse
import json
from pathlib import Path
import random

from dyadic_orbit_aggregation import affine_residue_pairs, unipotent_interval_sum


def literal(bits, depth, x_interval, y_interval, alpha, beta, gamma, delta):
    modulus = 1 << bits
    total = 0
    for k in range(1 << (bits - depth)):
        for t in range(*x_interval):
            z = (beta * t + gamma) % modulus
            y = (alpha * z * pow(1 + (1 << depth) * k * z, -1, modulus)
                 + delta) % modulus
            total += y_interval[0] <= y < y_interval[1]
    return total


def check():
    rng = random.Random(20261032)
    pair_checks = orbit_checks = 0
    for bits in range(7):
        modulus = 1 << bits
        for _ in range(128):
            a = 2 * rng.randrange(-modulus, modulus + 1) + 1
            b = rng.randrange(-modulus, modulus + 1)
            x0, x1 = sorted(rng.choices(range(3 * modulus + 1), k=2))
            y0, y1 = sorted(rng.choices(range(3 * modulus + 1), k=2))
            expected = sum((y - a * x - b) % modulus == 0
                           for x in range(x0, x1) for y in range(y0, y1))
            assert affine_residue_pairs(modulus, a, b, (x0, x1), (y0, y1)) == expected
            pair_checks += 1
    for bits in range(1, 9):
        modulus = 1 << bits
        for depth in range(1, bits + 1):
            for _ in range(48):
                x_interval = tuple(sorted(rng.choices(range(modulus + 1), k=2)))
                y_interval = tuple(sorted(rng.choices(range(modulus + 1), k=2)))
                alpha, beta = (rng.randrange(1, modulus, 2) for _ in range(2))
                gamma, delta = (rng.randrange(modulus) for _ in range(2))
                expected = literal(bits, depth, x_interval, y_interval,
                                   alpha, beta, gamma, delta)
                actual = unipotent_interval_sum(bits, depth, x_interval, y_interval,
                                                alpha, beta, gamma, delta)
                assert actual == expected, (bits, depth, actual, expected)
                orbit_checks += 1
    # Exhaust every interval pair for the two smallest nontrivial groups.
    for bits, depth in ((2, 1), (3, 1), (3, 2)):
        modulus = 1 << bits
        intervals = [(a, b) for a in range(modulus + 1)
                     for b in range(a, modulus + 1)]
        for x_interval in intervals:
            for y_interval in intervals:
                assert unipotent_interval_sum(bits, depth, x_interval, y_interval) == literal(
                    bits, depth, x_interval, y_interval, 1, 1, 0, 0)
                orbit_checks += 1
    large_checks = 0
    for bits, depth in ((64, 1), (128, 16), (256, 7)):
        modulus, size = 1 << bits, 1 << (bits - depth)
        interval = (modulus // 7, modulus // 3 + 17)
        for x_interval, y_interval in ((interval, (0, modulus)),
                                        ((0, modulus), interval)):
            actual = unipotent_interval_sum(bits, depth, x_interval, y_interval,
                                            7, 13, 30, 51)
            assert actual == size * (interval[1] - interval[0])
            large_checks += 1
    return {
        'scope': 'Uniform complete cyclic-family sums, with affine sandwiches; no factoring improvement.',
        'seed': 20261032,
        'affine_residue_pair_controls': pair_checks,
        'literal_orbit_interval_controls': orbit_checks,
        'large_permutation_controls': large_checks,
        'largest_group_size': '2^249',
        'remaining_obligation': 'Actual inverse sections have correlated varying affine parameters; uniform orbit sums alone do not aggregate them.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
