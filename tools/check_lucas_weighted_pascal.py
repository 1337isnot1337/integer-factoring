#!/usr/bin/env python3
"""Independent exact controls for gap-skipping weighted Pascal profiles."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random

from check_all_field_pascal_crossings import normalized_infinity, weighted_orders
from check_pascal_complement_convolution import modular_complement
from lucas_weighted_pascal import (greedy_complement, lucas_binomial,
                                  weighted_constant_profile, weighted_pascal_profile)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    lucas_checks = 0
    complements = 0
    for prime in (2, 3, 5, 7):
        for top in range(80):
            for bottom in range(top + 2):
                expected = comb(top, bottom) % prime if bottom <= top else 0
                assert lucas_binomial(top, bottom, prime) == expected
                lucas_checks += 1
        for size in range(1, 7):
            for modes in combinations(range(6), size):
                expected = modular_complement(modes, prime)
                assert list(greedy_complement(modes, prime)) == list(expected)
                complements += 1

    large = []
    for prime in (2, 3, 5, 7):
        for exponent in (32, 64, 128):
            high = prime ** exponent
            modes = (0, 1, high, high + 1)
            stats = {}
            actual = greedy_complement(modes, prime, statistics=stats)
            # A separate explicit basis: first two available low points,
            # then the corresponding translates by the high digit.
            expected = (2, 3, high + 2, high + 3)
            assert actual == expected, (prime, exponent, actual)
            assert stats['blocks_visited'] < 40 * prime * (exponent + 2)
            large.append({'prime': prime, 'high_digit_exponent': exponent,
                          'largest_index_bits': (high + 1).bit_length(), **stats})

    rng = random.Random(20261031)
    profiles = []
    for prime, modes in ((2, (1, 2)), (2, (1, 2, 4)), (2, (3, 5, 6)),
                         (3, (1, 3)), (3, (2, 4, 6)), (5, (1, 5)),
                         (2, (0, 1, 2)), (3, (0, 1, 2))):
        antichain = all(a == b or lucas_binomial(b, a, prime) == 0
                        for a in modes for b in modes)
        for _ in range(12):
            if antichain:
                heights = [rng.randrange(-8, 9) for _ in modes]
            else:
                slope = rng.randrange(1, 4)
                heights = [-slope * a for a in modes]
            expected, complement, _ = weighted_pascal_profile(modes, heights, prime)
            actual, precision = weighted_orders(normalized_infinity(modes, prime),
                                                heights + heights, [0] * (2 * len(modes)), prime)
            assert expected == actual
            profiles.append({'prime': prime, 'modes': list(modes), 'heights': heights,
                             'complement': list(complement), 'orders': list(actual),
                             'precision': precision})

    modes, heights, prime = (0, 1), (-2, -1), 2
    try:
        weighted_pascal_profile(modes, heights, prime)
    except ValueError:
        pass
    else:
        raise AssertionError('an irregular height vector was accepted')
    complement = greedy_complement(modes, prime)
    constant = [[lucas_binomial(c, a, prime) for c in complement] for a in modes]
    false_core, _ = weighted_constant_profile(constant, heights, complement, prime)
    false_prediction = tuple(sorted([a + h for a, h in zip(modes, heights)] + list(false_core)))
    actual, _ = weighted_orders(normalized_infinity(modes, prime), heights + heights, [0] * 4, prime)
    assert false_prediction != actual

    # Ramification multiplies the mode and complement exponents too.
    core = [[1, 1], [1, 0]]
    rational_core, _ = weighted_constant_profile(core, (Fraction(-3, 2), Fraction(5, 3)), (3, 5), 2)
    integral_core, _ = weighted_constant_profile(core, (-9, 10), (18, 30), 2)
    assert tuple(6 * x for x in rational_core) == integral_core

    report = {
        'status': 'Exact constructive profile theorem and finite controls; no post-half pressure witness.',
        'source_predecessor': '18cec06deac47c8f142ea67a60e2278ac2234d67',
        'lucas_coefficient_comparisons': lucas_checks,
        'independent_complement_comparisons': complements,
        'large_gap_cases': large,
        'independent_polynomial_profile_cases': profiles,
        'irregular_height_negative_control': {'modes': list(modes), 'heights': list(heights),
                                              'false_prediction': list(false_prediction),
                                              'actual_orders': list(actual)},
        'rational_ramification_control': True,
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'module_sha256': sha256(Path(__file__).with_name('lucas_weighted_pascal.py').read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'lucas_checks': lucas_checks, 'complement_cases': complements,
                      'large_gap_cases': len(large), 'weighted_profile_cases': len(profiles),
                      'negative_control': True}))


if __name__ == '__main__':
    main()
