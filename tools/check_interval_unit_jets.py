#!/usr/bin/env python3
"""Literal and long-range recurrence controls for compressed unit jets."""

import argparse
import json
from pathlib import Path
import random

from interval_unit_jets import QuadraticJets, interval_elementaries, interval_product_jet


def literal(ring, interval, argument, number):
    result = ring.constant(1)
    for a in range(*interval):
        factor = (ring.add(argument, ring.constant(-a)) if number is None
                  else ring.add(ring.constant(number), ring.scale(argument, -a)))
        result = ring.multiply(result, factor)
    return result


def check():
    rng = random.Random(20261034)
    elementary_checks = literal_checks = inverse_checks = 0
    for _ in range(256):
        lower = rng.randrange(12)
        upper = lower + rng.randrange(17)
        degree = rng.randrange(21)
        expected = [1] + [0] * degree
        for a in range(lower, upper):
            for k in range(degree, 0, -1):
                expected[k] += a * expected[k - 1]
        assert interval_elementaries(lower, upper, degree) == tuple(expected)
        elementary_checks += 1
    for precision in range(1, 13):
        for degree in range(6):
            ring = QuadraticJets(precision, degree)
            w = ring.constant(0, 1)
            assert ring.add(ring.add(ring.multiply(w, w), w), ring.constant(1)) == ring.constant()
            for _ in range(8):
                argument = ring.argument()
                argument = ring.add(argument, ring.constant(rng.randrange(17)))
                assert ring.multiply(argument, ring.inverse(argument)) == ring.constant(1)
                inverse_checks += 1
                lower = rng.randrange(1, 17)
                interval = (lower, lower + rng.randrange(25))
                number = 2 * rng.randrange(128) + 1
                for selected in (None, number):
                    actual = interval_product_jet(ring, interval, argument, selected)
                    assert actual == literal(ring, interval, argument, selected)
                    literal_checks += 1
    large = []
    for precision, degree, input_bits, length_bits in ((32, 6, 120, 110),
                                                       (48, 8, 240, 200)):
        ring = QuadraticJets(precision, degree)
        x = ring.argument()
        number = (1 << (input_bits + 7)) - 1
        lower = (1 << input_bits) + 7
        upper = lower + (1 << length_bits) + 12
        interval = (lower, upper)
        middle = (lower + upper) // 2
        p = interval_product_jet(ring, interval, x)
        q = interval_product_jet(ring, interval, x, number)
        for selected, value in ((None, p), (number, q)):
            split = ring.multiply(interval_product_jet(ring, (lower, middle), x, selected),
                                  interval_product_jet(ring, (middle, upper), x, selected))
            assert split == value
        shifted = interval_product_jet(ring, interval, ring.add(x, ring.constant(1)))
        assert ring.multiply(ring.add(x, ring.constant(-upper + 1)), shifted) == ring.multiply(
            ring.add(x, ring.constant(-lower + 1)), p)
        denominator = ring.add(x, ring.constant(number))
        moved = ring.scale(ring.multiply(x, ring.inverse(denominator)), number)
        q_moved = interval_product_jet(ring, interval, moved, number)
        left = ring.multiply(ring.power(denominator, upper - lower),
                             ring.multiply(ring.add(ring.constant(number),
                                                     ring.scale(x, 1 - upper)), q_moved))
        right = ring.scale(ring.multiply(ring.add(ring.constant(number),
                                                 ring.scale(x, 1 - lower)), q),
                            pow(number, upper - lower, ring.modulus))
        assert left == right
        large.append({'precision': precision, 'jet_degree': degree,
                      'endpoint_bits': upper.bit_length(),
                      'interval_length': str(upper - lower),
                      'number': str(number), 'P_jet': p, 'Q_jet': q,
                      'split_and_boundary_checks': 4})
    return {
        'scope': 'Polynomial bit-cost local jets of long interval products; global gcd extraction remains open.',
        'seed': 20261034,
        'elementary_symmetric_controls': elementary_checks,
        'literal_product_jet_controls': literal_checks,
        'unit_inverse_controls': inverse_checks,
        'large_cases': large,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'large_cases'}, indent=2))
    print('Passed two long-range jet cases, with four split/boundary identities each.')
