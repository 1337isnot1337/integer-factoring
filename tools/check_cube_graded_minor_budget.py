#!/usr/bin/env python3
"""Exact graded cube exclusion and a larger-alphabet profile target."""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import argparse
import json
from pathlib import Path
import random

from cube_graded_minor_budget import (
    cumulative, envelope, integer_hash, integer_relaxation,
    positive_interval_certificate, power, rational_blocks,
)


def compressed_controls():
    rng = random.Random(915071)
    count = 0
    for alphabet in range(2, 7):
        for _ in range(8):
            letters = [0] + [rng.randrange(1, 8) for _ in range(alphabet - 1)]
            result = envelope(Counter(letters))
            for wires, factor in ((1, Fraction(2, 3)), (2, Fraction(1)), (4, Fraction(2))):
                energies = sorted(sum(word) for word in product(letters, repeat=wires))
                total = 0
                for prefix, energy in enumerate(energies, 1):
                    total += energy
                    assert cumulative(result['blocks'], prefix) <= factor * total
                    count += 1
    binary = envelope(Counter({0: 1, 3: 1}))
    expanded = [int(order) for multiplicity, order in binary['blocks'] for _ in range(multiplicity)]
    assert expanded == [0, 2, 4, 6] + [12] * 7 + [18] * 4 + [24]
    return count


def four_letter_certificate():
    result = envelope(Counter({0: 1, 3: 1, 9: 1, 12: 1}))
    assert all(order.denominator == 1 for _, order in result['blocks'])
    partition = Counter({7 * int(order): multiplicity for multiplicity, order in result['blocks']})
    polynomial = power(partition, 3)
    for degree, value in power(Counter({0: 1, 24: 1, 72: 1, 96: 1}), 7).items():
        polynomial[degree] -= value
    certificate = positive_interval_certificate(polynomial)
    assert (certificate['leading_degree'], certificate['quotient_degree'],
            certificate['positive_coefficients'], certificate['minimum_positive_coefficient']) == (14, 2002, 2003, 3)
    assert 3 * 7 ** 2 > 12 ** 2
    return {'letter_heights': [0, 1, 3, 4], 'homogeneous_energy_multiplier': 3,
            'profile_blocks': rational_blocks(result), 'strict_pressure_floor': '7/12',
            'all_scale_certificate': certificate}


def large_abstract_witness():
    alphabet = 256
    result = envelope(Counter({0: 1, 3: alphabet - 1}))
    profile = integer_relaxation(result)
    expected = {0: 1, 2: 255, 4: 255, 6: 65025, 13: 262650,
                14: 62985, 18: 66325500, 24: 4228250625}
    assert dict(profile) == expected
    assert sum(profile.values()) == alphabet ** 4
    scale = Fraction(6, 13)
    occupation = Fraction(19, 20)
    pressure = Fraction(72, 125)
    partition = sum(count * scale ** order for order, count in profile.items())
    a, b = occupation.numerator, occupation.denominator
    c, d = pressure.numerator, pressure.denominator
    left = partition ** (b * d) * scale ** (-8 * a * d)
    entropy_base = Fraction(b ** b * (alphabet - 1) ** a,
                            a ** a * (b - a) ** (b - a))
    right = entropy_base ** (4 * c)
    lhs = left.numerator * right.denominator
    rhs = right.numerator * left.denominator
    assert lhs < rhs
    assert 3 * c ** 2 < d ** 2
    return {'alphabet': alphabet, 'nonzero_height': 1,
            'homogeneous_energy_multiplier': 3,
            'profile_order_multiplicity': sorted(profile.items()),
            'scale': [6, 13], 'nonzero_fraction': [19, 20],
            'nonzero_letter_distribution': 'uniform on all 255 nonzero letters',
            'strict_pressure_upper_bound': '72/125',
            'partition_numerator': str(partition.numerator),
            'partition_denominator': str(partition.denominator),
            'comparison_left_bits': lhs.bit_length(), 'comparison_right_bits': rhs.bit_length(),
            'comparison_left_sha256': integer_hash(lhs),
            'comparison_right_sha256': integer_hash(rhs),
            'satisfies_all_shared_prefix_budgets': True,
            'actual_tensor_constructed': False,
            'scope': 'Feasibility of the abstract necessary profile budgets only.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls = compressed_controls()
    four = four_letter_certificate()
    large = large_abstract_witness()
    result = {'status': 'passed', 'compressed_prefix_controls': controls,
              'four_letter_obstruction': four, 'large_abstract_profile': large,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in ('cube_graded_minor_budget.py', 'check_cube_graded_minor_budget.py')}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'compressed_prefix_controls': controls,
                      'q4_positive_coefficients': 2003, 'q4_floor': '7/12',
                      'q256_abstract_pressure': '72/125', 'actual_q256_tensor': False}, sort_keys=True))


if __name__ == '__main__':
    main()
