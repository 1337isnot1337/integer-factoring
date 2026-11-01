#!/usr/bin/env python3
"""Exact rational natural-pressure certificate for the stabilized 24-letter code."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from affine_place_profile_oracle import primitive_seed
from check_implicit_reflection_codes import matrix
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_physical_code_gauges import start_params


def integer_digest(value):
    return sha256(value.to_bytes(max(1, (value.bit_length() + 7) // 8), 'big')).hexdigest()


def homogeneous_sum(orders, multiplier, numerator, denominator):
    minimum, maximum = min(orders), max(orders)
    degree = multiplier * (maximum - minimum)
    value = sum(count * pow(numerator, multiplier * (order - minimum))
                * pow(denominator, multiplier * (maximum - order))
                for order, count in Counter(orders).items())
    return value, degree, minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    weights = physical_bit_energies(16)
    energy = lambda symbol: sum(value for bit, value in enumerate(weights) if symbol >> bit & 1)
    letters = sorted(range(256), key=lambda symbol: (energy(symbol), symbol))[:24]
    raw = restricted_tensor(16, letters)
    seed, common = primitive_seed(raw)
    q, energies = len(letters), [energy(symbol) for symbol in letters]
    first = RectangularOracle(seed['values'], q, q).prefixes(start_params(np.asarray(energies), 1), 1)['orders']
    profiles = []
    for cut, infinity in ((1, True), (0, False)):
        if infinity:
            row = [energies[b] // 2 for a in range(q) for b in range(q)]
            column = row
        else:
            row = [(energies[a] + energies[b]) // 2 for a in range(q) for b in range(q)]
            column = [0] * (q * q)
        orders = MatrixOracle(matrix(seed, cut), q * q).profile(row, column, 1, infinity)['orders']
        assert all(value == round(value) for value in orders)
        expected = first[:q * q] if infinity else first[q * q:]
        assert sorted(orders) == sorted(expected)
        profiles.append([int(value) for value in sorted(orders)])
    # R=A/B=10/7, z=54/55, t=A log(1/z), and type weights z^(B E).
    A, B, p, denominator = 10, 7, 54, 55
    sums = [homogeneous_sum(orders, A, p, denominator) for orders in profiles]
    assert sum(row[2] for row in sums) == 0
    physical, physical_degree, physical_minimum = homogeneous_sum(energies, B, p, denominator)
    assert physical_minimum == 0
    left = pow(sums[0][0] * sums[1][0], B)
    right = pow(physical, 2 * A)
    exponent_difference = 2 * A * physical_degree - B * sum(row[1] for row in sums)
    if exponent_difference >= 0:
        left *= pow(denominator, exponent_difference)
    else:
        right *= pow(denominator, -exponent_difference)
    assert left < right
    type_counts = [pow(p, B * value) * pow(denominator, B * (max(energies) - value))
                   for value in energies]
    assert sum(type_counts) == physical and min(type_counts) > 0
    names = (Path(__file__).name, 'restrict_implicit_reflection.py', 'check_implicit_reflection_codes.py',
             'reflection_dual_rail_gaussian.py', 'affine_place_profile_oracle.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    result = {'scope': 'Exact local integer profiles and strict rational pressure certificate rho<5/7; this remains above the factoring threshold 1/sqrt(2).',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'interval': 16, 'labels': letters, 'energies': energies,
              'primitive_degree': seed['degree'], 'scalar_content': common,
              'tensor_sha256': digest(seed['values']),
              'integer_profiles_infinity_then_zero': profiles,
              'profile_channel_count': sum(map(len, profiles)),
              'rational_base': [p, denominator], 'support_parameter': [A, B],
              'balanced_parameter_strict_upper_bound': [A, 2 * B],
              'cleared_denominator_exponent_difference': exponent_difference,
              'left_integer_bits': left.bit_length(), 'right_integer_bits': right.bit_length(),
              'left_integer_sha256': integer_digest(left), 'right_integer_sha256': integer_digest(right),
              'strict_difference_sha256': integer_digest(right - left),
              'positive_type_counts': [str(value) for value in type_counts],
              'type_count_sum': str(physical), 'passed': True}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': True, 'degree': seed['degree'], 'profile_channels': result['profile_channel_count'],
                      'rho_upper': [A, 2 * B], 'integer_bits': [left.bit_length(), right.bit_length()]}))


if __name__ == '__main__':
    main()
