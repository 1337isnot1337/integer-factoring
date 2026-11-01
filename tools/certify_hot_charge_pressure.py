#!/usr/bin/env python3
"""Exact rational pressure certificates from complete natural charge profiles."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import time

from hot_code_charge_profiles import hot_energy, natural_profiles


def integer_digest(value):
    return sha256(value.to_bytes(max(1, (value.bit_length() + 7) // 8), 'big')).hexdigest()


def homogeneous_sum(exponents, numerator, denominator):
    counts = Counter(exponents)
    degree = max(counts)
    value = sum(count * pow(numerator, exponent) * pow(denominator, degree - exponent)
                for exponent, count in counts.items())
    return value, degree


def exact_comparison(profile, numerator, denominator, support_numerator, support_denominator):
    """Certify P(z^A)^(2B) < (Q_X(z^B) Q_Y(z^B))^A.

    The squared comparison covers unequal paired codes. R=A/B is the
    support parameter and rho<R/2 is the resulting balanced bound.
    """
    if (not 0 < numerator < denominator or gcd(numerator, denominator) != 1
            or support_numerator <= 0 or support_denominator <= 0
            or gcd(support_numerator, support_denominator) != 1):
        raise ValueError('Use a reduced rational base in (0,1) and a positive reduced support parameter.')
    a, b = support_numerator, support_denominator
    pp, dp = homogeneous_sum([a * value for value in profile['orders_w']], numerator, denominator)
    qx, dx = homogeneous_sum([b * value for value in profile['energies_x']], numerator, denominator)
    qy, dy = homogeneous_sum([b * value for value in profile['energies_y']], numerator, denominator)
    left, right = pow(pp, 2 * b), pow(qx * qy, a)
    denominator_difference = 2 * b * dp - a * (dx + dy)
    if denominator_difference >= 0:
        right *= pow(denominator, denominator_difference)
    else:
        left *= pow(denominator, -denominator_difference)
    reduced = gcd(a, 2 * b)
    return {'passed': left < right, 'rational_base': [numerator, denominator],
            'support_parameter': [a, b], 'balanced_parameter_strict_upper_bound': [a // reduced, 2 * b // reduced],
            'comparison': 'P(z^A)^(2B) < (Q_X(z^B) Q_Y(z^B))^A',
            'cleared_denominator_exponent_difference': denominator_difference,
            'left_integer_bits': left.bit_length(), 'right_integer_bits': right.bit_length(),
            'left_integer_sha256': integer_digest(left), 'right_integer_sha256': integer_digest(right),
            'strict_difference_sha256': integer_digest(abs(right - left)),
            'positive_type_count_formula': 'numerator^first_exponent * denominator^second_exponent',
            'type_count_exponents_x': [[b * value, dx - b * value] for value in profile['energies_x']],
            'type_count_exponents_y': [[b * value, dy - b * value] for value in profile['energies_y']],
            'type_count_sum_x_hex': hex(qx), 'type_count_sum_y_hex': hex(qy)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=64)
    parser.add_argument('--modes', type=int, default=12)
    parser.add_argument('--base', type=int, nargs=2, default=[459, 460])
    parser.add_argument('--support', type=int, nargs=2, default=[64, 45])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 2 <= args.size <= 1 << args.modes or not 1 <= args.modes <= 16:
        raise ValueError('Use a nontrivial physical prefix in at most sixteen hot modes.')
    began = time.monotonic()
    labels = sorted(range(1 << args.modes), key=lambda value: (hot_energy(value), value))[:args.size]
    profile = natural_profiles(labels)
    comparison = exact_comparison(profile, *args.base, *args.support)
    if not comparison['passed']:
        raise ArithmeticError('The proposed strict rational pressure bound did not pass.')
    interval = 4
    while interval // 2 < profile['modes'] or 4 * interval <= max(profile['energies_x'] + profile['energies_y']):
        interval *= 2
    names = (Path(__file__).name, 'hot_code_charge_profiles.py', 'hot_code_charge_profiles.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'restrict_implicit_reflection.py')
    report = {'scope': 'Exact complete natural profiles and strict rational support bound; no factoring exponent improvement unless its separate threshold is passed.',
              'arguments': {'size': args.size, 'modes': args.modes, 'base': args.base, 'support': args.support},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'eligible_reflection_interval': interval, 'profile': profile, 'certificate': comparison,
              'elapsed_seconds': time.monotonic() - began, 'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'size': args.size,
                      'rho_strict_bound': comparison['balanced_parameter_strict_upper_bound'],
                      'integer_bits': [comparison['left_integer_bits'], comparison['right_integer_bits']],
                      'elapsed_seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()
