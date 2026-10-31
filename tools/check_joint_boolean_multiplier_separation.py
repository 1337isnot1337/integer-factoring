#!/usr/bin/env python3
"""Exact joint-filter separation from every scalar Boolean multiplier on three parities."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from boolean_multiplier_profiles import (
    complete_multiplier_automaton, homogeneous_layer_patterns,
    layer_pattern_excess, multiply,
)


Q = Fraction(3, 4)
SINGLE = Fraction(3970975, 3981312)
JOINT = Fraction(1261129, 1327104)


def fraction_pair(value):
    return [value.numerator, value.denominator]


def square_root_lower(value, bits=100):
    root = isqrt((value.numerator << (2 * bits)) // value.denominator)
    lower = Fraction(root, 1 << bits)
    upper = Fraction(root + 1, 1 << bits)
    assert lower * lower <= value < upper * upper
    return lower


def profile_rows(automaton):
    rows = [{'kind': 'continuous', **row}
            for row in automaton['continuous_profiles']]
    for tail_index, tail in enumerate(automaton['singleton_tails']):
        for kind in ('startup', 'cycle'):
            for row in tail[kind]:
                rows.append({
                    'kind': kind, 'active_patterns': (tail['pattern'],),
                    'outside_exponents': tail['outside_exponents'],
                    'divergence': row['divergence'],
                    'degree_offset': row['degree_offset'],
                    'period': tail['period'] if kind == 'cycle' else 0,
                    'tail_index': tail_index,
                })
    return rows


def exact_piece_lower(row, distribution, support):
    active = set(row['active_patterns'])
    a = sum((distribution[h] * Q ** (-row['outside_exponents'][h])
             for h in range(16) if h not in active), Fraction(0))
    b = sum((distribution[h] * Q ** (-(h.bit_length() - 1))
             for h in active), Fraction(0))
    offset = row.get('degree_offset', 0)
    constant = a * (1 - support) + b * support * Q ** (-offset)
    u, v = a * support * Q ** (-offset), b * (1 - support)

    def at(position):
        x = Q ** (-position)
        return constant + u * x + v / x

    if row['kind'] == 'continuous':
        lower = row['lower_degree']
        upper = row['upper_degree']
        if u * Q ** (-2 * lower) >= v:
            return at(lower), 'lower_endpoint', lower
        if upper is not None and u * Q ** (-2 * upper) <= v:
            return at(upper), 'upper_endpoint', upper
        if u == 0:
            assert upper is None
            return constant, 'infinite_limit', None
        return constant + 2 * square_root_lower(u * v), 'square_root_lower', None

    h, period = row['divergence'], row['period']
    if period == 0:
        return at(h), 'startup', h
    if u == 0:
        return constant, 'infinite_limit', None
    index = 0
    while u * Q ** (-2 * (h + period * index)) < v:
        index += 1
    choices = [(at(h + period * k), h + period * k)
               for k in {max(0, index - 1), index}]
    value, position = min(choices)
    return value, 'arithmetic_progression', position


def joint_polynomials():
    coefficients = [1 + (1 << 2) + (1 << 4), 1 + 2 + 4, 1, 1 + 2 + 8]
    original = [coefficients[k.bit_count()] for k in range(8)]
    beta = []
    for k in range(8):
        value, subset = 0, k
        while True:
            power = k.bit_count() - subset.bit_count()
            c = sum(1 << a for a in range(power + 1) if a & ~power == 0)
            value ^= multiply(c, original[subset])
            if subset == 0:
                break
            subset = (subset - 1) & k
        beta.append(value)
    orders = [(value & -value).bit_length() - 1 for value in beta]
    degrees = [value.bit_length() - 1 for value in original]
    assert orders == [0, 3, 3, 6, 3, 6, 6, 4]
    assert degrees == [4, 2, 2, 0, 2, 0, 0, 3]
    assert all(value.bit_count() % 2 for value in original + beta)
    zero = sum((Q ** (orders[k] - k.bit_count()) for k in range(8)), Fraction(0)) / 8
    infinity = sum((Q ** (-degree) for degree in degrees), Fraction(0)) / 8
    assert zero * infinity == JOINT
    return {'coefficient_bit_polynomials': original, 'beta_bit_polynomials': beta,
            'beta_orders': orders, 'infinity_degrees': degrees,
            'zero_ratio': fraction_pair(zero), 'infinity_ratio': fraction_pair(infinity)}


def check():
    automaton = complete_multiplier_automaton(3)
    rows = profile_rows(automaton)
    assert len(rows) == 42
    branches, supports, kinds = [], [], Counter()
    global_lower = Fraction(1)
    minimizers = []
    for table in range(1, 255):
        truth = [(table >> k) & 1 for k in range(8)]
        patterns = homogeneous_layer_patterns(truth, 3)
        support = Fraction(sum(truth), 8)
        distribution = [Fraction(patterns.count(h), 8) for h in range(16)]
        local_lower = Fraction(1)
        for index, row in enumerate(rows):
            value, kind, position = exact_piece_lower(row, distribution, support)
            assert value >= SINGLE, (table, index, value)
            local_lower = min(local_lower, value)
            global_lower = min(global_lower, value)
            kinds[kind] += 1
            record = [table, index, kind, position, fraction_pair(value)]
            branches.append(record)
            if value == SINGLE:
                minimizers.append(record)
        for row in automaton['fixed_profiles']:
            value = (sum((distribution[h] * Q ** row['excesses'][h]
                          for h in range(16)), Fraction(0))
                     * (1 - support + support * Q ** (-row['degree'])))
            assert value >= SINGLE
            local_lower = min(local_lower, value)
        supports.append({'truth_table': table, 'support_size': sum(truth),
                         'certified_lower': fraction_pair(local_lower)})
    assert global_lower == SINGLE and minimizers
    multiplier = (1 << 2) | (1 << 4)
    pattern = homogeneous_layer_patterns([1] + [0] * 7, 3)
    attained = sum((Q ** layer_pattern_excess(h, multiplier)
                    for h in pattern), Fraction(0)) / 8
    attained *= Fraction(7, 8) + Fraction(1, 8) * Q ** (-4)
    assert attained == SINGLE
    witness = joint_polynomials()
    epsilon = Fraction(1, 7 ** 5)
    separation = (1 - epsilon) ** 6 * SINGLE - (1 + epsilon) ** 6 * JOINT
    assert separation > 0
    assert (1 - epsilon) ** 6 * JOINT > Fraction(9, 10)
    assert Fraction(9, 10) * Fraction(7, 4) ** 16 > 1
    digest = sha256(json.dumps(branches, separators=(',', ':')).encode()).hexdigest()
    return {
        'status': 'Exact all-support, all-degree, all-extension and all-ramification separation on three parity groups; no factoring improvement.',
        'q': fraction_pair(Q), 'nonconstant_boolean_supports': 254,
        'unbounded_envelope_pieces_per_support': len(rows),
        'exact_piece_certificates': len(branches), 'piece_kinds': dict(kinds),
        'piece_certificate_sha256': digest,
        'automaton_sha256': sha256(json.dumps(automaton, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        'uniform_single_multiplier_minimum': fraction_pair(SINGLE),
        'uniform_joint_ratio': fraction_pair(JOINT),
        'uniform_minimizers': minimizers, 'support_lower_bounds': supports,
        'joint_literal_witness': witness,
        'minimum_physical_length': 16,
        'maximum_three_group_parity_bias': fraction_pair(epsilon),
        'physical_separation_lower': fraction_pair(separation),
        'critical_pressure_positive_lower_after_exponentiation_at_n16': fraction_pair(Fraction(9, 10) * Fraction(7, 4) ** 16),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in (
        'exact_piece_certificates', 'piece_kinds',
        'uniform_single_multiplier_minimum', 'uniform_joint_ratio',
        'minimum_physical_length', 'physical_separation_lower')}))


if __name__ == '__main__':
    main()
