#!/usr/bin/env python3
"""Exact finite checks of the complete Boolean Cauchy Smith profile."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path
import random

from boolean_cauchy_operator import (
    boundary_rows, inverse_group_rank, kernel_matrices, monomial_operator,
    monomial_to_value_operator, multiply_matrices, normalized_inverse_mod_two,
    rotate_subset, smith_orders, unit_inverse, value_operator,
)
from cauchy_packet_polynomials import modular_packet_polynomial


def check():
    precision = 24
    modulus = 1 << precision
    records = []
    normalized_entries = 0
    factorization_entries = 0
    smith_profiles = 0
    inverse_products = 0
    contrasts = 0
    group_ranks = 0
    permanent_coefficients = 0
    rng = random.Random(29061)
    cases = [(length, center, point) for length in (2, 4)
             for center in range(-length, length) if center for point in (0, 1, 2, 3)]
    cases += [(8, center, 1) for center in (-8, -3, -1, 1, 3)] + [(8, 1, 0), (8, 1, 2)]
    for length, center, point in cases:
        count = 1 << length
        boundary = sum(1 << x for x in boundary_rows(length, center))
        b = abs(center)
        original, normalized = kernel_matrices(length, center, point, precision)
        for x in range(length):
            for y in range(length):
                assert normalized[x][y] % 2 == int(y == (x + center) % length)
                normalized_entries += 1
        operator = monomial_operator(original, modulus)
        unit = monomial_operator(normalized, modulus)
        powers = [(subset & boundary).bit_count() for subset in range(count)]
        for subset in range(count):
            factor = pow(2 * point, powers[subset], modulus)
            expected_column = rotate_subset(subset, length, center)
            for target in range(count):
                assert unit[subset][target] % 2 == int(target == expected_column)
                assert operator[subset][target] == factor * unit[subset][target] % modulus
                factorization_entries += 1
        actual_values = value_operator(original, modulus)
        assert actual_values == monomial_to_value_operator(operator, modulus)
        if point:
            step = 1 + ((point & -point).bit_length() - 1)
            expected_orders = sorted(step * power for power in powers)
        else:
            expected_orders = [0] * (1 << (length - b)) + [None] * (count - (1 << (length - b)))
        if length <= 4:
            assert smith_orders(actual_values, precision) == expected_orders
            smith_profiles += 1
            inverse_unit = unit_inverse(unit, modulus)
            normalized_inverse = [[value * pow(2 * point, b - powers[column], modulus) % modulus
                                   for column, value in enumerate(row)] for row in inverse_unit]
            expected_product = [[pow(2 * point, b, modulus) * int(i == j)
                                 for j in range(count)] for i in range(count)]
            assert multiply_matrices(operator, normalized_inverse, modulus) == expected_product
            assert multiply_matrices(normalized_inverse, operator, modulus) == expected_product
            inverse_products += 2
            value_inverse = monomial_to_value_operator(normalized_inverse, modulus)
            bit_rows = [sum((value & 1) << column for column, value in enumerate(row)) for row in value_inverse]
            assert bit_rows == normalized_inverse_mod_two(length, center)
        for _ in range(3):
            row = [rng.randrange(-3, 8) for _ in range(count)]
            top_coefficient = sum((-1)**(length - state.bit_count()) * row[state] for state in range(count))
            character = [(-1)**state.bit_count() for state in range(count)]
            observed = sum(row[y] * actual_values[y][x] * character[x]
                           for y in range(count) for x in range(count)) % modulus
            # The full-squarefree coefficient of the monomial operator is
            # independently built by Boolean finite differences above.
            permanent = operator[-1][-1]
            assert observed == (-1)**length * top_coefficient * permanent % modulus
            contrasts += 1
        records.append({'L': length, 't': center, 'zeta': point, 'boundary_rows': b,
                        'orders': {str(key): value for key, value in Counter(expected_orders).items()}})

    for length in (2, 4, 8):
        count = 1 << length
        character = [(-1)**state.bit_count() for state in range(count)]
        for center in range(-length, length):
            if not center:
                continue
            b = abs(center)
            polynomial = modular_packet_polynomial(length, center, b + 1, 0, character)
            expected = [0] * (b + 1)
            expected[b] = 1 << b
            while len(polynomial) > 1 and not polynomial[-1]:
                polynomial.pop()
            assert polynomial == expected, (length, center, polynomial)
            permanent_coefficients += length * (2 * length - 1) + 1
            if b < length:
                assert inverse_group_rank(length, center) == 1 << (length - b)
                group_ranks += 1

    root = Path(__file__).resolve().parent.parent
    sources = ['tools/boolean_cauchy_operator.py', 'tools/check_boolean_cauchy_smith.py',
               'tools/cauchy_packet_polynomials.py', 'research/boolean-cauchy-smith-profile.md']
    return {'status': 'pass', 'operator_cases': len(records),
            'normalized_kernel_entries_checked': normalized_entries,
            'boolean_factorization_entries_checked': factorization_entries,
            'independent_smith_profiles': smith_profiles,
            'normalized_inverse_products': inverse_products,
            'arbitrary_row_character_controls': contrasts,
            'permanent_polynomial_coefficients_checked': permanent_coefficients,
            'source_group_crossing_ranks_checked': group_ranks,
            'records': records,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in sources},
            'scope': 'Finite exact controls for an all-size Boolean substitution Smith theorem, '
                     'parity-character identity, normalized inverse, and specified group-cut '
                     'obstruction. The inverse-entry complexity proof is separate; no scalar '
                     'decoder or improved factoring score is claimed.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.write_text(json.dumps({'status': 'fail', 'error': repr(error)}, indent=2) + '\n')
        raise
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
