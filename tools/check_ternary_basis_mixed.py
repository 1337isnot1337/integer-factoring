#!/usr/bin/env python3
"""Exact basis, local profiles, and pressure for a ternary mixed crossing.

This verifies an ordinary uniform construction at rho < 821/1000.
It does not improve the established factoring exponent one half.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import local_polynomial, mixing_seed, rank
from check_weighted_multipoint_hermite import (
    flatten_polynomials, smith_orders, typed_power, type_words,
)


def inverse(matrix, prime=3):
    n = len(matrix)
    rows = [[int(x) % prime for x in row] + [int(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scalar = pow(rows[j][j], -1, prime)
        rows[j] = [scalar * value % prime for value in rows[j]]
        for i in range(n):
            if i != j:
                scalar = rows[i][j]
                rows[i] = [(a - scalar * b) % prime for a, b in zip(rows[i], rows[j])]
    assert all(rows[i][j] == int(i == j) for i in range(n) for j in range(n))
    return [row[n:] for row in rows]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def transform_seed(basis_a, basis_b):
    matrices = (basis_a, basis_b, transpose(inverse(basis_b)), transpose(inverse(basis_a)))
    source, output = mixing_seed(), {}
    for key in product(range(3), repeat=4):
        value = [0, 0, 0]
        for old_key, old_value in source.items():
            scalar = 1
            for leg in range(4):
                scalar *= matrices[leg][key[leg]][old_key[leg]]
            for exponent, coefficient in enumerate(old_value):
                value[exponent] = (value[exponent] + scalar * coefficient) % 3
        if any(value):
            output[key] = value
    for x, y, u, v in product(range(3), repeat=4):
        assert sum(output.get((x, y, u, v), ())) % 3 == int(x == v and y == u)
    return output


def ramified_profile(seed, point, cut, slopes):
    """Independent complete Toeplitz ranks after t=tau^2 and integer gauges."""
    twice_slopes = [int(2 * slope) for slope in slopes]
    assert all(2 * slope == value for slope, value in zip(slopes, twice_slopes))
    terms = {}
    for key, value in seed.items():
        local = local_polynomial(value, point, 2, 3)
        height = sum(key[leg] * twice_slopes[leg] for leg in range(4))
        terms[key] = {2 * j + height: a for j, a in enumerate(local) if a}
    minimum = min(j for terms_at_entry in terms.values() for j in terms_at_entry)
    clearing = max(0, -minimum)
    table = {}
    for key, values in terms.items():
        if not values:
            continue
        coefficients = [0] * (max(values) + clearing + 1)
        for j, value in values.items():
            coefficients[j + clearing] = value
        table[key] = coefficients
    orders = smith_orders(flatten_polynomials(table, (3,) * 4, cut), 3)
    return [j - clearing for j in orders], clearing


def typed_controls(seed, profiles):
    counts = (1, 1, 1)
    words = type_words(counts)
    table, dimensions = typed_power(seed, 3, 3, 3, counts, counts)
    for x, y, u, v in product(range(6), repeat=4):
        assert sum(table.get((x, y, u, v), ())) % 3 == int(x == v and y == u)
    controls = []
    for record in profiles:
        point, cut = record['point'], record['cut']
        slopes = [Fraction(*pair) for pair in record['slopes']]
        shift_twice = sum(2 * slope for slope in slopes)
        assert shift_twice.denominator == 1
        # Every retained physical word has sum of symbols exactly three.
        for keys in product(words, repeat=4):
            actual = sum(2 * slopes[leg] * sum(keys[leg]) for leg in range(4))
            assert actual == 3 * shift_twice
        local = {key: local_polynomial(value, point, 6, 3) for key, value in table.items()}
        coefficients = flatten_polynomials(local, dimensions, cut)
        effective_orders = [value - int(shift_twice) for value in record['ramified_orders']]
        channel_sums = Counter(sum(word) for word in product(effective_orders, repeat=3))
        ranks = []
        for j in range(7):
            actual_rank = rank(coefficients[j], 3)
            channels = sum(count for total, count in channel_sums.items() if total <= 2 * j)
            cleared_index = 2 * j + 3 * (int(shift_twice) + record['clearing'])
            assert cleared_index >= 0
            bound = min(36, (cleared_index + 1) * channels)
            assert actual_rank <= bound
            ranks.append({'jet': j, 'rank': actual_rank, 'bound': bound})
        controls.append({'point': point, 'cut': cut, 'complete_type_entries': 6 ** 4,
                         'scalar_height_checks': 6 ** 4, 'ranks': ranks})
    return controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    basis = [[1, 0, 0], [2, 1, 0], [1, 1, 1]]
    seed = transform_seed(basis, basis)
    cases = [(2, 0, (0, 0, Fraction(1, 2), Fraction(1, 2)),
              [2, 3, 3, 4, 4, 4, 5, 5, 6]),
             (None, 0, (Fraction(-1, 2), Fraction(-1, 2), 0, 0),
              [-2, -1, -1, 0, 0, 0, 1, 1, 2]),
             (0, 1, (-1, -1, 1, 1), [0, 2, 2, 4, 4, 4, 6, 6, 8])]
    profiles = []
    for point, cut, slopes, expected in cases:
        orders, clearing = ramified_profile(seed, point, cut, slopes)
        assert orders == expected, (point, orders, expected)
        profiles.append({'point': point, 'cut': cut,
                         'slopes': [[Fraction(s).numerator, Fraction(s).denominator] for s in slopes],
                         'ramified_orders': orders, 'clearing': clearing,
                         'mean_shift': str(sum(slopes))})
    scale = Fraction(23, 10)
    partition = (1 + 1 / scale + 1 / scale ** 2) ** 2
    pressure = scale ** 2 * partition ** 2
    assert pressure == Fraction(859 ** 4, 100 * 23 ** 6)
    numerator, denominator = pressure.numerator, pressure.denominator
    assert numerator ** 250 < 3 ** 821 * denominator ** 250
    # This witness is strictly stronger than the earlier rho=9/10 one.
    assert Fraction(821, 1000) < Fraction(9, 10)
    typed = typed_controls(seed, profiles)
    seed_json = json.dumps([[list(key), value] for key, value in sorted(seed.items())], separators=(',', ':'))
    report = {'status': 'pass', 'prime': 3, 'alphabet': 3, 'degree': 2,
              'scope': 'Exact uniform mixed-place construction below rho=821/1000; factoring score unchanged.',
              'basis_a': basis, 'basis_b': basis, 'inverse_basis': inverse(basis),
              'seed_sha256': sha256(seed_json.encode()).hexdigest(),
              'nonzero_seed_entries': len(seed), 'endpoint_entries': 81,
              'paired_types': [[1, 3]] * 3, 'ramification': 2,
              'profiles': profiles, 'typed_controls': typed,
              'rational_scale': [23, 10], 'rho': [821, 1000],
              'pressure': [numerator, denominator],
              'integer_certificate': {'left_base': numerator, 'right_base': denominator,
                                      'power': 250, 'alphabet_power': 821,
                                      'left_bits': (numerator ** 250).bit_length(),
                                      'right_bits': (3 ** 821 * denominator ** 250).bit_length()},
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'rho': [821, 1000], 'exact_profiles': len(profiles),
                      'typed_rank_checks': sum(len(record['ranks']) for record in typed),
                      'complete_type_entries': 6 ** 4, 'pressure': str(pressure)}, sort_keys=True))


if __name__ == '__main__':
    main()
