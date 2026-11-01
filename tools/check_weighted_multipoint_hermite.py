#!/usr/bin/env python3
"""Exact type restrictions for weighted same-place Hermite estimates.

Finite controls for the separate uniform transfer and rho > 3/4 floor.
All ranks and pressure comparisons use exact field or rational arithmetic.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import (
    affine_seed, local_polynomial, mixing_seed, multiply, rank,
)


def smith_orders(coefficients, prime):
    """Smith orders from ranks of complete lower block-Toeplitz matrices."""
    size = len(coefficients[0])
    assert size and all(len(row) == size for block in coefficients for row in block)
    degree = len(coefficients) - 1
    previous_rank = previous_increment = 0
    orders = []
    for length in range(1, size * degree + 2):
        matrix = [[0] * (size * length) for _ in range(size * length)]
        for i in range(length):
            for j in range(max(0, i - degree), i + 1):
                block = coefficients[i - j]
                for r in range(size):
                    matrix[i * size + r][j * size:(j + 1) * size] = block[r]
        current_rank = rank(matrix, prime)
        increment = current_rank - previous_rank
        assert previous_increment <= increment <= size
        orders.extend([length - 1] * (increment - previous_increment))
        previous_rank, previous_increment = current_rank, increment
        if increment == size:
            break
    assert len(orders) == size, (size, degree, orders)
    return orders


def flatten_polynomials(table, dimensions, cut):
    dx, dy, du, dv = dimensions
    assert dx == dv and dy == du
    rows = list(product(range(dx), range(dy if cut == 0 else du)))
    columns = list(product(range(du if cut == 0 else dy), range(dv)))
    degree = max((len(value) - 1 for value in table.values()), default=0)
    coefficients = []
    for k in range(degree + 1):
        block = []
        for a, b in rows:
            row = []
            for c, d in columns:
                key = (a, b, c, d) if cut == 0 else (a, c, b, d)
                value = table.get(key, ())
                row.append(value[k] if k < len(value) else 0)
            block.append(row)
        coefficients.append(block)
    return coefficients


def local_gauged_profile(seed, alphabet, degree, point, prime, cut, heights):
    local = {key: local_polynomial(value, point, degree, prime)
             for key, value in seed.items()}
    terms = {}
    for key, values in local.items():
        height = sum(heights[leg][letter] for leg, letter in enumerate(key))
        terms[key] = {k + height: value for k, value in enumerate(values) if value}
    minimum = min(k for values in terms.values() for k in values)
    clearing = max(0, -minimum)
    cleared = {}
    for key, values in terms.items():
        if not values:
            continue
        coefficients = [0] * (max(values) + clearing + 1)
        for k, value in values.items():
            coefficients[k + clearing] = value
        cleared[key] = coefficients
    orders = smith_orders(flatten_polynomials(cleared, (alphabet,) * 4, cut), prime)
    return [value - clearing for value in orders]


def type_words(counts):
    return [word for word in product(range(len(counts)), repeat=sum(counts))
            if tuple(word.count(i) for i in range(len(counts))) == tuple(counts)]


def typed_power(seed, alphabet, copies, prime, xcounts, ycounts):
    assert sum(xcounts) == sum(ycounts) == copies
    xwords, ywords = type_words(xcounts), type_words(ycounts)
    word_sets = (xwords, ywords, ywords, xwords)
    dimensions = tuple(map(len, word_sets))
    result = {}
    for key in product(*(range(d) for d in dimensions)):
        words = [word_sets[leg][index] for leg, index in enumerate(key)]
        value = [1]
        for position in range(copies):
            factor = seed.get(tuple(word[position] for word in words), ())
            if not factor:
                value = []
                break
            value = multiply(value, factor, prime)
        if any(value):
            result[key] = value
    return result, dimensions


def partition_sum(orders, scale):
    return sum((scale ** (-order) for order in orders), Fraction(0))


def check_case(seed, alphabet, degree, prime, copies, xcounts, ycounts):
    table, dimensions = typed_power(seed, alphabet, copies, prime, xcounts, ycounts)
    for key in product(*(range(d) for d in dimensions)):
        x, y, u, v = key
        assert sum(table.get(key, ())) % prime == int(x == v and y == u)
    points = [value for value in range(prime) if value != 1] + [None]
    profiles_checked = comparisons = negative_orders = 0
    records = []
    scales = (Fraction(3, 2), Fraction(2), Fraction(5, 2))
    gauges = [tuple(tuple(0 for _ in range(alphabet)) for _ in range(4)),
              tuple(tuple(((2 * leg + 3 * a + 1) % 5) - 2
                          for a in range(alphabet)) for leg in range(4)),
              tuple(tuple(((-1) ** (leg + a)) * (a + 1)
                          for a in range(alphabet)) for leg in range(4))]
    for point in points:
        local = {key: local_polynomial(value, point, degree * copies, prime)
                 for key, value in table.items()}
        actual = [smith_orders(flatten_polynomials(local, dimensions, cut), prime)
                  for cut in (0, 1)]
        for gauge in gauges:
            shift = (sum(xcounts[a] * (gauge[0][a] + gauge[3][a])
                         for a in range(alphabet))
                     + sum(ycounts[a] * (gauge[1][a] + gauge[2][a])
                           for a in range(alphabet)))
            for cut in (0, 1):
                seed_orders = local_gauged_profile(
                    seed, alphabet, degree, point, prime, cut, gauge)
                negative_orders += sum(order < 0 for order in seed_orders)
                profiles_checked += 1
                for scale in scales:
                    left = partition_sum(actual[cut], scale)
                    right = scale ** shift * partition_sum(seed_orders, scale) ** copies
                    assert left <= right, (point, cut, gauge, scale, left, right)
                    comparisons += 1
        records.append({'point': point, 'typed_profiles': actual})
    return {'prime': prime, 'alphabet': alphabet, 'copies': copies,
            'xcounts': xcounts, 'ycounts': ycounts, 'dimensions': dimensions,
            'complete_endpoint_entries': dimensions[0] ** 2 * dimensions[1] ** 2,
            'local_gauged_profiles': profiles_checked, 'rational_comparisons': comparisons,
            'negative_order_channels': negative_orders, 'profiles': records}


def negative_control():
    seed = affine_seed(2)
    counts = (1, 1)
    table, dimensions = typed_power(seed, 2, 2, 2, counts, counts)
    actual = smith_orders(flatten_polynomials(table, dimensions, 0), 2)
    heights = ((10, 10), (0, 0), (0, 0), (0, 0))
    local = local_gauged_profile(seed, 2, 1, 0, 2, 0, heights)
    scale = Fraction(2)
    actual_sum = partition_sum(actual, scale)
    unshifted = partition_sum(local, scale) ** 2
    assert actual_sum > unshifted
    assert actual_sum <= scale ** 20 * unshifted
    return {'actual_sum': str(actual_sum), 'incorrect_unshifted_bound': str(unshifted),
            'corrected_bound': str(scale ** 20 * unshifted), 'missing_shift_rejected': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [check_case(affine_seed(2), 2, 1, 2, 2, (1, 1), (1, 1)),
             check_case(affine_seed(3), 2, 1, 3, 3, (1, 2), (2, 1)),
             check_case(affine_seed(5), 2, 1, 5, 2, (2, 0), (1, 1)),
             check_case(mixing_seed(), 3, 2, 3, 2, (1, 1, 0), (0, 1, 1))]
    control = negative_control()
    result = {'schema_version': 1,
              'claim': 'Finite controls for weighted same-place Hermite estimates; no score improvement.',
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'cases': cases, 'negative_control': control,
              'totals': {key: sum(case[key] for case in cases) for key in
                         ('complete_endpoint_entries', 'local_gauged_profiles',
                          'rational_comparisons', 'negative_order_channels')}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals'], sort_keys=True))


if __name__ == '__main__':
    main()
