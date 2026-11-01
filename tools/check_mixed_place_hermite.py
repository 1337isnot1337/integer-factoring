#!/usr/bin/env python3
"""Exact opposite-flag crossing reconstruction at several rational places.

This proves finite identities and a rational rho < 9/10 seed certificate.
It does not improve the established factoring exponent one half.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path

from check_multipoint_hermite_crossings import (
    affine_seed, flatten, local_polynomial, local_profiles,
    low_order_word_count, mixing_seed, multiply, rank, solve, tensor_power,
)


def nullspace(matrix, prime):
    """A deterministic basis of the right kernel over a prime field."""
    work = [[value % prime for value in row] for row in matrix]
    columns = len(work[0])
    pivot_columns = []
    r = 0
    for column in range(columns):
        pivot = next((i for i in range(r, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        inverse = pow(work[r][column], -1, prime)
        work[r] = [value * inverse % prime for value in work[r]]
        for i in range(len(work)):
            if i == r:
                continue
            factor = work[i][column]
            if factor:
                work[i] = [(a - factor * b) % prime
                           for a, b in zip(work[i], work[r])]
        pivot_columns.append(column)
        r += 1
        if r == len(work):
            break
    free_columns = [column for column in range(columns) if column not in pivot_columns]
    result = []
    for column in free_columns:
        vector = [0] * columns
        vector[column] = 1
        for i, pivot in enumerate(pivot_columns):
            vector[pivot] = -work[i][column] % prime
        assert all(sum(a * b for a, b in zip(row, vector)) % prime == 0
                   for row in matrix)
        result.append(vector)
    return result


def jet_vector(degree, point, order, prime):
    if point is None:
        return [int(k == degree - order) for k in range(degree + 1)]
    return [(comb(k, order) * pow(point, k - order, prime) % prime if k >= order else 0)
            for k in range(degree + 1)]


def jet_flag(degree, pattern, prime):
    assert pattern
    counts = {point: 0 for point in pattern}
    result = []
    for index in range(degree + 1):
        point = pattern[index % len(pattern)]
        order = counts[point]
        counts[point] += 1
        result.append({'point': point, 'order': order,
                       'vector': jet_vector(degree, point, order, prime)})
    return result


def combine(coefficients, rows, prime):
    return [sum(coefficient * row[k] for coefficient, row in zip(coefficients, rows)) % prime
            for k in range(len(rows[0]))]


def opposite_basis(degree, row_pattern, column_pattern, prime):
    assert set(row_pattern).isdisjoint(column_pattern)
    assert 1 not in row_pattern and 1 not in column_pattern
    row_flag = jet_flag(degree, row_pattern, prime)
    column_flag = jet_flag(degree, column_pattern, prime)
    row_vectors = [jet['vector'] for jet in row_flag]
    column_vectors = [jet['vector'] for jet in column_flag]
    assert rank(row_vectors, prime) == rank(column_vectors, prime) == degree + 1
    dimensions_checked = 0
    for r in range(degree + 2):
        for c in range(degree + 2):
            intersection_dimension = r + c - rank(row_vectors[:r] + column_vectors[:c], prime)
            assert intersection_dimension == max(0, r + c - degree - 1)
            dimensions_checked += 1
    basis = []
    for k in range(degree + 1):
        rows = row_vectors[:k + 1]
        columns = column_vectors[:degree - k + 1]
        matrix = [[row[j] for row in rows] + [-column[j] % prime for column in columns]
                  for j in range(degree + 1)]
        kernel = nullspace(matrix, prime)
        assert len(kernel) == 1
        coefficients = kernel[0]
        left, right = coefficients[:k + 1], coefficients[k + 1:]
        vector = combine(left, rows, prime)
        assert vector == combine(right, columns, prime)
        pivot = next(value for value in vector if value)
        inverse = pow(pivot, -1, prime)
        vector = [value * inverse % prime for value in vector]
        left = [value * inverse % prime for value in left]
        right = [value * inverse % prime for value in right]
        assert left[-1] != 0
        basis.append({'vector': vector, 'row_coefficients': left,
                      'column_coefficients': right})
    vectors = [entry['vector'] for entry in basis]
    assert rank(vectors, prime) == degree + 1
    weights = solve(list(map(list, zip(*vectors))), [1] * (degree + 1), prime)
    assert combine(weights, vectors, prime) == [1] * (degree + 1)
    return row_flag, column_flag, basis, weights, dimensions_checked


def apply_functional(vector, polynomial, prime):
    return sum(a * b for a, b in zip(vector, polynomial)) % prime


def check_power(seed, alphabet, degree, copies, prime, row_pattern, column_pattern):
    full_degree = degree * copies
    row_flag, column_flag, basis, weights, dimension_checks = opposite_basis(
        full_degree, row_pattern, column_pattern, prime)
    table = tensor_power(seed, alphabet, copies, prime)
    dimension = alphabet ** copies
    all_places = list(dict.fromkeys(row_pattern + column_pattern))
    profiles = {point: local_profiles(seed, alphabet, degree, point, prime)
                for point in all_places}
    locals_by_place = {
        point: {key: local_polynomial(value, point, full_degree, prime)
                for key, value in table.items()} for point in all_places}
    pieces = []
    piece_ranks = []
    representation_entries = 0
    for k, (entry, weight) in enumerate(zip(basis, weights)):
        piece = {}
        for key, value in table.items():
            expected = weight * apply_functional(entry['vector'], value, prime) % prime
            left = sum(coefficient * locals_by_place[jet['point']][key][jet['order']]
                       for coefficient, jet in zip(entry['row_coefficients'], row_flag)) % prime
            right = sum(coefficient * locals_by_place[jet['point']][key][jet['order']]
                        for coefficient, jet in zip(entry['column_coefficients'], column_flag)) % prime
            assert expected == weight * left % prime == weight * right % prime
            representation_entries += 2
            if expected:
                piece[key] = [expected]
        pieces.append(piece)
        ranks, bounds = [], []
        for cut, flag, coefficients in ((0, row_flag, entry['row_coefficients']),
                                        (1, column_flag, entry['column_coefficients'])):
            actual = rank(flatten(piece, dimension, cut), prime)
            bound = sum((jet['order'] + 1) * low_order_word_count(
                profiles[jet['point']][cut], copies, jet['order'])
                for coefficient, jet in zip(coefficients, flag) if coefficient)
            if not weight:
                bound = 0
            assert actual <= bound
            ranks.append(actual)
            bounds.append(bound)
        piece_ranks.append({'index': k, 'weight': weight, 'ranks': ranks,
                            'profile_rank_bounds': bounds})
    for key in product(range(dimension), repeat=4):
        x, y, u, v = key
        expected = int(x == v and y == u)
        assert sum(table.get(key, ())) % prime == expected
        assert sum(piece.get(key, [0])[0] for piece in pieces) % prime == expected
    return {'prime': prime, 'alphabet': alphabet, 'seed_degree_bound': degree,
            'copies': copies, 'row_pattern': row_pattern, 'column_pattern': column_pattern,
            'complete_endpoint_entries': dimension ** 4,
            'representation_entries': representation_entries,
            'flag_dimension_checks': dimension_checks,
            'exact_rank_checks': 2 * len(piece_ranks), 'pieces': piece_ranks,
            'basis': basis, 'weights': weights}


def ternary_certificate():
    seed = mixing_seed()
    expected = {0: [[0] * 9, [0, 1, 1, 1, 1, 2, 4, 4, 4]],
                2: [[0, 0, 0, 1, 1, 1, 2, 2, 2], [0] * 9],
                None: [[0, 0, 0, 1, 1, 1, 2, 2, 2], [0] * 9]}
    actual = {point: local_profiles(seed, 3, 2, point, 3) for point in expected}
    assert actual == expected
    # These local determinant orders consume the whole homogeneous degree.
    assert sum(actual[2][0]) + sum(actual[None][0]) == 18
    assert sum(actual[0][1]) == 18
    scale = Fraction(5, 3)
    row = 3 * (1 + scale ** -2 + scale ** -4)
    column = 1 + 4 * scale ** -1 + scale ** -2 + 3 * scale ** -4
    bound = scale ** 2 * row * column
    assert row == Fraction(2793, 625)
    assert column == Fraction(2593, 625)
    assert bound == Fraction(2414083, 46875)
    assert bound ** 5 < 3 ** 18
    assert 54 ** 5 > 3 ** 18
    assert sum([0, 1, 1, 1, 1, 2]) == 6
    assert sum([0, 0, 0, 1, 1, 1]) == 3
    return {'scale': str(scale), 'row_bound': str(row), 'column_bound': str(column),
            'complete_bound': str(bound), 'rho': '9/10',
            'strict_power_inequality': {
                'left': str(bound.numerator ** 5),
                'right': str(3 ** 18 * bound.denominator ** 5)},
            'old_same_place_lower_bound_base': 54,
            'old_criterion_at_9_over_10_rejected': True,
            'profiles': [{'point': point, 'cuts': actual[point]} for point in expected]}


def overlap_negative_control():
    degree, prime = 3, 3
    rows = [entry['vector'] for entry in jet_flag(degree, [0], prime)]
    columns = [entry['vector'] for entry in jet_flag(degree, [0], prime)]
    r, c = 2, 3
    actual = r + c - rank(rows[:r] + columns[:c], prime)
    claimed = max(0, r + c - degree - 1)
    assert actual == 2 and claimed == 1
    return {'overlapping_point': 0, 'actual_intersection_dimension': actual,
            'disjoint_formula_dimension': claimed, 'rejected': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    seed_with_content = {key: multiply(value, [2, 4], 5)
                         for key, value in affine_seed(5).items()}
    cases = [check_power(affine_seed(2), 2, 1, 3, 2, [None], [0]),
             check_power(mixing_seed(), 3, 2, 1, 3, [2, None], [0]),
             check_power(mixing_seed(), 3, 2, 2, 3, [2, None], [0]),
             check_power(seed_with_content, 2, 2, 2, 5, [0, 2], [None, 3]),
             check_power(affine_seed(5), 2, 3, 2, 5, [0, 2], [None, 3])]
    certificate = ternary_certificate()
    negative = overlap_negative_control()
    result = {'schema_version': 1,
              'claim': 'Mixed-place opposite-flag reconstruction and a strict rho < 9/10 ternary seed; factoring score unchanged.',
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'cases': cases, 'ternary_certificate': certificate,
              'overlap_negative_control': negative,
              'totals': {key: sum(case[key] for case in cases) for key in
                         ('complete_endpoint_entries', 'representation_entries',
                          'flag_dimension_checks', 'exact_rank_checks')}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals'], sort_keys=True))


if __name__ == '__main__':
    main()
