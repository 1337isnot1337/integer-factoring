#!/usr/bin/env python3
"""Exact controls for matching-sector crossings and both complete profiles."""

from collections import Counter
from itertools import product
import json
from math import comb

from perfect_matching_crossing import (
    common_block, dense_seed, disjoint_pair_count, matching_count,
    perfect_matchings, polynomial_entry, reflected_parts,
    row_block_multiplicities, row_zero_counts, column_zero_counts,
    toeplitz_profile,
)
from neutral_reflected_crossing import (
    crossing_entry, normalize, polynomial_multiply,
)


def reverse(value, degree):
    return sum(((value >> i) & 1) << (degree - i)
               for i in range(degree + 1))


def inverse_controls(pairs):
    """Literal polynomial G(w) w^r G(w^-1)=w^r I."""
    matrix = common_block(pairs)
    reciprocal = [[reverse(value, pairs) for value in row] for row in matrix]
    n = len(matrix)
    count = 0
    for i in range(n):
        for j in range(n):
            value = 0
            for k in range(n):
                value ^= polynomial_multiply(matrix[i][k], reciprocal[k][j])
            assert value == ((1 << pairs) if i == j else 0)
            count += 1
    return count


def literal_cut(seed, column):
    q = seed['dimensions'][0]
    matrix = [[0] * (q * q) for _ in range(q * q)]
    for x, y, u, v in product(range(q), repeat=4):
        row, col = ((q * x + u, q * y + v) if column
                    else (q * x + y, q * u + v))
        matrix[row][col] = seed['values'][x + q * y + q * q * u + q ** 3 * v][0]
    return matrix


def verify():
    source_checks = endpoint_checks = inverse_checks = 0
    source_records = []
    for pairs in (3,):
        period = 1 << (2 * pairs + 1)
        labels = perfect_matchings(tuple(range(1, 2 * pairs + 1)))
        parts = tuple(reflected_parts(label, period) for label in labels)
        h = (pairs - 1) * period - pairs + 2
        w = period + pairs - 2
        assert all(sum(label) == h for label in parts)
        assert all(sum(period - k for k in label) == w for label in parts)
        for x, y, u, v in product(range(len(labels)), repeat=4):
            value = polynomial_entry(labels[x], labels[y], labels[u], labels[v])
            expected = normalize(-period * pairs, sum(
                ((value >> j) & 1) << (period * (pairs - j))
                for j in range(pairs + 1)))
            actual = crossing_entry(period, parts[x], parts[y], parts[v], parts[u])
            assert actual == expected, (x, y, u, v)
            source_checks += 1
        source_records.append({'pairs': pairs, 'period': period, 'state_weight': h,
                               'dual_weight': w, 'alphabet': len(labels)})

    # The restriction pairs>=3 in the source collapse is necessary.
    labels = perfect_matchings((1, 2, 3, 4))
    part = reflected_parts(labels[0], 32)
    assert crossing_entry(32, part, part, part, part) != (-64, 1)

    profiles = []
    toeplitz_controls = 0
    for pairs in (2, 3):
        seed = dense_seed(pairs)
        q = seed['dimensions'][0]
        endpoint_checks += q ** 4
        matrices = [literal_cut(seed, False), literal_cut(seed, True)]
        expected = [row_zero_counts(pairs), column_zero_counts(pairs)]
        for side, matrix in enumerate(matrices):
            actual, ranks = toeplitz_profile(matrix, 2 * pairs)
            assert actual == expected[side]
            toeplitz_controls += len(ranks) - 1
            reciprocal = [[reverse(value, pairs) for value in row] for row in matrix]
            reversed_profile, ranks = toeplitz_profile(reciprocal, 2 * pairs)
            predicted = ({pairs - order: count for order, count in expected[0].items()}
                         if side == 0 else {0: q * q})
            assert reversed_profile == predicted
            toeplitz_controls += len(ranks) - 1
        profiles.append({'pairs': pairs, 'alphabet': q,
                         'row_zero': expected[0],
                         'row_infinity': {-order: count for order, count in expected[0].items()},
                         'column_zero': expected[1], 'column_infinity': {-pairs: q * q}})

    incidence_checks = 0
    for pairs in range(5):
        inverse_checks += inverse_controls(pairs)
        labels = tuple(map(set, perfect_matchings(tuple(range(2 * pairs)))))
        actual = disjoint_pair_count(pairs)
        assert actual == sum(not (x & y) for x, y in product(labels, repeat=2))
        assert sum(count * matching_count(r) for r, count in
                   row_block_multiplicities(pairs).items()) == matching_count(pairs) ** 2
        if pairs >= 2:
            assert actual * (2 * pairs - 1) >= matching_count(pairs) ** 2 * (pairs - 1)
        for x in labels:
            for y in labels:
                common = len(x & y)
                recovered = 0
                for k in range(pairs + 1):
                    if comb(common, k) % 2:
                        from neutral_reflected_crossing import binomial_polynomial
                        recovered ^= binomial_polynomial(pairs - k)
                assert recovered == polynomial_entry(x, x, y, y)
                incidence_checks += 1

    return {'status': 'exact source, endpoint, reciprocal, and complete profile controls passed',
            'source_coefficient_entries': source_checks,
            'source_cases': source_records,
            'strict_endpoint_entries': endpoint_checks,
            'reciprocal_block_entries': inverse_checks,
            'partial_matching_incidence_entries': incidence_checks,
            'independent_block_toeplitz_ranks': toeplitz_controls,
            'profiles': profiles,
            'four_pair_profiles': {'row_zero': row_zero_counts(4),
                                   'column_zero': column_zero_counts(4)},
            'negative_controls': 1, 'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
