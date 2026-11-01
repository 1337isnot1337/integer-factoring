#!/usr/bin/env python3
"""Exact controls for the non-Gaussian binary cube profile floor."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random

import cube_homogeneous_gaussian as cube


def convolution(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        if a:
            for j, b in enumerate(second):
                if b:
                    result[i + j] += a * b
    return result


def power(polynomial, exponent):
    result = [1]
    for _ in range(exponent):
        result = convolution(result, polynomial)
    return result


def positive_certificate():
    base = [0] * 73
    for degree, coefficient in ((0, 1), (6, 1), (12, 1), (18, 1), (36, 7), (54, 4), (72, 1)):
        base[degree] = coefficient
    polynomial = power(base, 5)
    for j in range(13):
        polynomial[10 * j] -= comb(12, j)
    assert not any(polynomial[:6]) and polynomial[6] == 5
    quotient = polynomial[6:]
    degree = len(quotient) - 1
    transformed = [sum(value * comb(degree - i, k - i)
                       for i, value in enumerate(quotient[:k + 1]))
                   for k in range(degree + 1)]
    assert degree == 354 and min(transformed) == 5
    assert all(value > 0 for value in transformed)
    assert 3 * 3 ** 2 > 5 ** 2
    return {'leading_degree': 6, 'quotient_degree': degree,
            'polynomial': polynomial, 'positive_coefficients': transformed,
            'minimum_positive_coefficient': min(transformed),
            'coefficient_count': len(transformed)}


def endpoint(mask):
    return int(all((mask >> i & 1) == (mask >> (i ^ 7) & 1)
                   for i in range(4)))


def random_tensor(seed, height):
    rng = random.Random(seed)
    result = []
    for mask in range(256):
        degree = height * mask.bit_count()
        if not degree:
            result.append(cube.ONE)
            continue
        terms = set()
        for _ in range(3):
            a = rng.randrange(degree + 1)
            b = rng.randrange(degree - a + 1)
            term = (a, b, degree - a - b)
            terms.symmetric_difference_update((term,))
        if len(terms) % 2 != endpoint(mask):
            terms.symmetric_difference_update(((degree, 0, 0),))
        result.append(frozenset(terms))
    return tuple(result)


def selected_minor(tensor, representatives, allowed=None):
    opposite = tuple(i ^ 7 for i in representatives)
    def occupation(vertices, bits):
        return sum((1 << vertex) for j, vertex in enumerate(vertices) if bits >> j & 1)
    allowed = tuple(range(1 << len(representatives))) if allowed is None else tuple(allowed)
    count = len(allowed)
    matrix = [[tensor[occupation(representatives, row) | occupation(opposite, column)]
               for column in allowed] for row in allowed]
    determinant = cube.minor(matrix, range(count), range(count))
    return matrix, determinant


def check_tensor(tensor, height, counts):
    for mask, polynomial in enumerate(tensor):
        assert all(sum(term) == height * mask.bit_count() for term in polynomial)
        assert len(polynomial) % 2 == endpoint(mask)
        counts['homogeneous_endpoint_entries'] += 1
    assert any(tensor[mask] for mask in range(256) if mask.bit_count() % 2)
    for axis in range(3):
        matrix = cube.axis_matrix(tensor, axis)
        first, second = cube.axis_parts(axis)
        permutation = [second.index(vertex ^ 7) for vertex in first]
        for row in range(16):
            paired = sum((row >> i & 1) << permutation[i] for i in range(4))
            for column in range(16):
                assert len(matrix[row][column]) % 2 == int(column == paired)
                counts['complete_cut_endpoint_entries'] += 1
        row_occupation = sum(word.bit_count() for word in range(16))
        assert row_occupation == 32
        assert 2 * height * row_occupation == 64 * height
        counts['full_determinant_degree_identities'] += 1
        occupations = sorted(word.bit_count() for word in range(16))
        for size in range(1, 17):
            degree = 2 * height * sum(occupations[:size])
            assert degree >= 0 and degree <= 64 * height
            counts['all_prefix_minor_degree_identities'] += 1
    matrix, value = selected_minor(tensor, (0,))
    assert value and len(value) % 2 == 1
    assert all(sum(term) == 2 * height for term in value)
    orders = [cube.valuation(value, axis) for axis in range(3)]
    assert sum(orders) <= 2 * height
    for axis in range(3):
        full = cube.axis_matrix(tensor, axis)
        first, second = cube.axis_parts(axis)
        rows = (0, 1 << first.index(0))
        columns = (0, 1 << second.index(7))
        assert matrix == [[full[i][j] for j in columns] for i in rows]
        counts['shared_two_by_two_selections'] += 1
    for a, b in combinations(range(3), 2):
        omitted = 3 - a - b
        representatives = (0, 1 << omitted)
        selected, value = selected_minor(tensor, representatives)
        assert value and len(value) % 2 == 1
        assert all(sum(term) == 8 * height for term in value)
        assert cube.valuation(value, a) + cube.valuation(value, b) <= 8 * height
        for axis in (a, b):
            full = cube.axis_matrix(tensor, axis)
            first, second = cube.axis_parts(axis)
            rows = [sum((bits >> j & 1) << first.index(vertex)
                        for j, vertex in enumerate(representatives)) for bits in range(4)]
            columns = [sum((bits >> j & 1) << second.index(vertex ^ 7)
                           for j, vertex in enumerate(representatives)) for bits in range(4)]
            assert selected == [[full[i][j] for j in columns] for i in rows]
            counts['shared_four_by_four_selections'] += 1
        counts['exact_four_by_four_determinants'] += 1
        smaller, determinant = selected_minor(tensor, representatives, (0, 1, 2))
        assert determinant and len(determinant) % 2 == 1
        assert all(sum(term) == 4 * height for term in determinant)
        assert cube.valuation(determinant, a) + cube.valuation(determinant, b) <= 4 * height
        assert smaller == [row[:3] for row in selected[:3]]
        counts['exact_shared_three_by_three_determinants'] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    counts = Counter()
    for height in (1, 2, 3, 5):
        for seed in range(4):
            check_tensor(random_tensor(100 * height + seed, height), height, counts)
    sigma = [0, 2, 4, 6] + [12] * 7 + [18] * 4 + [24]
    assert sum(sigma[:2]) == 2 and sum(sigma[:4]) == 12 and sum(sigma) == 192
    assert sum(sigma[:3]) == 6
    occupations = sorted(word.bit_count() for word in range(16))
    assert all(sum(sigma[:k]) == 6 * sum(occupations[:k]) for k in range(5, 17))
    counts['seeded_non_gaussian_tensors'] = 16
    result = {'theorem': 'Every homogeneous binary cube tensor has common-scale '
              'unweighted Smith-profile pressure strictly greater than 3/5.',
              'scope': 'A common binary physical type; no lower bound on actual coefficient ranks.',
              'optimistic_orders_at_height_three': sigma,
              'counts': dict(counts), 'positive_certificate': positive_certificate(),
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', **dict(counts),
                      'positive_coefficients': 355, 'minimum_positive_coefficient': 5}, sort_keys=True))


if __name__ == '__main__':
    main()
