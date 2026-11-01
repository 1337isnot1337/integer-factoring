#!/usr/bin/env python3
"""Independent finite checks for the intersection matching crossing."""

import argparse
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import random

from check_permutation_homogeneous_labels import rank_two
from fermionic_piece_polynomial import polynomial_product
from intersection_matching_lift import lifted_polynomial, matching_labels


def multiply(a, b):
    """F_16 with modulus X^4+X+1."""
    answer = 0
    while b:
        if b & 1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19
    return answer


def power(a, exponent):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = multiply(answer, a)
        a = multiply(a, a)
        exponent >>= 1
    return answer


def pfaffian_eliminate(matrix):
    """Alternating Gaussian elimination over F_16."""
    matrix = [row[:] for row in matrix]
    size = len(matrix)
    if size % 2:
        return 0
    result = 1
    for first in range(0, size, 2):
        partner = next((j for j in range(first + 1, size)
                        if matrix[first][j]), None)
        if partner is None:
            return 0
        second = first + 1
        if partner != second:
            matrix[second], matrix[partner] = matrix[partner], matrix[second]
            for row in matrix:
                row[second], row[partner] = row[partner], row[second]
        pivot = matrix[first][second]
        result = multiply(result, pivot)
        inv = power(pivot, 14)
        for i in range(second + 1, size):
            for j in range(i + 1, size):
                update = (
                    multiply(matrix[first][i], matrix[second][j])
                    ^ multiply(matrix[first][j], matrix[second][i])
                )
                matrix[i][j] ^= multiply(inv, update)
                matrix[j][i] = matrix[i][j]
    return result


def pfaffian_match(matrix):
    @lru_cache(None)
    def recurse(remaining):
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        result = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            second = bit.bit_length() - 1
            result ^= multiply(matrix[first][second], recurse(rest ^ bit))
        return result
    return recurse((1 << len(matrix)) - 1)


def point_matrix(labels, vertex_count, value_l, values_z, kernel, scalar=1):
    a, b, u, y = labels
    particles = (
        [("U", e) for e in u] + [("B", e) for e in b]
        + [("A", e) for e in a] + [("Y", e) for e in y]
    )

    def z(edge):
        value = value_l
        for j in range(vertex_count):
            if edge >> j & 1:
                value = multiply(value, values_z[j])
        return value

    def t(edge):
        value = 1 ^ value_l
        for j in range(vertex_count):
            if edge >> j & 1:
                value = multiply(value, 1 ^ values_z[j])
        return value

    def k(e, f):
        if kernel == "disjointness":
            return int(not e & f)
        if kernel == "cross_copy":
            return 1
        return int(e == f) ^ multiply(scalar, (e & f).bit_count() % 2)

    count = len(particles)
    matrix = [[0] * count for _ in range(count)]
    for i, (left, e) in enumerate(particles):
        for j in range(i + 1, count):
            right, f = particles[j]
            pair = left, right
            value = 0
            if kernel != "cross_copy":
                if pair == ("U", "B") and e == f:
                    value = t(f)
                elif pair == ("U", "A") and e == f:
                    value = multiply(z(e), t(f))
                elif pair == ("U", "Y") and e == f:
                    value = z(f)
            if pair == ("B", "A"):
                value = multiply(k(e, f), z(f))
            elif pair == ("B", "Y"):
                value = multiply(k(e, f), multiply(z(f), t(e)))
            elif pair == ("A", "Y"):
                value = multiply(k(e, f), t(e))
            matrix[i][j] = matrix[j][i] = value
    return matrix


def evaluate_packed(polynomial, vertices, value_l, values_z):
    radix = 2 * vertices + 1
    answer = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        code = bit.bit_length() - 1
        value = power(value_l, code % radix)
        code //= radix
        for variable in values_z:
            value = multiply(value, power(variable, code % 5))
            code //= 5
        assert code == 0
        answer ^= value
    return answer


def determinant_field(matrix):
    matrix = [row[:] for row in matrix]
    size = len(matrix)
    result = 1
    for column in range(size):
        pivot = next((r for r in range(column, size) if matrix[r][column]), None)
        if pivot is None:
            return 0
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        result = multiply(result, matrix[column][column])
        inv = power(matrix[column][column], 14)
        for r in range(column + 1, size):
            factor = multiply(matrix[r][column], inv)
            for c in range(column, size):
                matrix[r][c] ^= multiply(factor, matrix[column][c])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    rng = random.Random(831742)
    gram_count = 0
    parameter_gram_count = 0
    for vertices in range(2, 9, 2):
        labels = matching_labels(vertices)
        b = vertices // 2
        for a, v in product(labels, repeat=2):
            rows = [sum(int(bool(e & f)) << j for j, f in enumerate(v)) for e in a]
            assert int(rank_two(rows) == b) == int(a == v)
            gram_count += 1
        for _ in range(100):
            a, v = rng.choice(labels), rng.choice(labels)
            for scalar in range(16):
                matrix = [
                    [int(e == f) ^ multiply(scalar, (e & f).bit_count() % 2)
                     for f in v] for e in a
                ]
                assert determinant_field(matrix) == int(a == v)
                parameter_gram_count += 1

    incidence_cases = 0
    for vertices in range(2, 25, 2):
        edges = tuple((1 << i) | (1 << j) for i, j in combinations(range(vertices), 2))
        q = [sum(((e & f).bit_count() % 2) << j for j, f in enumerate(edges))
             for e in edges]
        assert rank_two(q) == vertices - 2
        for row in q:
            square = 0
            while row:
                bit = row & -row
                row ^= bit
                square ^= q[bit.bit_length() - 1]
            assert square == 0
        incidence_cases += 1

    generic_folds = 0
    for size in (2, 4, 6, 8):
        for _ in range(50):
            m = [[0] * size for _ in range(size)]
            n = [[0] * size for _ in range(size)]
            for i in range(size):
                for j in range(i + 1, size):
                    m[i][j] = m[j][i] = rng.randrange(16)
                    n[i][j] = n[j][i] = rng.randrange(16)
            full = [m[i] + n[i] for i in range(size)] + [n[i] + m[i] for i in range(size)]
            folded = [[m[i][j] ^ n[i][j] for j in range(size)] for i in range(size)]
            assert pfaffian_eliminate(full) == power(pfaffian_eliminate(folded), 2)
            assert pfaffian_eliminate(folded) == pfaffian_match(folded)
            generic_folds += 1

    nonalternating = [[0, 0, 1, 0], [0, 0, 0, 1],
                      [1, 0, 0, 0], [0, 1, 0, 0]]
    assert pfaffian_eliminate(nonalternating) == 1
    assert pfaffian_eliminate([[0, 0], [0, 0]]) == 0

    point_checks = 0
    structured_folds = 0
    coefficient_crossings = 0
    for vertices in (2, 4, 6, 8):
        labels = matching_labels(vertices)
        for trial in range(32):
            external = tuple(rng.choice(labels) for _ in range(4))
            value_l = rng.randrange(16)
            values_z = [rng.randrange(16) for _ in range(vertices)]
            disjoint = point_matrix(external, vertices, value_l, values_z, "disjointness")
            cross = point_matrix(external, vertices, value_l, values_z, "cross_copy")
            overlap = point_matrix(external, vertices, value_l, values_z, "intersection")
            size = len(disjoint)
            assert overlap == [[disjoint[i][j] ^ cross[i][j]
                                for j in range(size)] for i in range(size)]
            full = [disjoint[i] + cross[i] for i in range(size)] + [
                cross[i] + disjoint[i] for i in range(size)]
            assert pfaffian_eliminate(full) == power(pfaffian_eliminate(overlap), 2)
            structured_folds += 1
            if vertices <= 6:
                assert pfaffian_eliminate(overlap) == pfaffian_match(overlap)
                polynomial = lifted_polynomial(*external, vertices)
                assert evaluate_packed(polynomial, vertices, value_l, values_z) == pfaffian_eliminate(overlap)
                expected = int(external[0] == external[1] and external[2] == external[3])
                assert polynomial.bit_count() % 2 == expected
                coefficient_crossings += 1
                polynomial_product.cache_clear()
            for scalar in (0, 1, 2, 7, 15):
                at_crossing = point_matrix(external, vertices, 1, [1] * vertices,
                                           "intersection", scalar)
                assert pfaffian_eliminate(at_crossing) == int(
                    external[0] == external[1] and external[2] == external[3])
                point_checks += 1

    result = {
        "full_matching_gram_pairs": gram_count,
        "scalar_parameter_gram_pairs_over_F16": parameter_gram_count,
        "complete_edge_incidence_rank_and_square_cases": incidence_cases,
        "generic_alternating_copy_folds": generic_folds,
        "structured_copy_folds": structured_folds,
        "independent_polynomial_and_pfaffian_checks": coefficient_crossings,
        "parameter_crossing_evaluations": point_checks,
        "negative_control": "nonalternating cross-copy block breaks the fold",
        "scope": "exact crossing constructor and scalar deformation; no improved joint-width bound",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
