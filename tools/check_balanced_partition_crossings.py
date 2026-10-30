#!/usr/bin/env python3
"""Independent finite checks for balanced-partition crossing constructors."""

import argparse
from itertools import combinations, permutations, product
import json
from math import factorial
from pathlib import Path
import random

from check_permutation_homogeneous_labels import rank_two
from even_partition_crossing import (
    coefficient, coefficient_terms, kernel, lifted_polynomial,
    partition_labels,
)
from fermionic_piece_polynomial import polynomial_product
from intersection_matching_lift import lifted_polynomial as matching_polynomial


F4 = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))


def determinant(matrix, multiply=lambda x, y: x & y):
    """Leibniz expansion, valid even when coefficients are nonunits."""
    result = 0
    for order in permutations(range(len(matrix))):
        value = 1
        for i, j in enumerate(order):
            value = multiply(value, matrix[i][j])
        result ^= value
    return result


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


def product_rows(left, right):
    result = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            row ^= bit
            value ^= right[bit.bit_length() - 1]
        result.append(value)
    return result


def expand_terms(terms, vertices, particles):
    radix = 4 * particles + 1
    places = tuple(radix * 5**i for i in range(vertices))
    result = 0
    for z, t, z_counts, t_counts in terms:
        polynomial = 1 << (z + sum(p * c for p, c in zip(places, z_counts)))
        for _ in range(t):
            polynomial ^= polynomial << 1
        for place, exponent in zip(places, t_counts):
            for _ in range(exponent):
                polynomial ^= polynomial << place
        result ^= polynomial
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    rng = random.Random(381772)
    gram_checks = scalar_checks = incidence_checks = 0
    counts = []
    for vertices, block_size in ((4, 2), (6, 2), (8, 2), (8, 4), (12, 4), (12, 6)):
        labels = partition_labels(vertices, block_size)
        q = vertices // block_size
        expected = factorial(vertices) // (factorial(block_size)**q * factorial(q))
        assert len(labels) == expected
        cases = product(labels, repeat=2) if len(labels) <= 35 else (
            (rng.choice(labels), rng.choice(labels)) for _ in range(1000))
        for a, b in cases:
            matrix = [[kernel(e, f) for f in b] for e in a]
            assert determinant(matrix) == int(a == b)
            gram_checks += 1
        for _ in range(48):
            a = rng.choice(labels)
            b = a if rng.randrange(3) == 0 else rng.choice(labels)
            s = rng.randrange(4)
            matrix = [[int(e == f) ^ F4[s][(e & f).bit_count() & 1]
                       for f in b] for e in a]
            assert determinant(matrix, lambda x, y: F4[x][y]) == int(a == b)
            scalar_checks += 1
        modes = tuple(sum(1 << i for i in selected)
                      for selected in combinations(range(vertices), block_size))
        incidence_gram = [sum(((e & f).bit_count() & 1) << j
                              for j, f in enumerate(modes)) for e in modes]
        assert rank_two(incidence_gram) == vertices - 2
        assert not any(product_rows(incidence_gram, incidence_gram))
        incidence_checks += 1
        counts.append([vertices, block_size, len(labels)])

    polynomial_checks = matching_checks = coefficient_checks = 0
    for vertices, block_size in ((4, 2), (6, 2), (8, 4)):
        labels = partition_labels(vertices, block_size)
        q = vertices // block_size
        cases = []
        for _ in range(24):
            a, b, u, y = (rng.choice(labels) for _ in range(4))
            cases.append((a, b, u, y))
            cases.append((a, a, u, u))
        for a, b, u, y in cases:
            terms = coefficient_terms(a, b, u, y, vertices, block_size)
            direct = lifted_polynomial(a, b, u, y, vertices, block_size)
            assert direct == expand_terms(terms, vertices, q)
            assert direct.bit_count() % 2 == int(a == b and u == y)
            polynomial_checks += 1
            if block_size == 2:
                assert direct == matching_polynomial(a, b, u, y, vertices)
                matching_checks += 1
            for _ in range(16):
                degrees = tuple(rng.randrange(5) for _ in range(vertices))
                l_degree = rng.randrange(4 * q + 1)
                code = l_degree + (4 * q + 1) * sum(v * 5**i for i, v in enumerate(degrees))
                assert coefficient(terms, l_degree, degrees) == (direct >> code) & 1
                coefficient_checks += 1
            polynomial_product.cache_clear()

    multipartite_checks = rank_one_checks = 0
    for part_size in (2, 3, 4):
        ps = tuple(permutations(range(part_size)))
        labels = tuple(tuple(
            (1 << i) | (1 << (part_size + sigma[i])) | (1 << (2 * part_size + tau[i]))
            for i in range(part_size)) for sigma, tau in product(ps, repeat=2))
        weights = (1,) * part_size + (2,) * part_size + (3,) * part_size
        cases = product(labels, repeat=2) if len(labels) <= 36 else (
            (rng.choice(labels), rng.choice(labels)) for _ in range(1000))
        for a, b in cases:
            matrix = [[int(e == f) ^ xor_sum(weights[v] for v in range(3 * part_size)
                                             if (e & f) >> v & 1)
                       for f in b] for e in a]
            assert determinant(matrix, lambda x, y: F4[x][y]) == int(a == b)
            multipartite_checks += 1
        for _ in range(64):
            a = rng.choice(labels)
            b = a if rng.randrange(3) == 0 else rng.choice(labels)
            arbitrary_weights = [rng.randrange(4) for _ in range(3 * part_size - 1)]
            arbitrary_weights.append(xor_sum(arbitrary_weights))
            values = {f: xor_sum(arbitrary_weights[v] for v in range(3 * part_size)
                                 if f >> v & 1) for f in b}
            matrix = [[int(e == f) ^ values[f] for f in b] for e in a]
            assert determinant(matrix, lambda x, y: F4[x][y]) == int(a == b)
            rank_one_checks += 1

    # A nonbalanced block already destroys the claimed diagonal identity.
    assert determinant([[1 ^ 1]]) == 0
    result = {
        "label_counts": counts,
        "binary_gram_checks": gram_checks,
        "F4_scalar_gram_checks": scalar_checks,
        "incidence_rank_and_square_checks": incidence_checks,
        "independent_polynomial_checks": polynomial_checks,
        "matching_source_comparisons": matching_checks,
        "coefficient_extraction_checks": coefficient_checks,
        "weighted_multipartite_gram_checks": multipartite_checks,
        "rank_one_deformation_gram_checks": rank_one_checks,
        "negative_control_nonbalanced_block": "passed",
        "scope": "exact constructors and effectivity; no improved asymptotic rank profile",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
