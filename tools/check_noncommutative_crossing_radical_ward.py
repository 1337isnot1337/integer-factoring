#!/usr/bin/env python3
"""Exact Ward identities and homogeneous codes from commuting symplectic forms."""

from itertools import combinations, product
from math import comb
from pathlib import Path
import argparse
import hashlib
import json
import random

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local
from check_noncommutative_crossing_symplectic_radicals import (
    closed_radical, omega, contraction,
)


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def transpose(matrix, size):
    return [sum(((matrix[i] >> j) & 1) << i for i in range(size)) for j in range(size)]


def compose_binary(first, second):
    return [xor_values(second[j] for j in bits(row)) for row in first]


def xor_values(values):
    result = 0
    for value in values:
        result ^= value
    return result


def left_binary(first, second):
    n = len(second)
    return [[xor_values(second[k][j] for k in bits(first[i])) for j in range(n)]
            for i in range(n)]


def right_binary(first, second):
    columns = transpose(second, len(second))
    return [[xor_values(first[i][k] for k in bits(columns[j])) for j in range(len(second))]
            for i in range(len(first))]


def local_clifford():
    matrices = []
    for i in range(2):
        columns = []
        for state in range(4):
            value = 0
            if not state >> i & 1:
                value ^= 1 << (state | (1 << i))
            if state >> (i ^ 1) & 1:
                value ^= 1 << (state ^ (1 << (i ^ 1)))
            columns.append(value)
        matrices.append(transpose(columns, 4))
    first, second = matrices
    assert compose_binary(first, first) == [0] * 4
    assert compose_binary(second, second) == [0] * 4
    assert compose_binary(first, second) == compose_binary(second, first)
    assert transpose(first, 4) == second
    q = compose_binary(first, second)
    assert q == [9, 0, 0, 9]
    return matrices


def verify_local_ward():
    stride = 8
    p, q, r = 1 << 1, 1 << stride, 1 << (stride * stride)
    eta = q ^ r
    local = evaluate_local((p, q, r, 1, 1, 1))
    family = {(x, u): [[local[((x * 4 + y) * 4 + u) * 4 + v] for v in range(4)]
                       for y in range(4)] for x, u in product(range(4), repeat=2)}
    count = 0
    for i, clifford in enumerate(local_clifford()):
        for x, u in product(range(4), repeat=2):
            matrix = family[x, u]
            left, right = left_binary(clifford, matrix), right_binary(matrix, clifford)
            expected = family[x ^ (1 << i), u] if x >> i & 1 else [[0] * 4 for _ in range(4)]
            for y, v in product(range(4), repeat=2):
                assert left[y][v] ^ right[y][v] == multiply(eta, expected[y][v])
                count += 1
    return {"exact_generic_polynomial_equations": count,
            "identity": "[c_i,M]=(q+r) A_i M", "Clifford_relations_checked": True}


def q_on_even_vector(vector, weight, width):
    result = 0
    parity = weight.bit_count() & 1
    for source in bits(vector):
        if parity:
            result ^= 1 << source
        for i in bits(weight):
            result ^= 1 << (source ^ (1 << i))
    return result


def image_intersection_dimensions(weights, width):
    columns = []
    for source in range(1 << width):
        vector = 1 << source
        for weight in weights:
            vector = q_on_even_vector(vector, weight, width)
        columns.append(vector)
    image_rank = rank_two(columns)
    result = []
    for grade in range(width + 1):
        inside = sum(1 << mask for mask in range(1 << width) if mask.bit_count() == grade)
        outside_rank = rank_two([vector & ~inside for vector in columns])
        result.append(image_rank - outside_rank)
    return image_rank, result


def css_dimension_audit(pairs, columns):
    forms = max(columns, default=0).bit_length()
    assert forms and rank_two(list(columns)) == forms
    totals = [0] * (2 * pairs + 1)
    syndromes = 0
    for sigma in range(1 << pairs):
        active = [i for i in range(pairs) if not sigma >> i & 1]
        t, s = len(active), sigma.bit_count()
        weights = [sum(((columns[i] >> row) & 1) << j for j, i in enumerate(active))
                   for row in range(forms)]
        independent = rank_two(weights) == forms
        gram_zero = all((a & b).bit_count() % 2 == 0 for a, b in product(weights, repeat=2))
        image_rank, actual = image_intersection_dimensions(weights, t)
        assert image_rank == ((1 << (t - forms)) if independent else 0)
        expected = [comb(t - 2 * forms, h - forms)
                    if independent and gram_zero and 0 <= h - forms <= t - 2 * forms
                    else 0 for h in range(t + 1)]
        assert actual == expected, (pairs, columns, sigma, actual, expected)
        for h, dimension in enumerate(actual):
            totals[s + 2 * h] += (1 << s) * dimension
        syndromes += 1
    fixed = [0] * len(totals)
    full_weights = [sum(((column >> row) & 1) << i for i, column in enumerate(columns))
                    for row in range(forms)]
    if all((a & b).bit_count() % 2 == 0 for a, b in product(full_weights, repeat=2)):
        for length in range(0, 2 * pairs + 1, 2):
            h = length // 2 - forms
            if 0 <= h <= pairs - 2 * forms:
                fixed[length] = comb(pairs - 2 * forms, h)
    assert all((d - f) % 2 == 0 for d, f in zip(totals, fixed))
    return {"pairs": pairs, "form_count": forms, "column_types": list(columns),
            "syndromes_checked": syndromes, "full_homogeneous_dimensions": totals,
            "J_fixed_dimensions": fixed,
            "J_half_dimensions": [(d - f) // 2 for d, f in zip(totals, fixed)]}


def q_full_vector(vector, pair_weight, pairs):
    result = 0
    for source in bits(vector):
        for i in bits(pair_weight):
            pair = (source >> (2 * i)) & 3
            if pair in (0, 3):
                result ^= 1 << source
                result ^= 1 << (source ^ (3 << (2 * i)))
    return result


def full_q_product(weights, pairs):
    size = 1 << (2 * pairs)
    columns = []
    for source in range(size):
        value = 1 << source
        for weight in weights:
            value = q_full_vector(value, weight, pairs)
        columns.append(value)
    rows = transpose(columns, size)
    assert rows == columns
    assert compose_binary(rows, rows) == [0] * size
    return rows


def verify_power_sandwiches():
    # p,r,eta occupy independent base-eight exponent digits; q=r+eta.
    p, r, eta = 2, 1 << 8, 1 << 64
    local = evaluate_local((p, r ^ eta, r, 1, 1, 1))
    rng = random.Random(20261030)
    reports = []
    for pairs, weights, count in ((2, (3,), 32), (2, (1, 2), 32),
                                  (3, (3, 6), 16)):
        size = 1 << (2 * pairs)
        projector = full_q_product(weights, pairs)
        nonzero = checked = 0
        for _ in range(count):
            x, u = rng.randrange(size), rng.randrange(size)
            matrix = [[0] * size for _ in range(size)]
            for y, v in product(range(size), repeat=2):
                value = 1
                for i in range(pairs):
                    a, b, c, d = ((word >> (2 * i)) & 3 for word in (x, y, u, v))
                    value = multiply(value, local[((a * 4 + b) * 4 + c) * 4 + d])
                matrix[y][v] = value
            result = right_binary(left_binary(projector, matrix), projector)
            forbidden = (1 << (64 * len(weights))) - 1
            for row in result:
                for value in row:
                    assert not value & forbidden
                    nonzero += bool(value)
                    checked += 1
        reports.append({"pairs": pairs, "form_weights": list(weights),
                        "external_coefficients_tested": count,
                        "polynomial_matrix_entries_checked": checked,
                        "nonzero_sandwich_entries": nonzero,
                        "guaranteed_q_plus_r_power": len(weights)})
    return reports


def verify_nested_collapse():
    controls = []
    for inner_pairs, inner_length in ((2, 2), (3, 3), (4, 4), (5, 3)):
        ambient, vectors, _, _ = closed_radical(inner_pairs, inner_length)
        pairs = inner_pairs + 2
        q_mask = 3 << (2 * (pairs - 1))
        for value in vectors[:8]:
            inner = {ambient[i] << 2 for i in bits(value)}
            actual = omega(contraction({mask | 3 for mask in inner}, pairs), pairs)
            expected = {mask | 3 for mask in inner} ^ {mask | q_mask for mask in inner}
            assert actual == expected
        controls.append({"inner_pairs": inner_pairs, "inner_grade": inner_length,
                         "basis_vectors_checked": min(8, len(vectors))})
    p, q, r = 2, 1 << 8, 1 << 64
    local = evaluate_local((p, q, r, 1, 1, 1))
    out, inputs = (0, 3), ((3, 0), (0, 3))
    scalar = 0
    for y, v in product(inputs, repeat=2):
        value = 1
        for a, b, c, d in zip(out, y, out, v):
            value = multiply(value, local[((a * 4 + b) * 4 + c) * 4 + d])
        scalar ^= value
    assert scalar == multiply(q ^ r, q ^ r)
    return {"vector_controls": controls, "two_pair_scalar_is_q_plus_r_squared": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    configurations = ((2, (1, 1)), (4, (1, 1, 1, 1)),
                      (4, (1, 1, 2, 2)), (5, (1, 1, 2, 3, 3)),
                      (6, (1, 1, 2, 2, 3, 3)), (6, (1, 1, 2, 2, 4, 4)),
                      (7, (1, 1, 1, 2, 2, 3, 3)),
                      (8, (1, 1, 1, 2, 2, 2, 3, 3)),
                      (7, (1, 2, 3, 4, 5, 6, 7)),
                      (6, (0, 1, 1, 2, 2, 3)))
    result = {
        "status": "Exact Ward identities, power sandwiches, homogeneous image dimensions, and nested collapse all passed.",
        "scope": "All-size characteristic-two crossing construction and factor theorem; no rank-sum improvement claimed.",
        "local_Ward": verify_local_ward(),
        "power_sandwich_controls": verify_power_sandwiches(),
        "weighted_form_dimension_controls": [css_dimension_audit(P, cols) for P, cols in configurations],
        "nested_scalar_collapse": verify_nested_collapse(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()

