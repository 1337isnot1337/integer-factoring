#!/usr/bin/env python3
"""Exact finite controls for the all-field two-product complement theorem.

The proof is the common-input/output-flag induction and triangular kernel
bound in two-product-complement.md. Finite controls do not prove all d.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import random


def rank(rows, p):
    a = [list(row) for row in rows]
    if not a:
        return 0
    answer = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(answer, len(a)) if a[i][col] % p), None)
        if pivot is None:
            continue
        a[answer], a[pivot] = a[pivot], a[answer]
        scale = pow(a[answer][col] % p, -1, p)
        a[answer] = [x * scale % p for x in a[answer]]
        for i in range(answer + 1, len(a)):
            scale = a[i][col] % p
            if scale:
                a[i] = [(x - scale * y) % p for x, y in zip(a[i], a[answer])]
        answer += 1
        if answer == len(a):
            break
    return answer


def rank2(rows):
    pivots = {}
    for row in rows:
        while row:
            position = row.bit_length() - 1
            if position in pivots:
                row ^= pivots[position]
            else:
                pivots[position] = row
                break
    return len(pivots)


def identity(d):
    return [[int(i == j) for j in range(d)] for i in range(d)]


def matrix_product(a, b, p):
    return [[sum(x * y for x, y in zip(row, column)) % p
             for column in zip(*b)] for row in a]


def inverse(a, p):
    d = len(a)
    b = [row[:] + unit for row, unit in zip(a, identity(d))]
    for col in range(d):
        pivot = next(i for i in range(col, d) if b[i][col] % p)
        b[col], b[pivot] = b[pivot], b[col]
        scale = pow(b[col][col], -1, p)
        b[col] = [x * scale % p for x in b[col]]
        for i in range(d):
            if i != col:
                scale = b[i][col]
                b[i] = [(x - scale * y) % p for x, y in zip(b[i], b[col])]
    result = [row[d:] for row in b]
    assert matrix_product(a, result, p) == identity(d)
    return result


def random_matrix(d, p, rng, triangular=False):
    return [[rng.randrange(p) if not triangular or i <= j else 0
             for j in range(d)] for i in range(d)]


def random_invertible(d, p, rng):
    while True:
        a = random_matrix(d, p, rng)
        if rank(a, p) == d:
            return a


def kronecker(a, b, p):
    return [[a[i][j] * b[u][v] % p
             for j in range(len(a[0])) for v in range(len(b[0]))]
            for i in range(len(a)) for u in range(len(b))]


def add_matrices(matrices, p):
    return [[sum(entries) % p for entries in zip(*rows)] for rows in zip(*matrices)]


def subtract(a, b, p):
    return [[(x - y) % p for x, y in zip(row, other)] for row, other in zip(a, b)]


def swap(d):
    return [[int(i == v and u == j) for j in range(d) for v in range(d)]
            for i in range(d) for u in range(d)]


def other_cut(a, d):
    return [[a[x * d + y][u * d + v] for y in range(d) for v in range(d)]
            for x in range(d) for u in range(d)]


def bit_rows(a):
    return [sum((x & 1) << j for j, x in enumerate(row)) for row in a]


def matrix_from_word(word, d):
    return [[(word >> (i * d + j)) & 1 for j in range(d)] for i in range(d)]


def exhaustive_two_label():
    d = 2
    e = bit_rows(swap(d))
    counts = Counter()
    minima = {}
    examples = {}
    for word in range(1 << 16):
        a = matrix_from_word(word, 4)
        c = rank2(bit_rows(other_cut(a, d)))
        if c <= 2:
            rr = rank2([x ^ y for x, y in zip(e, bit_rows(a))])
            assert rr >= d * (d - 1) // 2
            counts[c] += 1
            if c not in minima or rr < minima[c]:
                minima[c], examples[c] = rr, word
    assert sum(counts.values()) == 7576
    assert counts == {0: 1, 1: 225, 2: 7350}
    return dict(field=2, labels=d, all_tensors=1 << 16,
                tested_column_rank_at_most_two=sum(counts.values()),
                counts=dict(counts), minimum_complement_ranks=minima,
                minimizing_encoded_tensors=examples)


def exhaustive_singular_pencil():
    # Kronecker right and left minimal-index blocks, embedded in a square
    # pencil. Every combination is singular, but common kernel is zero
    # and the joint image is the whole space.
    d = 3
    p = [[1, 0, 0], [0, 0, 1], [0, 0, 0]]
    r = [[0, 1, 0], [0, 0, 0], [0, 0, 1]]
    assert all(rank2(bit_rows(add_matrices([[[a * x for x in row] for row in p],
                                           [[b * x for x in row] for row in r]], 2))) < d
               for a, b in product(range(2), repeat=2))
    assert rank2(bit_rows(p) + bit_rows(r)) == d
    assert rank2(bit_rows(list(map(list, zip(*p)))) +
                 bit_rows(list(map(list, zip(*r))))) == d
    table_p, table_r = [], []
    for word in range(1 << 9):
        matrix = matrix_from_word(word, d)
        table_p.append(bit_rows(kronecker(p, matrix, 2)))
        table_r.append(bit_rows(kronecker(r, matrix, 2)))
    e = bit_rows(swap(d))
    histogram = Counter()
    best, witness = d * d, None
    for q, left in enumerate(table_p):
        for s, right in enumerate(table_r):
            rr = rank2([x ^ y ^ z for x, y, z in zip(e, left, right)])
            assert rr >= d * (d - 1) // 2
            histogram[rr] += 1
            if rr < best:
                best, witness = rr, [q, s]
    assert sum(histogram.values()) == 1 << 18
    return dict(field=2, labels=d, all_second_factor_pairs=1 << 18,
                all_pencil_members_singular=True, common_kernel_dimension=0,
                joint_image_dimension=d, minimum_complement_rank=best,
                minimizing_second_factors=witness, rank_histogram=dict(histogram))


def random_checks():
    rng = random.Random(590011)
    arbitrary = triangular = transforms = sharp = 0
    minimum_margins = {}
    for p in (2, 3, 5, 7):
        for d in range(1, 8):
            threshold = d * (d - 1) // 2
            minimum = d * d
            for _ in range(50):
                ps, qs, rs, ss = [random_matrix(d, p, rng) for _ in range(4)]
                approximation = add_matrices([kronecker(ps, qs, p),
                                              kronecker(rs, ss, p)], p)
                assert rank(other_cut(approximation, d), p) <= 2
                rr = rank(subtract(swap(d), approximation, p), p)
                assert rr >= threshold
                minimum = min(minimum, rr)
                arbitrary += 1
            minimum_margins[f'{p}:{d}'] = minimum - threshold
            for count in (1, 2, 3, 4):
                for _ in range(10):
                    pairs = [(random_matrix(d, p, rng, triangular=True),
                              random_matrix(d, p, rng)) for _ in range(count)]
                    approximation = add_matrices([kronecker(a, b, p) for a, b in pairs], p)
                    assert rank(subtract(swap(d), approximation, p), p) >= threshold
                    triangular += 1
            assert rank(subtract(swap(d), kronecker(identity(d), identity(d), p), p), p) == threshold
            sharp += 1
        for d in (2, 3, 4):
            for _ in range(12):
                u, v = random_invertible(d, p, rng), random_invertible(d, p, rng)
                ui, vi = inverse(u, p), inverse(v, p)
                left, right = kronecker(u, v, p), kronecker(vi, ui, p)
                assert matrix_product(matrix_product(left, swap(d), p), right, p) == swap(d)
                ps, qs = random_matrix(d, p, rng), random_matrix(d, p, rng)
                transformed_p = matrix_product(matrix_product(u, ps, p), vi, p)
                transformed_q = matrix_product(matrix_product(v, qs, p), ui, p)
                transformed_product = matrix_product(matrix_product(left, kronecker(ps, qs, p), p), right, p)
                assert transformed_product == kronecker(transformed_p, transformed_q, p)
                transforms += 1
    return dict(arbitrary_two_product_cases=arbitrary,
                triangular_first_factor_cases=triangular,
                swap_preserving_equivalences=transforms,
                sharpness_cases=sharp, minimum_random_margin=minimum_margins)


def negative_three_product_controls():
    controls = []
    # For d=3, B=E-|Omega><Omega| has opposite-cut rank three and a
    # rank-one complement. Thus the same half-dimensional theorem cannot
    # be extended to three arbitrary products.
    d = 3
    omega = [int(i == j) for i in range(d) for j in range(d)]
    residual = [[x * y for y in omega] for x in omega]
    for p in (2, 3, 5, 7):
        approximation = subtract(swap(d), residual, p)
        products = []
        for i in range(d):
            for j in range(i + 1, d):
                skew = [[int(a == i and b == j) - int(a == j and b == i)
                         for b in range(d)] for a in range(d)]
                products.append(kronecker(skew, skew, p))
        three_products = [[-x % p for x in row] for row in add_matrices(products, p)]
        assert approximation == three_products
        c = rank(other_cut(approximation, d), p)
        rr = rank(residual, p)
        assert c == 3 and rr == 1 < d * (d - 1) // 2
        controls.append(dict(field=p, labels=d, column_rank=c,
                             complement_rank=rr,
                             explicit_three_product_identity=True,
                             proposed_unrestricted_three_product_bound_fails=True))
    return controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='PASS',
                  scope='All-field two-product complement interface; no new factoring exponent.',
                  exhaustive_two_label=exhaustive_two_label(),
                  exhaustive_singular_pencil=exhaustive_singular_pencil(),
                  random_and_structural=random_checks(),
                  three_product_negative_controls=negative_three_product_controls(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
