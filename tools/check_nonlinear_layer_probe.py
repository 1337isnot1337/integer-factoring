#!/usr/bin/env python3
"""Exact nonlinear crossing-layer control over F2; no improved bound."""

from __future__ import annotations

import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path


def polynomial_product(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def entry_polynomial(width, labels, permutation=None):
    if permutation is not None:
        labels = tuple(permutation[x] for x in labels)
    answer = 1
    for j in range(width):
        x, y, u, v = ((x >> j) & 1 for x in labels)
        crossing = int(x == v and y == u)
        identity = int(x == u and y == v)
        answer = polynomial_product(answer, identity | ((crossing ^ identity) << 1))
    return answer


def pencil(width):
    d = 1 << width
    answer = [[0] * (d * d) for _ in range(width + 1)]
    for x in range(d):
        for y in range(d):
            for u in range(d):
                for v in range(d):
                    poly = entry_polynomial(width, (x, y, u, v))
                    while poly:
                        bit = poly & -poly
                        answer[bit.bit_length() - 1][x * d + y] |= 1 << (u * d + v)
                        poly ^= bit
    return answer


def xor_rows(a, b):
    return [x ^ y for x, y in zip(a, b)]


def multiply_rows(a, b):
    result = []
    for row in a:
        value = 0
        while row:
            bit = row & -row
            value ^= b[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return result


def multiply_pencils(a, b):
    result = [[0] * len(a[0]) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = xor_rows(result[i + j], multiply_rows(ai, bj))
    while len(result) > 1 and not any(result[-1]):
        result.pop()
    return result


def relabel(rows, permutation, d):
    pairs = [permutation[x] * d + permutation[y] for x in range(d) for y in range(d)]
    return [sum(((rows[pairs[i]] >> pairs[j]) & 1) << j for j in range(d * d))
            for i in range(d * d)]


def leg_permutation(rows, permutation, d):
    result = [0] * (d * d)
    for x in range(d):
        for y in range(d):
            for u in range(d):
                for v in range(d):
                    labels = (x, y, u, v)
                    a, b, c, e = (labels[j] for j in permutation)
                    result[x * d + y] |= ((rows[a * d + b] >> (c * d + e)) & 1) << (u * d + v)
    return result


def matrix_rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def encode(rows):
    n = len(rows)
    return sum(row << (n * i) for i, row in enumerate(rows))


def decode(vector, n):
    mask = (1 << n) - 1
    return [(vector >> (n * i)) & mask for i in range(n)]


def insert(basis, vector, coordinate):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot in basis:
            other, mask = basis[pivot]
            vector ^= other
            coordinate ^= mask
        else:
            basis[pivot] = (vector, coordinate)
            return True
    return False


def reduce_vector(basis, vector):
    for pivot in sorted(basis, reverse=True):
        if vector >> pivot & 1:
            vector ^= basis[pivot][0]
    return vector


def represent(basis, vector):
    coordinate = 0
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in basis:
            return None
        vector ^= basis[pivot][0]
        coordinate ^= basis[pivot][1]
    return coordinate


def non_gaussian_control(width, permutation, indices):
    d = 1 << width

    def value(subset):
        mask = sum(1 << j for j in subset)
        labels = tuple((mask >> (leg * width)) & (d - 1) for leg in range(4))
        return entry_polynomial(width, labels, permutation)

    assert value(()) == 1
    i, j, k, ell = indices
    pfaffian = (polynomial_product(value((i, j)), value((k, ell)))
                ^ polynomial_product(value((i, k)), value((j, ell)))
                ^ polynomial_product(value((i, ell)), value((j, k))))
    actual = value(indices)
    assert actual ^ pfaffian == 6  # z+z^2, nonzero over F2(z).
    return {"width": width, "indices": list(indices),
            "actual_polynomial_bits": actual, "pfaffian_polynomial_bits": pfaffian,
            "defect_polynomial_bits": 6}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cnot = [x ^ ((x & 1) << 1) for x in range(4)]
    toffoli = [x ^ (4 if x & 3 == 3 else 0) for x in range(8)]
    gaussian = [non_gaussian_control(2, cnot, (0, 1, 3, 4)),
                non_gaussian_control(3, toffoli, (0, 1, 6, 10))]
    d, n = 8, 64
    p = pencil(3)
    q = [relabel(rows, toffoli, d) for rows in p]
    identity = [1 << i for i in range(n)]
    assert multiply_pencils(p, p) == [identity]
    assert multiply_pencils(q, q) == [identity]
    pq, qp = multiply_pencils(p, q), multiply_pencils(q, p)
    assert (pq[2][3] ^ qp[2][3]) >> 10 & 1
    t = multiply_pencils(p, multiply_pencils(q, p))
    assert len(t) == 5
    correction = t[4]
    assert all(xor_rows(t[i], q[i] if i < len(q) else [0] * n)
               == (correction if i in (2, 4) else [0] * n) for i in range(5))
    dual = lambda rows: leg_permutation(rows, (0, 2, 1, 3), d)
    assert (matrix_rank(correction), matrix_rank(dual(correction))) == (2, 8)
    crossing = [1 << (y * d + x) for x in range(d) for y in range(d)]
    total = [0] * n
    for rows in t:
        total = xor_rows(total, rows)
    assert total == crossing
    symmetries = [perm for perm in permutations(range(4))
                  if {frozenset((perm[0], perm[3])), frozenset((perm[1], perm[2]))}
                  == {frozenset((0, 3)), frozenset((1, 2))}]
    assert len(symmetries) == 8
    raw = {}
    for rows in p + q:
        for symmetry in symmetries:
            insert(raw, encode(leg_permutation(rows, symmetry, d)), 0)
    assert len(raw) == 11
    assert reduce_vector(raw, encode(correction)) != 0
    for rows in t:
        for symmetry in symmetries:
            insert(raw, encode(leg_permutation(rows, symmetry, d)), 0)
    basis = [vector for vector, _ in raw.values()]
    assert len(basis) == 13
    coordinates = {}
    for i, vector in enumerate(basis):
        insert(coordinates, vector, 1 << i)
    target = represent(coordinates, encode(crossing))
    assert target is not None
    dual_basis = [encode(dual(decode(vector, n))) for vector in basis]
    row_vectors = [0] * (1 << 13)
    col_vectors = [0] * (1 << 13)
    strict, passing = {}, {}
    ranks = []
    for coordinate in range(1, 1 << 13):
        bit = coordinate & -coordinate
        j = bit.bit_length() - 1
        row_vectors[coordinate] = row_vectors[coordinate ^ bit] ^ basis[j]
        col_vectors[coordinate] = col_vectors[coordinate ^ bit] ^ dual_basis[j]
        r = matrix_rank(decode(row_vectors[coordinate], n))
        c = matrix_rank(decode(col_vectors[coordinate], n))
        ranks.append((r, c))
        if r * c < 224:
            insert(strict, coordinate, 0)
        if r * c <= 224:
            insert(passing, coordinate, 0)
    assert reduce_vector(strict, target) != 0
    assert reduce_vector(passing, target) == 0
    pivot = reduce_vector(strict, target).bit_length() - 1
    separator = sum(((reduce_vector(strict, 1 << j) >> pivot) & 1) << j for j in range(13))
    assert (separator & target).bit_count() % 2 == 1
    assert all((separator & coordinate).bit_count() % 2 == 0
               for coordinate, (r, c) in enumerate(ranks, 1) if r * c < 224)
    result = {
        "status": "PASS", "scope": "nonlinear-layer finite probe over F2 only",
        "improved_factoring_exponent": False,
        "non_gaussian_controls": gaussian,
        "toffoli_permutation": toffoli,
        "nonzero_degree2_commutator_entry": {"row_labels": [0, 3], "column_labels": [1, 2]},
        "palindromic_product_degree": 4,
        "correction_rank_pair": [2, 8],
        "identity": "P Q P = Q + (z^2+z^4) D",
        "source_and_relabeling_span_dimension": 11,
        "correction_is_new_to_that_span": True,
        "coefficient_symmetry_span_dimension": 13,
        "nonzero_vectors_ranked": 8191,
        "exact_optimal_maximum_rank_product": 224,
        "strictly_cheaper_span_dimension": len(strict),
        "passing_span_dimension": len(passing),
        "crossing_coordinates": target, "separator_coordinates": separator,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
