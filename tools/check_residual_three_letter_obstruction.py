#!/usr/bin/env python3
"""Exact controls for the all-field residual three-letter code obstruction."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_residual_quartic_codes import rank, residual_operator


ROOT = Path(__file__).resolve().parent.parent
MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
INV = (0, 1, 3, 2)


def field_rank(matrix):
    matrix = [list(row) for row in matrix]
    rows, columns = len(matrix), len(matrix[0])
    result = 0
    for column in range(columns):
        pivot = next((i for i in range(result, rows) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[result], matrix[pivot] = matrix[pivot], matrix[result]
        scale = INV[matrix[result][column]]
        matrix[result] = [MUL[value][scale] for value in matrix[result]]
        for i in range(rows):
            if i == result or not matrix[i][column]:
                continue
            scale = matrix[i][column]
            matrix[i] = [a ^ MUL[scale][b] for a, b in zip(matrix[i], matrix[result])]
        result += 1
        if result == rows:
            break
    return result


def matrix_space():
    positions = [(0, 5, 10), (15,), (1, 11), (8, 13), (9,), (6,), (7,), (14,)]
    return [sum(1 << p for p in support) for support in positions]


def transpose(matrix):
    return sum(((matrix >> (4 * i + j)) & 1) << (4 * j + i)
               for i, j in product(range(4), repeat=2))


def constraints(a, b):
    # Each entry functional A_ij acts on (r0,...,r3,s0,...,s3).
    def entry(i, j):
        result = [0] * 8
        result[j] = a[i]
        result[4 + i] = b[j]
        return result

    def difference(i, j, k, ell):
        return [x ^ y for x, y in zip(entry(i, j), entry(k, ell))]

    return [entry(0, 2), entry(0, 3), entry(1, 0), entry(3, 0),
            difference(0, 0, 1, 1), difference(0, 0, 2, 2),
            difference(0, 1, 2, 3), difference(2, 0, 3, 1)]


def projective_vectors(field):
    result = []
    for pivot in range(4):
        for tail in product(range(field), repeat=3 - pivot):
            result.append((0,) * pivot + (1,) + tail)
    return result


def check():
    operator = residual_operator()
    def swapped(vector):
        return sum(((vector >> (4 * i + j)) & 1) << (4 * j + i)
                   for i, j in product(range(4), repeat=2))
    def applied(vector):
        return sum(((row & vector).bit_count() & 1) << i for i, row in enumerate(operator))
    for j in range(16):
        assert swapped(applied(1 << j)) == applied(swapped(1 << j))
    quotient = [[(operator[5 * i] >> (5 * j)) & 1 for j in range(4)] for i in range(4)]
    assert quotient == [[1, 0, 0, 0], [0, 1, 1, 0], [0, 1, 1, 0], [0, 0, 0, 0]]
    assert field_rank(quotient) == 2
    crossing_cut = [0] * 16
    for x, y, u, v in product(range(4), repeat=4):
        if operator[4 * y + v] >> (4 * u + x) & 1:
            crossing_cut[4 * x + v] |= 1 << (4 * y + u)
    columns = [sum(((row >> j) & 1) << i for i, row in enumerate(crossing_cut)) for j in range(16)]
    W = matrix_space()
    WT = [transpose(matrix) for matrix in W]
    assert rank(crossing_cut) == rank(W) == rank(WT) == 8
    assert rank(W + columns) == 8
    assert rank(WT + crossing_cut) == 8
    permutation = (0, 2, 1, 3)
    transformed = [sum(((matrix >> (4 * permutation[i] + permutation[j])) & 1) << (4 * i + j)
                       for i, j in product(range(4), repeat=2)) for matrix in WT]
    assert rank(W + transformed) == 8
    records = []
    for field in (2, 4):
        vectors = projective_vectors(field)
        counts = Counter()
        for a, b in product(vectors, repeat=2):
            r = field_rank(constraints(a, b))
            assert r >= 4
            counts[r] += 1
        records.append({'field': field, 'projective_vectors': len(vectors),
                        'ordered_pairs_checked': len(vectors) ** 2,
                        'equation_rank_distribution': dict(sorted(counts.items())),
                        'minimum_compression_dimension': min(r + 1 for r in counts)})
    return {'status': 'exact controls for an all-characteristic-two-field theorem',
            'crossing_cut_rank': 8, 'matrix_space_dimension': 8,
            'minimum_compressed_crossing_cut_rank': 2,
            'factor_exchange_commutation': True,
            'invariant_alternating_quotient_matrix': quotient,
            'invariant_alternating_quotient_rank': 2,
            'projective_equation_checks': records,
            'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                              for name in ('tools/check_residual_three_letter_obstruction.py',
                                           'tools/check_residual_quartic_codes.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rendered = json.dumps(check(), indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
