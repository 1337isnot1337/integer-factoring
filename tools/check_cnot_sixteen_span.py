#!/usr/bin/env python3
"""Exact controls for the all-field sixteen-label CNOT span obstruction."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
import random

from check_mixed_cut_rank_budget import rank_binary


def flip_rows(q):
    return [1 << (q * j + i) for i in range(q) for j in range(q)]


def kron_binary(left, right, right_width):
    rows = []
    for a in left:
        for b in right:
            row = 0
            bits = a
            while bits:
                bit = bits & -bits
                row ^= b << ((bit.bit_length() - 1) * right_width)
                bits ^= bit
            rows.append(row)
    return rows


def multiply4(a, b):
    result = 0
    for _ in range(2):
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def rank4(matrix):
    a = [row[:] for row in matrix]
    pivot = 0
    for column in range(len(a[0])):
        found = next((i for i in range(pivot, len(a)) if a[i][column]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        inverse = [0, 1, 3, 2][a[pivot][column]]
        a[pivot] = [multiply4(x, inverse) for x in a[pivot]]
        for i in range(pivot + 1, len(a)):
            scalar = a[i][column]
            if scalar:
                a[i] = [x ^ multiply4(scalar, y) for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def field_control(matrix, scalar):
    q = len(matrix)
    full = [[int(i == l and j == k) ^
             multiply4(scalar, multiply4(matrix[i][k], matrix[j][l]))
             for k in range(q) for l in range(q)]
            for i in range(q) for j in range(q)]
    assert rank4(full) >= comb(q, 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(202610305164)
    counts = {'binary_minor_controls': 0, 'field4_minor_controls': 0,
              'field4_local_conics': 0, 'binary_four_site_conics': 0}

    for q in [2, 3, 4, 8, 16]:
        matrices = (range(1 << (q * q)) if q <= 3 else
                    [rng.getrandbits(q * q) for _ in range(64)])
        flip = flip_rows(q)
        for packed in matrices:
            matrix = [(packed >> (q * i)) & ((1 << q) - 1) for i in range(q)]
            tensor = kron_binary(matrix, matrix, q)
            assert rank_binary([a ^ b for a, b in zip(flip, tensor)]) >= comb(q, 2)
            counts['binary_minor_controls'] += 1

    for q in [2, 3, 4]:
        words = (product(range(4), repeat=4) if q == 2 else
                 [tuple(rng.randrange(4) for _ in range(q * q)) for _ in range(32)])
        for word in words:
            matrix = [list(word[q * i:q * (i + 1)]) for i in range(q)]
            for scalar in range(4):
                field_control(matrix, scalar)
                counts['field4_minor_controls'] += 1

    packed_b, packed_c, packed_d = [0, 12, 10, 6], [0, 0, 12, 12], [0, 10, 0, 10]
    for x, y in product(range(4), repeat=2):
        if x == y == 0:
            continue
        vector = [0, y, x, x ^ y]
        matrix = [[0, y], [x, x ^ y]]
        row = [[multiply4(vector[i], vector[j]) for j in range(4)] for i in range(4)]
        for i, j in product(range(4), repeat=2):
            expected = 0
            for coefficient, basis in [(multiply4(x, y), packed_b),
                                       (multiply4(x, x), packed_c),
                                       (multiply4(y, y), packed_d)]:
                if (basis[i] >> j) & 1:
                    expected ^= coefficient
            assert row[i][j] == expected
        column = [[row[2 * i + k][2 * j + l] for k in range(2) for l in range(2)]
                  for i in range(2) for j in range(2)]
        expected = [[multiply4(matrix[i][k], matrix[j][l])
                     for k in range(2) for l in range(2)]
                    for i in range(2) for j in range(2)]
        assert column == expected
        assert rank4(row) == 1
        assert rank4(column) == (4 if x and y else 1)
        counts['field4_local_conics'] += 1

    conic_column_ranks = {}
    for parameters in product([(1, 0), (0, 1), (1, 1)], repeat=4):
        matrix, width = [1], 1
        active = 0
        for x, y in parameters:
            local = [y << 1, x | ((x ^ y) << 1)]
            matrix = kron_binary(matrix, local, 2)
            width *= 2
            active += bool(x and y)
        assert width == 16
        column = kron_binary(matrix, matrix, width)
        observed = rank_binary(column)
        assert observed == 4 ** active
        assert rank_binary([a ^ b for a, b in zip(flip_rows(width), column)]) >= 120
        conic_column_ranks[observed] = conic_column_ranks.get(observed, 0) + 1
        counts['binary_four_site_conics'] += 1

    cases = {'neither_exception': 512 + 11 * 16,
             'one_exception': 16 * 13 + 120 + 255 + 11 * 15,
             'both_exceptions': 32 * 13 + 119 + 14 + 11 * 14}
    assert cases == {'neither_exception': 688, 'one_exception': 748, 'both_exceptions': 703}
    assert min(cases.values()) ** 70 > 16 ** 99 * 13 ** 70
    assert 99 ** 2 > 2 * 70 ** 2
    evidence = {'scope': 'All characteristic-two fields; full V tensor fourth power; no claim for changed coefficient spaces or larger powers.',
                'counts': counts, 'four_site_conic_column_rank_histogram': conic_column_ranks,
                'case_bound_numerators': cases, 'case_bound_denominator': 13,
                'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    rendered = json.dumps(evidence, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
