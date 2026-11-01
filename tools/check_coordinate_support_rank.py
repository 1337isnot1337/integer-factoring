#!/usr/bin/env python3
"""Exact finite checks of the coordinate-support rank-product theorem."""

from collections import Counter
from itertools import product
from math import comb
from random import Random
import json

from check_partition_channels import rank_mod, column_matrix


def supported(x, y, u, v, mask):
    if x ^ y ^ u ^ v:
        return False
    a, b = x ^ y, x ^ v
    return not (b & (a ^ mask))


def diagonal_count(matrix, m):
    return sum(bool(matrix[x*m+y][x*m+y])
               for x in range(m) for y in range(m))


def theorem_check(matrix, labels, mask, prime):
    m = len(labels)
    for xi, yi, vi, ui in product(range(m), repeat=4):
        if matrix[xi*m+yi][vi*m+ui] % prime:
            assert supported(labels[xi], labels[yi], labels[ui], labels[vi], mask)
    row_rank = rank_mod(matrix, prime)
    col_rank = rank_mod(column_matrix(matrix, m), prime)
    count = diagonal_count(matrix, m)
    assert m*row_rank*col_rank >= count*count
    return row_rank, col_rank, count


def local_polynomial(x, y, u, v, length):
    """Direct product of local R+z(E+R), packed over F_2."""
    value = 1
    for j in range(length):
        xx, yy, uu, vv = ((z >> j) & 1 for z in (x, y, u, v))
        rr = int(xx == yy and vv == uu)
        cc = int(xx == vv and yy == uu) ^ rr
        value = (value if rr else 0) ^ (value << 1 if cc else 0)
    return value


def coefficient(x, y, u, v, i):
    a, b = x ^ y, x ^ v
    k = i-a.bit_count()
    if x ^ y ^ u ^ v or a & b or not 0 <= k <= b.bit_count():
        return 0
    return comb(b.bit_count(), k) % 2


def family_check(labels, length):
    m = len(labels)
    matrices = [[[0]*(m*m) for _ in range(m*m)] for _ in range(length+1)]
    entries = 0
    for xi, yi, vi, ui in product(range(m), repeat=4):
        x, y, u, v = labels[xi], labels[yi], labels[ui], labels[vi]
        packed = local_polynomial(x, y, u, v, length)
        assert packed.bit_count() % 2 == int(x == v and y == u)
        for i, matrix in enumerate(matrices):
            entry = coefficient(x, y, u, v, i)
            assert entry == (packed >> i) & 1
            matrix[xi*m+yi][vi*m+ui] = entry
            entries += 1
    profile = [theorem_check(matrix, labels, 0, 2)[:2] for matrix in matrices]
    nonzero = sum(bool(a) for a, b in profile)
    maximum = max(a*b for a, b in profile)
    assert nonzero*nonzero*maximum >= m**3
    distances = Counter((a ^ b).bit_count() for a in labels for b in labels)
    for i, pairs in distances.items():
        a, b = profile[i]
        assert m*a*b >= pairs*pairs
    return {
        'labels': labels, 'length': length,
        'profile': profile, 'nonzero_pieces': nonzero,
        'maximum_rank_product': maximum,
        'distance_multiplicities': dict(sorted(distances.items())),
    }, matrices, entries


def verify():
    rng = Random(30102026)
    cases = [([0, 1], 1), ([0, 1, 2, 3], 2),
             ([0, 1, 3, 6], 3), ([0, 3, 5, 6, 9], 4),
             (list(range(8)), 3), ([0, 3, 12, 15, 48, 51, 60, 63], 6),
             ([0, 3, 12, 15, 48, 51, 60, 63, 64, 67, 76, 79, 112, 115, 124, 127], 7)]
    records, families, entries = [], [], 0
    for labels, length in cases:
        record, matrices, count = family_check(labels, length)
        records.append(record)
        families.append((labels, matrices))
        entries += count

    supported_cases = 0
    for labels, length in cases[:5]:
        m = len(labels)
        for prime in (2, 3, 5):
            for _ in range(6):
                mask = rng.randrange(1 << length)
                matrix = [[0]*(m*m) for _ in range(m*m)]
                for xi, yi, vi, ui in product(range(m), repeat=4):
                    if supported(labels[xi], labels[yi], labels[ui], labels[vi], mask):
                        matrix[xi*m+yi][vi*m+ui] = rng.randrange(prime)
                theorem_check(matrix, labels, mask, prime)
                supported_cases += 1

    recombinations = 0
    for labels, matrices in families[:6]:
        m = len(labels)
        for _ in range(8):
            selected = [matrix for matrix in matrices if rng.randrange(2)]
            matrix = [[sum(a[i][j] for a in selected) % 2
                       for j in range(m*m)] for i in range(m*m)]
            theorem_check(matrix, labels, 0, 2)
            recombinations += 1

    labels, m = [0, 1], 2
    mixed = [[int(not (x ^ y ^ u ^ v)) for v in labels for u in labels]
             for x in labels for y in labels]
    rr, cc = rank_mod(mixed, 2), rank_mod(column_matrix(mixed, m), 2)
    count = diagonal_count(mixed, m)
    assert (rr, cc, count) == (2, 2, 4)
    assert m*rr*cc < count*count
    assert any(mixed[x*m+y][v*m+u] and not supported(x, y, u, v, 0)
               for x, y, u, v in product(labels, repeat=4))
    assert any(mixed[x*m+y][v*m+u] and not supported(x, y, u, v, 1)
               for x, y, u, v in product(labels, repeat=4))

    return {
        'scope': 'exact finite support/rank and coefficient checks; the note proves the all-size bound',
        'supported_tensors': supported_cases,
        'within_support_recombinations': recombinations,
        'direct_product_coefficient_entries': entries,
        'families': records,
        'mixed_support_negative_control': {'ranks': [rr, cc], 'diagonal_entries': count},
        'general_bound': 'max rank_row*rank_col >= m^3/M^2',
        'asymptotic_restriction': 'rho >= 3/4 when log M=o(log m)',
        'new_factoring_exponent_claimed': False,
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
