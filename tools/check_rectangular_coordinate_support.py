#!/usr/bin/env python3
"""Exact finite controls for unequal-alphabet coordinate support and rank sums."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_coordinate_support_rank import supported
from check_partition_channels import rank_mod
from physical_split_crossings import F4_PRODUCT
from weighted_smith_oracle_f4 import polynomial_multiply


def local_support(x, y, u, v, mask, length):
    for bit in range(length):
        a, b, c, d = ((z >> bit) & 1 for z in (x, y, u, v))
        crossing = a == d and b == c
        alternative = (a == c and b == d) if mask >> bit & 1 else (a == b and c == d)
        if not crossing and not alternative:
            return False
    return True


def rank_f4(matrix):
    """Restriction of scalars, followed by binary bit-row elimination."""
    rows = []
    for row in matrix:
        low = high = 0
        for j, c in enumerate(row):
            omega = F4_PRODUCT[c][2]
            low |= (c & 1) << (2*j) | (omega & 1) << (2*j+1)
            high |= (c >> 1) << (2*j) | (omega >> 1) << (2*j+1)
        rows.extend((low, high))
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    assert len(pivots) % 2 == 0
    return len(pivots)//2


def check(matrix, alphabet_x, alphabet_y, mask, length, field):
    a, b = len(alphabet_x), len(alphabet_y)
    opposite = [[0]*(a*b) for _ in range(a*b)]
    for xi, yi, vi, ui in product(range(a), range(b), range(a), range(b)):
        value = matrix[xi*b+yi][vi*b+ui]
        if value:
            assert local_support(alphabet_x[xi], alphabet_y[yi], alphabet_y[ui], alphabet_x[vi], mask, length)
        opposite[xi*b+ui][vi*b+yi] = value
    rank = rank_f4 if field == 4 else lambda x: rank_mod(x, field)
    rr, cc = rank(matrix), rank(opposite)
    count = sum(bool(matrix[i][i]) for i in range(a*b))
    assert min(a, b)*rr*cc >= count*count
    return {'rank_row': rr, 'rank_column': cc, 'crossing_diagonal_count': count,
            'integer_slack': min(a, b)*rr*cc-count*count}


def evaluate_one(poly):
    result = 0
    while poly:
        result ^= poly & 3
        poly >>= 2
    return result


def diagonal_transpose_product(diagonals):
    """Literal operator product of I+(1+z) S_D; no profile formula is used."""
    q = len(diagonals[0])
    matrix = [[int(i == j) for j in range(q*q)] for i in range(q*q)]
    inverse = (0, 1, 3, 2)
    for diagonal in diagonals:
        old = matrix
        matrix = [[0]*(q*q) for _ in range(q*q)]
        for row, v, u in product(range(q*q), range(q), range(q)):
            c = F4_PRODUCT[diagonal[u]][inverse[diagonal[v]]]
            factor = c ^ (c << 2)
            matrix[row][v*q+u] = old[row][v*q+u] ^ polynomial_multiply(factor, old[row][u*q+v])
    for x, y, v, u in product(range(q), repeat=4):
        value = matrix[x*q+y][v*q+u]
        assert evaluate_one(value) == int(x == v and y == u)
        if value:
            assert (x == v and y == u) or (x == u and y == v)
    return matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, cases, families = Random(121220261102), [], []
    truth_checks = 0
    for mask, x, y, u, v in product(range(2), range(2), range(2), range(2), range(2)):
        assert supported(x, y, u, v, mask) == local_support(x, y, u, v, mask, 1)
        truth_checks += 1
    for length, field in product(range(1, 5), (2, 3, 5)):
        available = list(range(1 << length))
        for _ in range(50):
            a = rng.randrange(1, min(6, len(available))+1)
            b = rng.randrange(1, min(7, len(available))+1)
            if a == b:
                b = b-1 if b > 1 else 2
            ax, ay = sorted(rng.sample(available, a)), sorted(rng.sample(available, b))
            mask = rng.randrange(1 << length)
            matrix = [[0]*(a*b) for _ in range(a*b)]
            for x, y, v, u in product(range(a), range(b), range(a), range(b)):
                if supported(ax[x], ay[y], ay[u], ax[v], mask):
                    matrix[x*b+y][v*b+u] = rng.randrange(field)
            result = check(matrix, ax, ay, mask, length, field)
            cases.append(dict(result, alphabet_x=ax, alphabet_y=ay, mask=mask, field=field))
    for q, degree in ((2, 1), (3, 3), (4, 5), (5, 7)):
        diagonals = [[rng.randrange(1, 4) for _ in range(q)] for _ in range(degree)]
        source = diagonal_transpose_product(diagonals)
        length = (q-1).bit_length()
        ax, ay = list(range(max(1, q-2))), list(range(q))
        coefficient_checks, total = [], 0
        for power in range(degree+1):
            matrix = [[(source[x*q+y][v*q+u] >> (2*power)) & 3
                       for v in ax for u in ay] for x in ax for y in ay]
            result = check(matrix, ax, ay, (1 << length)-1, length, 4)
            coefficient_checks.append(result)
            total += result['crossing_diagonal_count']
        assert total >= len(ax)*len(ay)
        families.append({'q': q, 'degree': degree, 'diagonals': diagonals,
                         'alphabet_x': ax, 'alphabet_y': ay, 'coefficients': coefficient_checks,
                         'integer_crossing_diagonal_coverage': total})
    # Mixing the two alternative pairings violates the support hypothesis.
    mixed = [[int(not (x ^ y ^ u ^ v)) for v in range(2) for u in range(2)]
             for x in range(2) for y in range(2)]
    opposite = [[mixed[x*2+y][v*2+u] for v in range(2) for y in range(2)]
                for x in range(2) for u in range(2)]
    rr, cc, d = rank_mod(mixed, 2), rank_mod(opposite, 2), sum(bool(mixed[i][i]) for i in range(4))
    assert 2*rr*cc < d*d
    names = ['check_rectangular_coordinate_support.py', 'check_coordinate_support_rank.py',
             'check_partition_channels.py', 'physical_split_crossings.py', 'weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'one_bit_truth_cases': truth_checks,
              'arbitrary_unequal_supported_cases': len(cases), 'cases': cases,
              'regularized_diagonal_transpose_families': families,
              'mixed_support_negative_control': {'ranks': [rr, cc], 'diagonal': d},
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Exact finite controls. The companion proof gives the all-field, all-size '
                       'unequal-alphabet support/rank-sum theorem. Only support-preserving pieces are covered; '
                       'no all-angle or general crossing obstruction and no stronger factoring score.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'one_bit_truth_cases',
                     'arbitrary_unequal_supported_cases', 'mixed_support_negative_control')}))


if __name__ == '__main__':
    main()
