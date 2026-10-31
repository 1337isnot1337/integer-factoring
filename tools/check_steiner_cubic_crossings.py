#!/usr/bin/env python3
"""Exact constant-orbit and full-cut controls for Steiner cubic crossings."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_quantum_pascal_crossings import determinant


def table(q, kind):
    if kind == 'boolean':
        return [[x if x == y else ((x+1) ^ (y+1))-1 for y in range(q)] for x in range(q)]
    if kind == 'affine':
        m = 0
        while 3**m < q:
            m += 1
        assert 3**m == q
        return [[sum((-(x//3**i)-(y//3**i)) % 3 * 3**i for i in range(m))
                 for y in range(q)] for x in range(q)]
    raise ValueError(kind)


def curve(operation):
    q, n = len(operation), len(operation)**2
    assert all(sorted(row) == list(range(q)) for row in operation)
    assert all(operation[x][x] == x and operation[x][y] == operation[y][x]
               and operation[operation[x][y]][y] == x for x, y in product(range(q), repeat=2))
    left = [q*operation[u][v]+v for u, v in product(range(q), repeat=2)]
    right = [q*u+operation[u][v] for u, v in product(range(q), repeat=2)]
    lr, rl = [left[right[i]] for i in range(n)], [right[left[i]] for i in range(n)]
    swap = [q*v+u for u, v in product(range(q), repeat=2)]
    assert [left[right[left[i]]] for i in range(n)] == swap
    matrix = [[0]*n for _ in range(n)]
    for p, value in ((list(range(n)), 3), (right, 10), (lr, 12), (rl, 12), (swap, 8)):
        for j, i in enumerate(p):
            matrix[i][j] ^= value
    assert all(matrix[i][j].bit_count() % 2 == int(i == swap[j])
               for i, j in product(range(n), repeat=2))
    return matrix


def coefficient_rows(matrix, reverse=False):
    return [[sum(((value >> (3-k if reverse else k)) & 1) << j
                 for j, value in enumerate(row)) for row in matrix] for k in range(4)]


def toeplitz_rank(coefficients, length):
    n = len(coefficients[0])
    rows = []
    for i in range(length):
        for r in range(n):
            rows.append(sum(coefficients[i-j][r] << (n*j)
                            for j in range(max(0, i-3), i+1)))
    high = binary_rank(rows)
    assert high == binary_rank(rows, low=True)
    return high


def check():
    small = curve(table(3, 'affine'))
    off = [3*x+y for x, y in product(range(3), repeat=2) if x != y]
    orbit = [[small[i][j] for j in off] for i in off]
    maximum_degrees, minor_count = [], 0
    for k in range(1, 7):
        degrees = []
        for rows, columns in product(combinations(range(6), k), repeat=2):
            value = determinant([[orbit[i][j] for j in columns] for i in rows])
            degrees.append(value.bit_length()-1)
            minor_count += 1
        maximum_degrees.append(max(degrees))
    assert maximum_degrees == [3, 6, 7, 6, 3, 0]
    assert determinant(orbit) == 1
    records, ranks, endpoint_entries = [], 0, 0
    for q, kind in ((3, 'affine'), (7, 'boolean'), (9, 'affine'), (15, 'boolean')):
        operation = table(q, kind)
        matrix = curve(operation)
        column = [[matrix[q*x+y][q*u+v] for y, v in product(range(q), repeat=2)]
                  for x, u in product(range(q), repeat=2)]
        t = q*(q-1)//6
        infinity = [-3]*(2*t)+[-1]*t+[0]*q+[1]*t+[3]*(2*t)
        zero = [0]+[1]*(q-1)+[2]*(4*t)+[3]*(2*t)
        for cut, orders, maximum in ((coefficient_rows(matrix, True), [a+3 for a in infinity], 7),
                                      (coefficient_rows(column), zero, 4)):
            for length in range(1, maximum+1):
                assert toeplitz_rank(cut, length) == sum(max(0, length-a) for a in orders)
                ranks += 1
        # Literal diagonal/off-diagonal block structure, including the coupling.
        for x, u, y, v in product(range(q), repeat=4):
            if x != u and y == v:
                assert column[q*x+u][q*y+v] == 0
            if x == u and y == v:
                assert column[q*x+u][q*y+v] == (3 ^ (2 if x == y else 0))
            if x == u and y != v:
                assert column[q*x+u][q*y+v] == (10 if x == operation[y][v] else 0)
        endpoint_entries += q**4
        records.append(dict(alphabet=q, design=kind, triple_orbits=t,
                            row_infinity_counts={str(a): infinity.count(a) for a in sorted(set(infinity))},
                            column_zero_counts={str(a): zero.count(a) for a in sorted(set(zero))}))
        if q >= 7:
            assert 12*t*t >= q**3
    # q=3: s^6(B(s)-27) has positive coefficients after s=1+h.
    u, v = {3: 2, -3: 2, 1: 1, -1: 1, 0: 3}, {0: 1, -1: 2, -2: 4, -3: 2}
    polynomial = [0]*10
    for a, x in u.items():
        for b, y in v.items():
            polynomial[a+b+6] += x*y
    polynomial[6] -= 27
    shifted = [sum(value*comb(i, j) for i, value in enumerate(polynomial) if i >= j) for j in range(10)]
    assert shifted == [54, 180, 351, 533, 648, 568, 325, 113, 22, 2]
    return dict(complete_constant_orbit_minors=minor_count, maximum_minor_degrees=maximum_degrees,
                full_cut_toeplitz_ranks=ranks, endpoint_and_block_entries=endpoint_entries,
                small_alphabet_shifted_polynomial=shifted, records=records,
                scope='All Steiner quasigroups, the displayed three-involution cubic, characteristic two and its extensions; bare criterion only.',
                bare_rho_lower_bound='3/4', factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
