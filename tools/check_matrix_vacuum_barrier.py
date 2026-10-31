#!/usr/bin/env python3
"""Exact boundary identities and matrix-vacuum rank controls; optional catalogue."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_controlled_permutation_smith import binary_rank
from quantum_pascal_crossing import quantum_binomials
from weighted_smith_oracle import WeightedSmithOracle, polynomial_multiply


def rank(matrix, prime):
    rows = [[v % prime for v in row] for row in matrix]
    pivots = {}
    for row in rows:
        for j, value in enumerate(row):
            if not value:
                continue
            if j in pivots:
                factor = value
                row[:] = [(a-factor*b) % prime for a, b in zip(row, pivots[j])]
            else:
                inverse = pow(value, -1, prime)
                pivots[j] = [(a*inverse) % prime for a in row]
                break
    return len(pivots)


def zero_cut(a, b, prime=2):
    d = len(a)
    return [[sum(a[r-x][k-y]*b[k-u][r-v]
                 for r in range(max(x, v), d) for k in range(max(y, u), d)) % prime
             for y, v in product(range(d), repeat=2)]
            for x, u in product(range(d), repeat=2)]


def polynomial_basis(vacuum):
    d = len(vacuum)
    g = quantum_binomials(d)
    return [[polynomial_multiply(vacuum[r-x][k-y], g[k][y])
             if r >= x and k >= y else 0
             for r, k in product(range(d), repeat=2)]
            for x, y in product(range(d), repeat=2)]


def check():
    rng = Random(819932)
    identities = ranks = boundary_minors = trace_controls = determinants = 0
    for prime, d, sample in product((2, 3, 5), range(2, 6), range(6)):
        a, b = [[[rng.randrange(prime) for _ in range(d)] for _ in range(d)] for _ in range(2)]
        c = zero_cut(a, b, prime)
        rc = rank(c, prime)
        assert 4*rc >= rank(a, prime)*rank(b, prime)
        ranks += 1
        def entry(x, u, y, v):
            return c[d*x+u][d*y+v] if max(x, u, y, v) < d else 0
        for x, u, y, v in product(range(d), repeat=4):
            difference = (entry(x, u, y, v)-entry(x+1, u, y, v+1)
                          -entry(x, u+1, y+1, v)+entry(x+1, u+1, y+1, v+1)) % prime
            assert difference == a[d-1-x][d-1-y]*b[d-1-u][d-1-v] % prime
            identities += 1
        def order(vector):
            return next((i for i, value in enumerate(vector) if value), d)
        minor = [[c[d*x+d-1][d*(d-1)+v] for v in range(d)] for x in range(d)]
        expected = d-max(order([a[r][0] for r in range(d)]), order(b[0]))
        assert rank(minor, prime) == expected <= rc
        boundary_minors += 1
        if d <= 3:
            for x, u, y, v in product(range(d), repeat=4):
                left = [[a[r-x][k-y] if r >= x and k >= y else 0
                         for k in range(d)] for r in range(d)]
                right = [[b[r-u][k-v] if r >= u and k >= v else 0
                          for k in range(d)] for r in range(d)]
                assert sum(left[r][k]*right[k][r] for r, k in product(range(d), repeat=2)) % prime == entry(x, u, y, v)
                trace_controls += 1
    for d in (2, 3):
        for _ in range(5):
            p, q = [[[rng.randrange(1, 8) for _ in range(d)] for _ in range(d)] for _ in range(2)]
            p[0][0], q[0][0] = 3, 7
            first, second = polynomial_basis(p), polynomial_basis(q)
            tensor = [[0]*(d*d) for _ in range(d*d)]
            for i, j in product(range(d*d), repeat=2):
                for r, k in product(range(d), repeat=2):
                    tensor[i][j] ^= polynomial_multiply(first[i][d*r+k], second[j][d*k+r])
            _, flags = WeightedSmithOracle(tensor).profile([0]*(d*d), [0]*(d*d))
            expected = 1
            for _ in range(d*d):
                expected = polynomial_multiply(expected, polynomial_multiply(p[0][0], q[0][0]))
            assert flags[-1].polynomial == expected
            determinants += 1
    return dict(boundary_difference_entries=identities, rank_product_controls=ranks,
                boundary_minor_ranks=boundary_minors, literal_operator_traces=trace_controls,
                full_polynomial_row_determinants=determinants,
                fields=[2, 3, 5], factoring_exponent='1/2; unchanged')


def affine_rank_two_catalogue():
    import numpy as np
    subspaces = sorted({tuple(sorted((a, b, a ^ b)))
                        for a in range(1, 16) for b in range(a+1, 16)})
    assert len(subspaces) == 35
    gl = [bits for bits in range(16)
          if ((bits & 1)*((bits >> 3) & 1)) ^ (((bits >> 1) & 1)*((bits >> 2) & 1))]
    assert len(gl) == 6
    codes = set()
    for left, right, change in product(subspaces, subspaces, gl):
        code = 0
        for i, j in product(range(4), repeat=2):
            bit = 0
            for a, b in product(range(2), repeat=2):
                bit ^= ((left[a] >> i) & 1)*((change >> (2*a+b)) & 1)*((right[b] >> j) & 1)
            code |= bit << (4*i+j)
        codes.add(code)
    assert len(codes) == 7350
    features = [1 | (a << 1) for a in range(8)]
    evaluations = [sum(((c & f).bit_count() % 2) << j for j, f in enumerate(features)) for c in range(16)]
    column_masks = [sum(((255 << y) & 255) << (8*x) for x in range(8)) for y in range(8)]
    def transpose(a):
        result = 0
        while a:
            bit = a & -a
            r, c = divmod(bit.bit_length()-1, 8)
            result |= 1 << (8*c+r)
            a ^= bit
        return result
    histogram, minimum = Counter(), None
    for code in sorted(codes):
        packed = 0
        for x, f in enumerate(features):
            vector = 0
            for i in range(4):
                if f >> i & 1:
                    vector ^= (code >> (4*i)) & 15
            packed |= evaluations[vector] << (8*x)
        assert binary_rank([(packed >> (8*x)) & 255 for x in range(8)]) == 2
        operators = np.array([((packed << (8*x+y)) & ((1 << 64)-1)) & column_masks[y]
                              for x, y in product(range(8), repeat=2)], dtype=np.uint64)
        transposes = np.array([transpose(int(a)) for a in operators], dtype=np.uint64)
        row = (np.bitwise_count(operators[:, None] & transposes[None, :]) & 1).astype(np.uint8)
        column = row.reshape(8, 8, 8, 8).transpose(0, 2, 1, 3).reshape(64, 64)
        rows = list(map(int, np.packbits(column, axis=1, bitorder='little').copy().view('<u8').ravel()))
        r = binary_rank(rows)
        histogram[r] += 1
        if minimum is None or r < minimum['zero_cut_rank']:
            minimum = dict(code=code, packed_matrix=packed, zero_cut_rank=r)
    a = [[(minimum['packed_matrix'] >> (8*x+y)) & 1 for y in range(8)] for x in range(8)]
    literal = zero_cut(a, a)
    assert rank(literal, 2) == minimum['zero_cut_rank'] == 24
    assert binary_rank([sum(v << j for j, v in enumerate(row)) for row in literal], low=True) == 24
    assert (64*24)**2 > 8**7
    return dict(coefficient_matrices=len(codes), rank_histogram=dict(sorted(histogram.items())),
                minimum=minimum, exact_bare_cost_floor=64*24,
                bare_rho_lower_bound='strictly greater than 7/8',
                scope='Every polynomial matrix-vacuum path with P(0)=Q(0)=F C F^T and rank(C)=2; F(s)=(1,s0,s1,s2).')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--affine-catalogue', action='store_true', help='Complete 7,350-case catalogue; requires NumPy 2.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    if args.affine_catalogue:
        result['affine_rank_two'] = affine_rank_two_catalogue()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
