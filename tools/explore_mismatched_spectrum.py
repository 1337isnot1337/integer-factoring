#!/usr/bin/env python3
"""Exact constant-matrix spectra from mismatched Weyl block dumps."""
import argparse
from collections import Counter
import json
import struct
from time import monotonic

from explore_weyl_central_core import (
    minimal_polynomial, poly_degree, poly_divide, poly_gcd,
    poly_multiply, read_header,
)
from weighted_energy_pencil import rank_two
from witt_exchange_jets import multiply


def poly_remainder(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def squarefree_factors(f):
    """Factor a squarefree polynomial using the exact Berlekamp kernel."""
    n = poly_degree(f)
    if n <= 1:
        return [f]
    pivots = {}
    kernel = []
    for j in range(n):
        value = poly_remainder(1 << (2*j), f) ^ (1 << j)
        witness = 1 << j
        while value:
            k = value.bit_length()-1
            if k not in pivots:
                pivots[k] = (value, witness)
                break
            image, relation = pivots[k]
            value ^= image
            witness ^= relation
        if not value:
            kernel.append(witness)
    factors = [f]
    for v in kernel:
        refined = []
        for g in factors:
            divisor = poly_gcd(g, v)
            if divisor not in (1, g):
                refined.extend([divisor, poly_divide(g, divisor)])
            else:
                refined.append(g)
        factors = refined
    assert len(factors) == len(kernel)
    return factors


def factor_polynomial(f):
    if f == 1:
        return Counter()
    derivative = sum(1 << (i-1) for i in range(1, f.bit_length(), 2) if f >> i & 1)
    if not derivative:
        square_root = sum(1 << (i//2) for i in range(0, f.bit_length(), 2) if f >> i & 1)
        return Counter({p: 2*m for p,m in factor_polynomial(square_root).items()})
    common = poly_gcd(f, derivative)
    squarefree = poly_divide(f, common)
    return Counter(squarefree_factors(squarefree)) + factor_polynomial(common)


def poly_divmod(a, b):
    q = 0
    while a.bit_length() >= b.bit_length():
        shift = a.bit_length()-b.bit_length()
        q ^= 1 << shift
        a ^= b << shift
    return q, a


def poly_inverse(a, modulus):
    r0, r1, s0, s1 = modulus, a, 0, 1
    while r1:
        q, remainder = poly_divmod(r0, r1)
        r0, r1 = r1, remainder
        s0, s1 = s1, s0 ^ poly_multiply(q, s1)
    assert r0 == 1
    return poly_remainder(s0, modulus)


def polynomial_matrix(A, f):
    n = len(A)
    if not f:
        return [0]*n
    if f == 1:
        return [1 << i for i in range(n)]
    result = [0]*n
    for k in range(f.bit_length()-1, -1, -1):
        result = multiply(result, A)
        if f >> k & 1:
            result = [row ^ (1 << i) for i,row in enumerate(result)]
    return result


def load_blocks(path):
    with open(path, "rb") as f:
        h, d, count, weights, labels = read_header(f)
        blocks = []
        for _ in range(count):
            n = struct.unpack("<I", f.read(4))[0]
            pairs = struct.unpack("<%dI" % n, f.read(4*n))
            polynomials = struct.unpack("<%dQ" % (n*n), f.read(8*n*n))
            blocks.append((pairs, polynomials))
        assert not f.read()
    return h, d, weights, labels, blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix")
    parser.add_argument("--output")
    parser.add_argument("--projectors", action="store_true")
    args = parser.parse_args()
    started = monotonic()
    h, d, weights, labels, blocks = load_blocks(args.prefix + ".row.bin")
    histogram = Counter()
    rank_histogram = Counter()
    annihilator = 1
    total_rank = 0
    cache = {}
    for pairs, polynomials in blocks:
        n = len(pairs)
        matrix = tuple(sum((polynomials[i*n+j]&1) << j for j in range(n)) for i in range(n))
        if matrix not in cache:
            cache[matrix] = (minimal_polynomial(matrix), rank_two(matrix))
        polynomial, rank = cache[matrix]
        histogram[hex(polynomial)] += 1
        rank_histogram[(n, rank)] += 1
        annihilator = poly_multiply(poly_divide(annihilator, poly_gcd(annihilator, polynomial)), polynomial)
        total_rank += rank
    result = {
        "energy": h, "labels": d, "weights": weights,
        "blocks": len(blocks), "distinct_block_matrices": len(cache),
        "row_rank_A0": total_rank,
        "minimal_polynomial": hex(annihilator),
        "minimal_degree": poly_degree(annihilator),
        "minimal_polynomial_histogram": dict(sorted(histogram.items())),
        "primary_factors": [[hex(p), poly_degree(p), m] for p,m in sorted(factor_polynomial(annihilator).items())],
        "block_rank_histogram": [[n, r, count] for (n, r), count in sorted(rank_histogram.items())],
        "seconds": round(monotonic()-started, 3),
    }
    if args.projectors:
        primary = []
        projectors = []
        for p, m in sorted(factor_polynomial(annihilator).items()):
            primary_power = 1
            for _ in range(m):
                primary_power = poly_multiply(primary_power, p)
            primary.append(primary_power)
            other = poly_divide(annihilator, primary_power)
            e = poly_remainder(poly_multiply(other, poly_inverse(other, primary_power)), annihilator)
            assert poly_remainder(poly_multiply(e,e), annihilator) == e
            projectors.append(e)
        total = 0
        for e in projectors:
            total ^= e
        assert total == 1
        _, _, _, _, column_blocks = load_blocks(args.prefix + ".col.bin")
        location = {}
        for bi, (pairs, _) in enumerate(column_blocks):
            for j, code in enumerate(pairs):
                location[code] = (bi, j)
        column_matrices = [[[0]*len(pairs) for pairs,_ in column_blocks] for _ in projectors]
        projector_ranks = [0]*len(projectors)
        evaluated = {}
        for pairs, polys in blocks:
            n = len(pairs)
            A = tuple(sum((polys[i*n+j]&1) << j for j in range(n)) for i in range(n))
            if A not in evaluated:
                mp = cache[A][0]
                evaluated[A] = [polynomial_matrix(A, poly_remainder(e, mp)) for e in projectors]
            for pi, matrix in enumerate(evaluated[A]):
                projector_ranks[pi] += rank_two(matrix)
                for i, row in enumerate(matrix):
                    alpha, y = divmod(pairs[i], d)
                    while row:
                        bit = row & -row
                        j = bit.bit_length()-1
                        beta, u = divmod(pairs[j], d)
                        bi, r = location[alpha*d+u]
                        bj, c = location[y*d+beta]
                        assert bi == bj
                        column_matrices[pi][bi][r] ^= 1 << c
                        row ^= bit
        column_ranks = [sum(map(rank_two, matrices)) for matrices in column_matrices]
        assert sum(projector_ranks) == d*d
        result["spectral_projectors"] = [
            {"primary": hex(p), "idempotent": hex(e), "row_rank": r, "column_rank": c}
            for p,e,r,c in zip(primary,projectors,projector_ranks,column_ranks)
        ]
        result["projector_rank_sum"] = sum((r*c)**.5 for r,c in zip(projector_ranks,column_ranks))
        result["finite_target"] = d**(2**.5)
        result["seconds"] = round(monotonic()-started,3)
    if args.output:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    print(json.dumps({k:v for k,v in result.items() if "histogram" not in k}, indent=2))


if __name__ == "__main__":
    main()
