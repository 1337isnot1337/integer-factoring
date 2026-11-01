#!/usr/bin/env python3
"""Small exact transpose-closed rank-one orbit in odd characteristic."""

import argparse
import json
from math import comb
from time import monotonic

from check_odd_weyl_spectrum import labels_at


def rank(matrix, p):
    pivots = {}
    for raw in matrix:
        row = list(raw)
        for pivot in range(len(row)-1, -1, -1):
            if row[pivot]:
                if pivot in pivots:
                    value = row[pivot]
                    row = [(a-value*b) % p for a, b in zip(row, pivots[pivot])]
                else:
                    inv = pow(row[pivot], -1, p)
                    pivots[pivot] = [(a*inv) % p for a in row]
                    break
    return len(pivots)


def matrix_multiply(A, B, p):
    result = [[0]*len(B[0]) for _ in A]
    for i, row in enumerate(A):
        for j, value in enumerate(row):
            if value:
                result[i] = [(a+value*b) % p for a, b in zip(result[i], B[j])]
    return result


def inverse(A, p):
    n = len(A)
    rows = [list(row)+[int(i == j) for j in range(n)] for i, row in enumerate(A)]
    for i in range(n):
        pivot = next(j for j in range(i, n) if rows[j][i])
        rows[i], rows[pivot] = rows[pivot], rows[i]
        scale = pow(rows[i][i], -1, p)
        rows[i] = [(x*scale) % p for x in rows[i]]
        for j in range(n):
            if j != i and rows[j][i]:
                value = rows[j][i]
                rows[j] = [(x-value*y) % p for x, y in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


def calculate(p, h, profile=False):
    started = monotonic()
    m = h//p
    weights = [1]
    w = 1
    while w <= m:
        weights.append(w)
        w *= p
    labels = labels_at(m, weights)
    ids = {a: i for i, a in enumerate(labels)}
    d = len(labels)
    pairs = [(a, y) for a in range(d) for y in range(d)
             if all((x-v) % 2 == 0 for x, v in zip(labels[a], labels[y]))]
    positions = {pair: i for i, pair in enumerate(pairs)}
    n = len(pairs)
    edges = [[], []]

    def entry(alpha, yy, beta):
        polynomial = [1]+[0]*m
        for a, y, b, weight in zip(alpha, yy, beta, weights):
            following = [0]*(m+1)
            for k in range(min(a, y)+1):
                if a-k > b:
                    continue
                c = comb(y, k)*comb(b, a-k) % p
                if c:
                    for i in range(m-weight*k+1):
                        following[i+weight*k] = (following[i+weight*k]+c*polynomial[i]) % p
            polynomial = following
        return [(k, x) for k, x in enumerate(polynomial) if x]

    for a, alpha in enumerate(labels):
        for y, yy in enumerate(labels):
            for b, beta in enumerate(labels):
                uu = tuple(bb+v-aa for aa, v, bb in zip(alpha, yy, beta))
                u = ids.get(uu)
                if u is None or ((a, y) not in positions and (a, u) not in positions):
                    continue
                polynomial = entry(alpha, yy, beta)
                if not polynomial:
                    continue
                if (a, y) in positions:
                    edges[0].append((positions[a, y], positions[b, u], polynomial))
                if (a, u) in positions:
                    edges[1].append((positions[a, u], positions[b, y], polynomial))
    edges += [[(c, r, value) for r, c, value in family] for family in edges[:]]
    basis, depths, pivots = [], [], {}

    def insert(raw, depth):
        value = list(raw)
        for pivot in range(n-1, -1, -1):
            if value[pivot]:
                if pivot in pivots:
                    c = value[pivot]
                    value = [(a-c*b) % p for a, b in zip(value, pivots[pivot])]
                else:
                    inv = pow(value[pivot], -1, p)
                    value = [(a*inv) % p for a in value]
                    pivots[pivot] = value
                    basis.append(value)
                    depths.append(depth)
                    return

    insert([int(a == y) for a, y in pairs], 0)
    position = 0
    while position < len(basis):
        vector = basis[position]
        for family in edges:
            images = [[0]*n for _ in range(m+1)]
            for r, c, polynomial in family:
                value = vector[c]
                if value:
                    for k, coefficient in polynomial:
                        images[k][r] = (images[k][r]+value*coefficient) % p
            for image in images:
                insert(image, depths[position]+1)
        position += 1
    gram = [[sum(x*y for x, y in zip(a, b)) % p for b in basis] for a in basis]
    gram_rank = rank(gram, p)
    result = {"prime": p, "energy": h, "labels": d, "parity_zero_dimension": n,
              "orbit_dimension": len(basis), "pairing_rank": gram_rank,
              "maximum_word_length": max(depths)}
    if profile and gram_rank == len(basis):
        transpose = [list(column) for column in zip(*basis)]
        projection = matrix_multiply(matrix_multiply(transpose, inverse(gram, p), p), basis, p)
        assert matrix_multiply(projection, projection, p) == projection
        A0 = [[0]*n for _ in range(n)]
        for r, c, value in edges[0]:
            A0[r][c] = sum(coefficient for _, coefficient in value) % p
        projected = matrix_multiply(A0, projection, p)
        assert projected == matrix_multiply(projection, A0, p)
        parities = [sum((x & 1) << j for j, x in enumerate(a)) for a in labels]
        blocks = {}
        for a in range(d):
            for y in range(d):
                blocks.setdefault(parities[a] ^ parities[y], []).append((a, y))
        cut_positions = {pair: (key, i) for key, block in blocks.items()
                         for i, pair in enumerate(block)}

        def column_rank(local, original=False):
            matrices = {key: [[0]*len(block) for _ in block] for key, block in blocks.items()}
            def add(a, y, b, u, value):
                if not value:
                    return
                key, r = cut_positions[a, u]
                key2, c = cut_positions[b, y]
                assert key == key2
                matrices[key][r][c] = (matrices[key][r][c]+value) % p
            if original:
                for a, alpha in enumerate(labels):
                    for y, yy in enumerate(labels):
                        for b, beta in enumerate(labels):
                            uu = tuple(bb+v-aa for aa, v, bb in zip(alpha, yy, beta))
                            u = ids.get(uu)
                            if u is None:
                                continue
                            value = 1
                            for aa, vv in zip(alpha, uu):
                                value = value*comb(aa+vv, aa) % p
                            add(a, y, b, u, value)
            for r, (a, y) in enumerate(pairs):
                for c, (b, u) in enumerate(pairs):
                    add(a, y, b, u, -local[r][c] if original else local[r][c])
            return sum(rank(A, p) for A in matrices.values())

        result["central_projector_cut_ranks"] = [len(basis), column_rank(projection)]
        result["projected_constant_cut_ranks"] = [rank(projected, p), column_rank(projected)]
        result["remaining_constant_column_rank"] = column_rank(projected, True)
    result["seconds"] = round(monotonic()-started, 3)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("energies", type=int, nargs="*", default=[3, 9, 27])
    parser.add_argument("--profile", action="store_true")
    args = parser.parse_args()
    for h in args.energies:
        print(json.dumps(calculate(args.prime, h, args.profile)), flush=True)


if __name__ == "__main__":
    main()
