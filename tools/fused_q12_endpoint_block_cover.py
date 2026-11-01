#!/usr/bin/env python3
"""Generate exact all-place block covers for q12 fused F16 endpoint seeds."""
from itertools import product
from collections import Counter
from hashlib import sha256
import json

from compressed_fused_six_vertex import construct_compressed
from f16_identity_operator_profiles import MUL, INV, add, multiply, inverse_matrix
from fused_six_vertex_crossings import pdiv, pscale, power
from check_fused_multipoint_profiles import field_rank
from check_fused_closed_profiles import canonical, irreducible_factors

Q, M = 12, 11
PLACES = tuple(range(16))+(None,)
CHARGES = tuple(tuple(Q*a+b for a, b in product(range(Q), repeat=2) if a+b == c)
                for c in range(2*Q-1))


def digest(matrix):
    return sha256(json.dumps(matrix, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def evaluate(polynomial, point, degree):
    if point is None:
        return polynomial.get(degree, 0)
    result = 0
    for k, c in polynomial.items():
        result ^= MUL[c][power(point, k)]
    return result


def normalise_at(matrix, point, degree):
    n = Q*Q
    endpoint = [[evaluate(p, point, degree) for p in line] for line in matrix]
    output = [[{} for _ in range(n)] for _ in range(n)]
    for indices in CHARGES:
        inverse = inverse_matrix([[endpoint[i][j] for j in indices] for i in indices])
        for a, i in enumerate(indices):
            for b, j in enumerate(indices):
                for c, k in enumerate(indices):
                    if inverse[c][b]:
                        output[i][j] = add(output[i][j], matrix[i][k], inverse[c][b])
    return output, endpoint


def block_entries(matrix, difference):
    return [((x, y, u, v), matrix[Q*y+x][Q*u+v])
            for x, u in product(range(Q), repeat=2) if u-x == difference
            for y, v in product(range(Q), repeat=2) if y-v == difference]


def extended_gcd(a, b):
    x0, x1, y0, y1 = {0: 1}, {}, {}, {0: 1}
    while b:
        quotient, remainder = pdiv(a, b)
        a, b = b, remainder
        x0, x1 = x1, add(x0, multiply(quotient, x1))
        y0, y1 = y1, add(y0, multiply(quotient, y1))
    scalar = INV[a[max(a)]]
    return pscale(a, scalar), pscale(x0, scalar), pscale(y0, scalar)


def bezout(entries):
    current, chosen, coefficients = {}, [], []
    for coordinates, polynomial in entries:
        if not polynomial:
            continue
        new, left, right = extended_gcd(current, polynomial)
        coefficients = [multiply(left, c) for c in coefficients]+[right]
        chosen.append(coordinates)
        current = new
        if current == {0: 1}:
            break
    if not current:
        raise ArithmeticError('A charge block vanished identically.')
    return dict(polynomial=current,
                terms=[dict(coordinates=coordinates, coefficient=c)
                       for coordinates, c in zip(chosen, coefficients) if c])


def pivot_minor(matrix):
    work, labels, columns = [row[:] for row in matrix], list(range(len(matrix))), []
    rank = 0
    for j in range(len(matrix)):
        pivot = next((i for i in range(rank, len(matrix)) if work[i][j]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        labels[rank], labels[pivot] = labels[pivot], labels[rank]
        inverse = INV[work[rank][j]]
        work[rank] = [MUL[inverse][v] for v in work[rank]]
        for i in range(rank+1, len(matrix)):
            scalar = work[i][j]
            if scalar:
                work[i] = [a ^ MUL[scalar][b] for a, b in zip(work[i], work[rank])]
        columns.append(j)
        rank += 1
        if rank == len(matrix):
            break
    return dict(rows=labels[:rank], columns=columns)


def build_case(theta):
    source = construct_compressed(M, theta, normalize_joint=False)
    matrix, degree = source['matrix'], source['degree']
    records, regular = [], []
    for point in PLACES:
        endpoint = [[evaluate(p, point, degree) for p in row] for row in matrix]
        full = all(field_rank([[endpoint[i][j] for j in indices] for i in indices]) == len(indices)
                   for indices in CHARGES)
        regular.append(full)
        if not full:
            continue
        normalised, _ = normalise_at(matrix, point, degree)
        blocks, support, infinity, covered = [], Counter(), [], 0
        for difference in range(1-Q, Q):
            entries = block_entries(normalised, difference)
            identity = bezout(entries)
            factors = irreducible_factors(identity['polynomial'])
            for factor in {canonical(p) for p in factors}:
                support[factor] += 1
            blocks.append(dict(difference=difference, identity=identity,
                               factors=sorted([list(canonical(p)) for p in factors])))
            if covered < Q:
                rows = [(x, u) for x, u in product(range(Q), repeat=2) if u-x == difference]
                cols = [(y, v) for y, v in product(range(Q), repeat=2) if y-v == difference]
                leading = [[normalised[Q*y+x][Q*u+v].get(degree, 0) for y, v in cols] for x, u in rows]
                witness = pivot_minor(leading)
                if witness['rows']:
                    infinity.append(dict(difference=difference, **witness))
                    covered += len(witness['rows'])
        records.append(dict(endpoint=point, normalised_sha256=digest(normalised), blocks=blocks,
                            finite_block_lower=2*Q-1-max(support.values(), default=0),
                            infinity_minors=infinity, infinity_rank_lower=covered))
    return dict(theta=theta, degree=degree, source_sha256=digest(matrix), regular_fibres=regular,
                rectangular_symmetry_checks=source['rectangular_symmetry_checks'],
                output_symmetry_checks=source['output_symmetry_checks'], endpoints=records)
