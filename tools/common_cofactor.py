#!/usr/bin/env python3
"""Exact projected cofactor polynomials for ternary binary quadratics."""
from __future__ import annotations

from itertools import product

from check_joint_determinant_barrier import determinant,span_witness
from check_multipoint_hermite_crossings import flatten


def matrices(seed):
    result = []
    for cut in (0,1):
        coefficients = [flatten(seed,2,cut,k) for k in range(3)]
        result.append([[tuple(coefficients[k][i][j] for k in range(3))
                        for j in range(4)] for i in range(4)])
    return result


def compound(matrix):
    # Basis of three-vectors ordered by the omitted coordinate.
    return [[determinant([[matrix[a][b] for b in range(4) if b != j]
                         for a in range(4) if a != i],3)
             for j in range(4)] for i in range(4)]


def all_projections(cofactors):
    projective = [a for a in product(range(3),repeat=4)
                  if any(a) and next(x for x in a if x) == 1]
    result = {}
    for a in projective:
        row = [[sum(a[i]*cofactors[i][j][k] for i in range(4)) % 3
                for k in range(7)] for j in range(4)]
        for b in projective:
            p = tuple(sum(b[j]*row[j][k] for j in range(4)) % 3 for k in range(7))
            value = sum(p) % 3
            if value:
                p = tuple(x*value % 3 for x in p)
                result.setdefault(p,([x*value % 3 for x in a],list(b)))
    return result


def find_certificate(seed):
    row,col = map(compound,matrices(seed))
    for i,j in product(range(4),repeat=2):
        vectors = [list(p)+[sum(p) % 3] for p in row[i]]
        vectors += [[-x % 3 for x in p]+[0] for p in col[j]]
        weights = span_witness(vectors,[0]*7+[1],3)
        if weights is not None:
            p = [sum(weights[t]*row[i][t][k] for t in range(4)) % 3 for k in range(7)]
            return {'row_left':[int(t == i) for t in range(4)],
                    'row_right':weights[:4],
                    'column_left':[int(t == j) for t in range(4)],
                    'column_right':weights[4:],'polynomial':p}
    # The two exceptional fixed-row cases require a general left vector.
    first,second = all_projections(row),all_projections(col)
    common = set(first).intersection(second)
    if not common:
        return None
    p = min(common)
    return {'row_left':first[p][0],'row_right':first[p][1],
            'column_left':second[p][0],'column_right':second[p][1],
            'polynomial':list(p)}
