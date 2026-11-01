#!/usr/bin/env python3
"""Pascal kernels and complete Clifford matrices with constant GL(L,F4) flags."""
from __future__ import annotations

from cube_f4_bivariate_smith import F4_MUL, F4_INV, add, multiply, scale

ZERO = (0, 0)


def inverse_flag(matrix):
    n = len(matrix)
    if (not n or any(len(row) != n for row in matrix)
            or any(type(value) is not int or value not in range(4) for row in matrix for value in row)):
        raise ValueError('Supply a nonempty square matrix with F4 coefficients 0,1,2,3.')
    work = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = next((i for i in range(column, n) if work[i][column]), None)
        if pivot is None:
            raise ValueError('The F4 flag must be invertible.')
        work[column], work[pivot] = work[pivot], work[column]
        inverse = F4_INV[work[column][column]]
        work[column] = [F4_MUL[inverse][value] for value in work[column]]
        for i in range(n):
            if i != column and work[i][column]:
                coefficient = work[i][column]
                work[i] = [a ^ F4_MUL[coefficient][b] for a, b in zip(work[i], work[column])]
    return [row[n:] for row in work]


def matrix_product(first, second):
    out = [[ZERO for _ in second[0]] for _ in first]
    for i, row in enumerate(first):
        for k, value in enumerate(row):
            if value == ZERO:
                continue
            for j, other in enumerate(second[k]):
                if other != ZERO:
                    if max(other) <= 1:
                        term = scale(value, other[0] | (other[1] << 1))
                    elif max(value) <= 1:
                        term = scale(other, value[0] | (value[1] << 1))
                    else:
                        term = multiply(value, other)
                    out[i][j] = add(out[i][j], term)
    return out


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def change_flag(seed, constant):
    n = seed['n']
    if len(constant) != n:
        raise ValueError('Supply one F4 flag coordinate per Pascal mode.')
    inverse = inverse_flag(constant)
    original = [[tuple(value) if isinstance(value, tuple) else (value, 0) for value in row]
                for row in seed['matrix']]
    left = [[(value & 1, value >> 1) for value in row] for row in constant]
    right = [[(value & 1, value >> 1) for value in row] for row in inverse]
    first = matrix_product([row[:n] for row in original[:n]], right)
    second = matrix_product(left, [row[n:] for row in original[n:]])
    changed = dict(seed)
    changed.update(matrix=[first[i] + original[i][n:] for i in range(n)]
                   + [original[n + i][:n] + second[i] for i in range(n)],
                   coefficient_field=4, flag_changed=True,
                   constant_flag=[list(row) for row in constant])
    return changed


def inverse_kernel(seed):
    n, outer = seed['n'], seed['outer']
    matrix = seed['matrix']
    pascal = [[(value, 0) for value in row] for row in seed['pascal']]
    first, second = [row[:n] for row in matrix[:n]], [row[n:] for row in matrix[n:]]
    pa, bp = matrix_product(pascal, first), matrix_product(second, pascal)
    scaled = [[(int(i == j) << (2 * outer), 0) for j in range(n)] for i in range(n)]
    return [bp[i] + scaled[i] for i in range(n)] + [pascal[i] + pa[i] for i in range(n)]


def clifford_data(seed, place):
    n, outer, matrix = seed['n'], seed['outer'], seed['matrix']
    pascal = [[(value, 0) for value in row] for row in seed['pascal']]
    if place == 0:
        first, second = [row[:n] for row in matrix[:n]], [row[n:] for row in matrix[n:]]
        identity = [[(int(i == j), 0) for j in range(n)] for i in range(n)]
        zero = [[ZERO] * n for _ in range(n)]
        pa, bp = matrix_product(pascal, first), matrix_product(second, pascal)
        blocks = [[zero, identity, transpose(pa), zero],
                  [pascal, zero, zero, first],
                  [second, zero, zero, identity],
                  [zero, transpose(bp), transpose(pascal), zero]]
        shifts = [[0, -2 * outer, -2 * outer, 0], [0, 0, 0, 0],
                  [0, 0, 0, 0], [0, -2 * outer, -2 * outer, 0]]
        values = [[blocks[bi][bj][i][j] for bj in range(4) for j in range(n)]
                  for bi in range(4) for i in range(n)]
        shifts = [[shifts[bi][bj] for bj in range(4) for _ in range(n)]
                  for bi in range(4) for _ in range(n)]
        rows = list(range(n)) + list(range(2 * n, 3 * n))
        columns = list(range(n, 2 * n)) + list(range(3 * n, 4 * n))
    elif place == 1:
        transposed, inverse = transpose(matrix), inverse_kernel(seed)
        zero = [[ZERO] * (2 * n) for _ in range(2 * n)]
        values = [transposed[i] + zero[i] for i in range(2 * n)]
        values += [zero[i] + inverse[i] for i in range(2 * n)]
        shifts = [[0] * (4 * n) for _ in range(4 * n)]
        rows, columns = list(range(2 * n)), list(range(2 * n, 4 * n))
    else:
        raise ValueError('Choose place zero or one.')
    return values, shifts, rows, columns
