#!/usr/bin/env python3
"""Bounded regular Smith frames over F4 with the explicit coefficient guard."""
from __future__ import annotations

from weighted_smith_oracle_f4 import polynomial_multiply


def pack(pair):
    a, b = pair
    return sum((((a >> k) & 1) | (((b >> k) & 1) << 1)) << (2*k)
               for k in range(max(a.bit_length(), b.bit_length())))


def planes(value):
    low, high, power = 0, 0, 0
    while value:
        low |= (value & 1) << power
        high |= ((value >> 1) & 1) << power
        value >>= 2
        power += 1
    return low, high


def order(value):
    if not value:
        raise ValueError('The zero polynomial has no finite order.')
    return ((value & -value).bit_length()-1)//2


def matrix_product(first, second, precision=None):
    rows, shared, columns = len(first), len(second), len(second[0])
    assert all(len(row) == shared for row in first)
    assert all(len(row) == columns for row in second)
    mask = (1 << (2*precision))-1 if precision is not None else None
    result = [[0]*columns for _ in range(rows)]
    for i in range(rows):
        for k in range(shared):
            if not first[i][k]:
                continue
            for j in range(columns):
                if second[k][j]:
                    value = polynomial_multiply(first[i][k], second[k][j])
                    result[i][j] ^= value if mask is None else value & mask
    return result


class WindowUnresolved(ArithmeticError):
    pass


def eliminate_window(source, precision):
    n, mask = len(source), (1 << (2*precision))-1
    matrix = [[v & mask for v in row] for row in source]
    left = [[int(i == j) for j in range(n)] for i in range(n)]
    right = [[int(i == j) for j in range(n)] for i in range(n)]
    orders = []
    for k in range(n):
        candidates = [(order(matrix[i][j]), i, j)
                      for i in range(k, n) for j in range(k, n) if matrix[i][j]]
        if not candidates:
            raise WindowUnresolved('A promised F4 pivot is outside this finite window.')
        valuation, row, column = min(candidates)
        matrix[k], matrix[row] = matrix[row], matrix[k]
        left[k], left[row] = left[row], left[k]
        for values in (matrix, right):
            for line in values:
                line[k], line[column] = line[column], line[k]
        unit = matrix[k][k] >> (2*valuation)
        assert unit & 3
        for i in range(k+1, n):
            if not matrix[i][k]:
                continue
            quotient = matrix[i][k] >> (2*valuation)
            for values in (matrix, left):
                values[i] = [(polynomial_multiply(unit, a) ^ polynomial_multiply(quotient, b)) & mask
                             for a, b in zip(values[i], values[k])]
            assert matrix[i][k] == 0
        for j in range(k+1, n):
            if not matrix[k][j]:
                continue
            quotient = matrix[k][j] >> (2*valuation)
            for values in (matrix, right):
                for line in values:
                    line[j] = (polynomial_multiply(unit, line[j])
                               ^ polynomial_multiply(quotient, line[k])) & mask
            assert matrix[k][j] == 0
        orders.append(valuation)
    assert orders == sorted(orders)
    assert all(not matrix[i][j] for i in range(n) for j in range(n) if i != j)
    return tuple(orders), left, right, matrix


def regular_smith_frames(source, raise_cap):
    n = len(source)
    if (not n or any(len(row) != n for row in source)
            or any(type(v) is not int or v < 0 for row in source for v in row)
            or type(raise_cap) is not int or raise_cap < 1):
        raise ValueError('Use a square packed-F4 polynomial matrix and positive allowance.')
    if any(not any(row) for row in source):
        raise ArithmeticError('A zero row prevents full generic rank.')
    degree_bound = sum(max((v.bit_length()-1)//2 for v in row) for row in source)
    limit = 2*degree_bound+raise_cap+1
    precision, attempts = min(64, limit), []
    while True:
        try:
            result = eliminate_window(source, precision)
        except WindowUnresolved:
            result = None
        maximum = max(result[0]) if result is not None else None
        adequate = maximum is not None and precision > 2*maximum+raise_cap
        attempts.append({'precision': precision, 'all_pivots_resolved': result is not None,
                         'maximum_order': maximum, 'adequate': adequate})
        if adequate:
            orders, left, right, diagonal = result
            output_precision = maximum+raise_cap
            mask = (1 << (2*output_precision))-1
            left, right, diagonal = ([[v & mask for v in row] for row in matrix]
                                     for matrix in (left, right, diagonal))
            actual = matrix_product(matrix_product(left, source, output_precision), right, output_precision)
            assert actual == diagonal
            return orders, left, right, diagonal, {
                'field_size': 4, 'internal_precision': precision, 'output_precision': output_precision,
                'source_determinant_degree_bound': degree_bound, 'attempts': attempts}
        if precision == limit:
            raise ArithmeticError('F4 frames did not resolve within the proven determinant bound.')
        precision = min(2*precision, limit)
