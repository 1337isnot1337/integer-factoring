#!/usr/bin/env python3
"""Finite windows of the existing regular binary Smith frames.

The returned coefficients agree with the untruncated algorithm below the
output precision. Precision is increased until an explicit adequacy test
holds; no missing pivot or coefficient is treated as zero evidence.
"""
from __future__ import annotations

from weyl_lattice_filter_synthesis import matrix_multiply, polynomial_multiply, polynomial_order


class WindowUnresolved(ArithmeticError):
    pass


def eliminate_window(source, precision):
    n, mask = len(source), (1 << precision)-1
    matrix = [[value & mask for value in row] for row in source]
    left = [[int(i == j) for j in range(n)] for i in range(n)]
    right = [[int(i == j) for j in range(n)] for i in range(n)]
    orders = []
    for k in range(n):
        candidates = [(polynomial_order(matrix[i][j]), i, j)
                      for i in range(k, n) for j in range(k, n) if matrix[i][j]]
        if not candidates:
            raise WindowUnresolved('A promised pivot is outside this finite window.')
        order, row, column = min(candidates)
        matrix[k], matrix[row] = matrix[row], matrix[k]
        left[k], left[row] = left[row], left[k]
        for values in (matrix, right):
            for line in values:
                line[k], line[column] = line[column], line[k]
        unit = matrix[k][k] >> order
        assert unit & 1
        for i in range(k+1, n):
            if not matrix[i][k]:
                continue
            quotient = matrix[i][k] >> order
            for values in (matrix, left):
                values[i] = [(polynomial_multiply(unit, a) ^ polynomial_multiply(quotient, b)) & mask
                             for a, b in zip(values[i], values[k])]
            assert matrix[i][k] == 0
        for j in range(k+1, n):
            if not matrix[k][j]:
                continue
            quotient = matrix[k][j] >> order
            for values in (matrix, right):
                for line in values:
                    line[j] = (polynomial_multiply(unit, line[j])
                               ^ polynomial_multiply(quotient, line[k])) & mask
            assert matrix[k][j] == 0
        orders.append(order)
    assert orders == sorted(orders)
    assert all(not matrix[i][j] for i in range(n) for j in range(n) if i != j)
    return tuple(orders), left, right, matrix


def windowed_regular_smith_basis(source, raise_cap):
    n = len(source)
    if (not n or any(len(row) != n for row in source)
            or any(type(v) is not int or v < 0 for row in source for v in row)
            or type(raise_cap) is not int or raise_cap < 1):
        raise ValueError('Use a square binary polynomial matrix and positive integer allowance.')
    if any(not any(row) for row in source):
        raise ArithmeticError('A zero row prevents full generic rank.')
    # Every nonzero determinant term has degree at most this row bound.
    degree_bound = sum(max(v.bit_length()-1 for v in row) for row in source)
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
            mask = (1 << output_precision)-1
            left, right, diagonal = ([[v & mask for v in row] for row in matrix]
                                     for matrix in (left, right, diagonal))
            product = matrix_multiply(matrix_multiply(left, source), right)
            assert [[v & mask for v in row] for row in product] == diagonal
            return orders, left, right, diagonal, {
                'internal_precision': precision, 'output_precision': output_precision,
                'source_determinant_degree_bound': degree_bound, 'attempts': attempts}
        if precision == limit:
            raise ArithmeticError('Full-rank frames did not resolve within the proven determinant bound.')
        precision = min(2*precision, limit)
