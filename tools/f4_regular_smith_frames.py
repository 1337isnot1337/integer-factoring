#!/usr/bin/env python3
"""Exact regular Smith frames over F4 to a declared finite precision.

Polynomials are the repository's two binary coefficient bitplanes. The
returned polynomial matrices have invertible constant terms and diagonalize
the full source modulo t**precision. No physical gauge or tensor projection
is implicit in this matrix computation.
"""
from __future__ import annotations

from native_matrix_smith_f4 import MatrixOracle
from reflection_laurent_physical_frames import add, valuation
from weyl_lattice_filter_synthesis import polynomial_multiply


ZERO = (0, 0)
ONE = (1, 0)


def polynomial_product(first, second, precision=None):
    a, b = first
    c, d = second
    ac = polynomial_multiply(a, c)
    bd = polynomial_multiply(b, d)
    cross = polynomial_multiply(a ^ b, c ^ d)
    result = ac ^ bd, cross ^ ac
    if precision is not None:
        mask = (1 << precision) - 1
        result = result[0] & mask, result[1] & mask
    return result


def matrix_product(first, second, precision=None):
    rows, inner, columns = len(first), len(second), len(second[0])
    if any(len(row) != inner for row in first):
        raise ValueError('The polynomial matrix shapes do not compose.')
    result = [[ZERO for _ in range(columns)] for _ in range(rows)]
    for i in range(rows):
        for k in range(inner):
            if first[i][k] == ZERO:
                continue
            for j in range(columns):
                if second[k][j] != ZERO:
                    result[i][j] = add(result[i][j], polynomial_product(first[i][k], second[k][j], precision))
    return result


def regular_frame(source, precision):
    size = len(source)
    if not size or any(len(row) != size for row in source) or precision < 1:
        raise ValueError('Use a nonempty square polynomial matrix and positive precision.')
    if any(len(value) != 2 or any(not isinstance(plane, int) or plane < 0 for plane in value)
           for row in source for value in row):
        raise ValueError('F4 polynomial entries are pairs of nonnegative coefficient bitplanes.')
    mask = (1 << precision) - 1
    matrix = [[(value[0] & mask, value[1] & mask) for value in row] for row in source]
    left = [[ONE if i == j else ZERO for j in range(size)] for i in range(size)]
    right = [[ONE if i == j else ZERO for j in range(size)] for i in range(size)]
    orders = []
    for k in range(size):
        candidates = [(valuation(matrix[i][j]), i, j)
                      for i in range(k, size) for j in range(k, size) if matrix[i][j] != ZERO]
        if not candidates:
            raise ValueError('Precision does not expose every Smith channel; increase it.')
        order, row, column = min(candidates)
        matrix[k], matrix[row] = matrix[row], matrix[k]
        left[k], left[row] = left[row], left[k]
        for values in (matrix, right):
            for line in values:
                line[k], line[column] = line[column], line[k]
        unit = tuple(plane >> order for plane in matrix[k][k])
        assert (unit[0] | unit[1]) & 1
        for i in range(k + 1, size):
            if matrix[i][k] == ZERO:
                continue
            assert valuation(matrix[i][k]) >= order
            quotient = tuple(plane >> order for plane in matrix[i][k])
            for values in (matrix, left):
                values[i] = [add(polynomial_product(unit, a, precision),
                                 polynomial_product(quotient, b, precision))
                             for a, b in zip(values[i], values[k])]
            assert matrix[i][k] == ZERO
        for j in range(k + 1, size):
            if matrix[k][j] == ZERO:
                continue
            assert valuation(matrix[k][j]) >= order
            quotient = tuple(plane >> order for plane in matrix[k][j])
            for values in (matrix, right):
                for line in values:
                    line[j] = add(polynomial_product(unit, line[j], precision),
                                  polynomial_product(quotient, line[k], precision))
            assert matrix[k][j] == ZERO
        orders.append(order)
    assert orders == sorted(orders)
    assert all(matrix[i][j] == ZERO for i in range(size) for j in range(size) if i != j)
    assert matrix_product(matrix_product(left, source, precision), right, precision) == matrix
    return tuple(orders), left, right, matrix


def bounded_frame(source, allowance=8):
    if allowance < 1:
        raise ValueError('Use a positive row-target allowance.')
    size = len(source)
    observed = MatrixOracle(source, size).profile([0] * size, [0] * size)['orders']
    assert all(abs(value - round(value)) < 1e-8 for value in observed)
    orders = tuple(round(value) for value in observed)
    result = regular_frame(source, max(orders) + allowance)
    assert result[0] == orders
    return result
