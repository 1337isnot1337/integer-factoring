#!/usr/bin/env python3
"""Certified polynomial representatives of regular Smith frames to finite order.

Every matrix operation is exact modulo z**precision. The returned left and
right polynomial matrices are actual regular units. Their product with the
untruncated source agrees with the returned diagonal through that precision.
This is sufficient for every row target at most precision.
"""
from __future__ import annotations

from native_matrix_smith_f4 import MatrixOracle
from weyl_lattice_filter_synthesis import polynomial_order, polynomial_multiply


def product_mod(first, second, precision):
    mask = (1 << precision)-1
    return [[_xor(polynomial_multiply(first[i][k], second[k][j])
                  for k in range(len(second))) & mask
             for j in range(len(second[0]))] for i in range(len(first))]


def _xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def regular_smith_frame(source, precision):
    size = len(source)
    if not size or any(len(row) != size for row in source) or precision < 1:
        raise ValueError('Use a nonempty square polynomial matrix and positive precision.')
    if any(not isinstance(value, int) or value < 0 for row in source for value in row):
        raise ValueError('Binary polynomials use nonnegative coefficient bit masks.')
    mask = (1 << precision)-1
    matrix = [[value & mask for value in row] for row in source]
    left = [[int(i == j) for j in range(size)] for i in range(size)]
    right = [[int(i == j) for j in range(size)] for i in range(size)]
    orders = []
    for k in range(size):
        candidates = [(polynomial_order(matrix[i][j]), i, j)
                      for i in range(k, size) for j in range(k, size) if matrix[i][j]]
        if not candidates:
            raise ValueError('Precision does not expose every Smith channel; increase it.')
        order, row, column = min(candidates)
        matrix[k], matrix[row] = matrix[row], matrix[k]
        left[k], left[row] = left[row], left[k]
        for values in (matrix, right):
            for line in values:
                line[k], line[column] = line[column], line[k]
        unit = matrix[k][k] >> order
        assert unit & 1
        for i in range(k+1, size):
            if not matrix[i][k]:
                continue
            assert polynomial_order(matrix[i][k]) >= order
            quotient = matrix[i][k] >> order
            for values in (matrix, left):
                values[i] = [(polynomial_multiply(unit, a) ^ polynomial_multiply(quotient, b)) & mask
                             for a, b in zip(values[i], values[k])]
            assert matrix[i][k] == 0
        for j in range(k+1, size):
            if not matrix[k][j]:
                continue
            assert polynomial_order(matrix[k][j]) >= order
            quotient = matrix[k][j] >> order
            for values in (matrix, right):
                for line in values:
                    line[j] = (polynomial_multiply(unit, line[j]) ^
                               polynomial_multiply(quotient, line[k])) & mask
            assert matrix[k][j] == 0
        orders.append(order)
    assert orders == sorted(orders)
    assert all(not matrix[i][j] for i in range(size) for j in range(size) if i != j)
    assert product_mod(product_mod(left, source, precision), right, precision) == matrix
    return tuple(orders), left, right, matrix


def bounded_smith_frame(source, allowance=8):
    """Choose sufficient precision using the complete independent native profile."""
    size = len(source)
    if allowance < 1:
        raise ValueError('Use a positive row allowance.')
    oracle = MatrixOracle([[(value, 0) for value in row] for row in source], size)
    profile = oracle.profile([0]*size, [0]*size)['orders']
    assert all(abs(value-round(value)) < 1e-8 for value in profile)
    orders = tuple(round(value) for value in profile)
    result = regular_smith_frame(source, max(orders)+allowance)
    assert result[0] == orders
    return result


def build_lattice(values, q, degree, heights, raise_cap=8):
    """Build the frozen affine implementation with a bounded frame constructor.

    This temporary module binding is confined to construction in the calling
    Python process; this helper is intended for single-threaded search workers.
    The original implementation and historical source hashes remain intact.
    """
    import full_alphabet_two_cut_lattices as implementation
    previous = implementation.regular_smith_basis
    try:
        implementation.regular_smith_basis = lambda source: bounded_smith_frame(source, raise_cap)
        return implementation.FullAlphabetTwoCutLattice(values, q, degree, heights, raise_cap)
    finally:
        implementation.regular_smith_basis = previous
