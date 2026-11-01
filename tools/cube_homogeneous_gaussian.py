#!/usr/bin/env python3
"""Exact Laurent-polynomial arithmetic for a common three-cut cube family.

Polynomials are finite sets of three integer exponent tuples over F2.
Negative powers are used only in the exact Smith change-of-basis matrices.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations

ZERO = frozenset()
ONE = frozenset({(0, 0, 0)})


def monomial(*exponents):
    return frozenset({tuple(exponents)})


def add(first, second):
    return first ^ second


def multiply(first, second):
    result = set()
    for a in first:
        for b in second:
            term = tuple(x + y for x, y in zip(a, b))
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return frozenset(result)


def shift(polynomial, exponents):
    return frozenset(tuple(x + y for x, y in zip(term, exponents)) for term in polynomial)


def valuation(polynomial, axis):
    return min((term[axis] for term in polynomial), default=None)


def identity(order):
    return [[ONE if i == j else ZERO for j in range(order)] for i in range(order)]


def matrix_product(first, second):
    if len(first[0]) != len(second):
        raise ValueError('Incompatible matrix dimensions.')
    result = [[ZERO for _ in second[0]] for _ in first]
    for i, row in enumerate(first):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(second[k]):
                    if other:
                        result[i][j] ^= multiply(value, other)
    return result


def minor(matrix, rows, columns):
    """Determinant; signs disappear in characteristic two."""
    rows, columns = tuple(rows), tuple(columns)
    if len(rows) != len(columns):
        return ZERO
    if not rows:
        return ONE
    value = ZERO
    for j, column in enumerate(columns):
        value ^= multiply(matrix[rows[0]][column],
                          minor(matrix, rows[1:], columns[:j] + columns[j + 1:]))
    return value


def compound(matrix):
    order = len(matrix)
    subsets = [tuple(i for i in range(order) if mask >> i & 1)
               for mask in range(1 << order)]
    return [[minor(matrix, rows, columns) for columns in subsets] for rows in subsets]


def augmentation_coefficients():
    """c_S in the basis product_(i in S)(1+g_i)."""
    coefficients = [ZERO] * 8
    coefficients[1] = monomial(6, 0, 0)
    coefficients[2] = monomial(0, 6, 0)
    coefficients[4] = monomial(0, 0, 6)
    coefficients[3] = monomial(3, 3, 0)
    coefficients[5] = monomial(3, 0, 3)
    coefficients[6] = monomial(0, 3, 3)
    coefficients[7] = monomial(2, 2, 2)
    for coefficient in coefficients[1:]:
        coefficients[0] ^= coefficient
    return coefficients


def kernel():
    coefficients = augmentation_coefficients()
    steps = [ZERO] * 8
    for step in range(8):
        for subset, coefficient in enumerate(coefficients):
            if step & ~subset == 0:
                steps[step] ^= coefficient
    if steps[0]:
        raise ArithmeticError('The Cayley kernel is not alternating.')
    result = [[steps[i ^ j] for j in range(8)] for i in range(8)]
    if any(len(result[i][j]) % 2 != int(i ^ j == 7)
           for i in range(8) for j in range(8)):
        raise ArithmeticError('The one-particle endpoint failed.')
    return result


def gaussian_tensor(source):
    order = len(source)

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return ONE
        if mask.bit_count() & 1:
            return ZERO
        low = mask & -mask
        i = low.bit_length() - 1
        rest, result = mask ^ low, ZERO
        for j in range(i + 1, order):
            if rest >> j & 1:
                result ^= multiply(source[i][j], pfaffian(rest ^ (1 << j)))
        return result

    return tuple(pfaffian(mask) for mask in range(1 << order))


def axis_parts(axis):
    if axis not in range(3):
        raise ValueError('Use one of the three cube axes.')
    return ([i for i in range(8) if not (i >> axis & 1)],
            [i for i in range(8) if i >> axis & 1])


def axis_matrix(tensor, axis):
    first, second = axis_parts(axis)
    result = [[ZERO] * 16 for _ in range(16)]
    for mask, value in enumerate(tensor):
        row = sum((mask >> position & 1) << j for j, position in enumerate(first))
        column = sum((mask >> position & 1) << j for j, position in enumerate(second))
        result[row][column] = value
    return result


def gaussian_cut_factors(source, axis):
    first, second = axis_parts(axis)
    within_first = gaussian_tensor([[source[i][j] for j in first] for i in first])
    within_second = gaussian_tensor([[source[i][j] for j in second] for i in second])
    cross = [[source[i][j] for j in second] for i in first]
    left = [[within_first[row ^ selected] if selected & ~row == 0 else ZERO
             for selected in range(16)] for row in range(16)]
    right = [[within_second[column ^ selected] if selected & ~column == 0 else ZERO
              for column in range(16)] for selected in range(16)]
    return left, compound(cross), right, cross


def augmentation_crossblock(source, axis):
    first, second = axis_parts(axis)
    cross = [[source[i][j] for j in second] for i in first]
    zeta = [[ONE if row & ~column == 0 else ZERO for column in range(4)]
            for row in range(4)]
    return matrix_product(matrix_product(zeta, cross), zeta)


def exact_monomial_smith(matrix, axis):
    """Exact finite Smith operations when a least-valuation monomial is available.

The permitted ring is F2[other variables and their inverses][axis variable].
Every row or column multiplier is checked to have nonnegative axis powers.
The result includes exact left/right transformations and their operation log.
"""
    order = len(matrix)
    value = [row[:] for row in matrix]
    left, right = identity(order), identity(order)
    operations = []
    for pivot in range(order):
        entries = [(valuation(value[i][j], axis), i, j)
                   for i in range(pivot, order) for j in range(pivot, order) if value[i][j]]
        if not entries:
            raise ArithmeticError('The square matrix is singular.')
        lowest = min(item[0] for item in entries)
        candidates = [item for item in entries if item[0] == lowest and len(value[item[1]][item[2]]) == 1]
        if not candidates:
            raise ArithmeticError('A polynomial pivot requires a different exact normal-form algorithm.')
        _, row, column = min(candidates)
        if row != pivot:
            value[row], value[pivot] = value[pivot], value[row]
            left[row], left[pivot] = left[pivot], left[row]
            operations.append(('swap_rows', row, pivot))
        if column != pivot:
            for target in (value, right):
                for target_row in target:
                    target_row[column], target_row[pivot] = target_row[pivot], target_row[column]
            operations.append(('swap_columns', column, pivot))
        exponent = next(iter(value[pivot][pivot]))
        inverse_exponent = tuple(-power for power in exponent)
        for row in range(pivot + 1, order):
            if not value[row][pivot]:
                continue
            quotient = shift(value[row][pivot], inverse_exponent)
            if valuation(quotient, axis) < 0:
                raise ArithmeticError('A row multiplier has a forbidden axis pole.')
            for target in (value, left):
                target[row] = [a ^ multiply(quotient, b) for a, b in zip(target[row], target[pivot])]
            operations.append(('add_row', row, pivot, tuple(sorted(quotient))))
        for column in range(pivot + 1, order):
            if not value[pivot][column]:
                continue
            quotient = shift(value[pivot][column], inverse_exponent)
            if valuation(quotient, axis) < 0:
                raise ArithmeticError('A column multiplier has a forbidden axis pole.')
            for target in (value, right):
                for target_row in target:
                    target_row[column] ^= multiply(quotient, target_row[pivot])
            operations.append(('add_column', column, pivot, tuple(sorted(quotient))))
    if any(value[i][j] for i in range(order) for j in range(order) if i != j):
        raise ArithmeticError('The exact diagonalization left a nonzero off-diagonal entry.')
    if matrix_product(matrix_product(left, matrix), right) != value:
        raise ArithmeticError('The exact Smith transformation did not reconstruct.')
    orders = [valuation(value[i][i], axis) for i in range(order)]
    if orders != sorted(orders):
        raise ArithmeticError('The Smith orders are not ordered.')
    return {'orders': orders, 'diagonal': [value[i][i] for i in range(order)],
            'left': left, 'right': right, 'operations': operations}


def determinantal_orders(matrix, axis):
    order, previous, result, count = len(matrix), 0, [], 0
    for size in range(1, order + 1):
        valuations = []
        for rows in combinations(range(order), size):
            for columns in combinations(range(order), size):
                value = minor(matrix, rows, columns)
                if value:
                    valuations.append(valuation(value, axis))
                count += 1
        current = min(valuations)
        result.append(current - previous)
        previous = current
    return result, count


def cayley_orders(degree, pair_first, pair_second, triple):
    """Exact orders of a+x eps+b eta+c eps eta, with v(a)=degree."""
    e = min(degree, pair_first, pair_second)
    f = min(degree, triple)
    first = min(e, f)
    second = min(2 * e, degree + f) - first
    return [first, second, 2 * degree - second, 2 * degree - first]


def monomial_profile(degree, xy_x, xz_x, yz_y, triple_x, triple_y):
    triple_z = degree - triple_x - triple_y
    if min(degree, xy_x, xz_x, yz_y, triple_x, triple_y, triple_z) < 0:
        raise ValueError('Invalid homogeneous monomial exponents.')
    if max(xy_x, xz_x, yz_y) > degree:
        raise ValueError('A pair exponent exceeds the homogeneous degree.')
    return [cayley_orders(degree, xy_x, xz_x, triple_x),
            cayley_orders(degree, degree - xy_x, yz_y, triple_y),
            cayley_orders(degree, degree - xz_x, degree - yz_y, triple_z)]


def serializable(value):
    if isinstance(value, frozenset):
        return [list(term) for term in sorted(value)]
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    return value
