#!/usr/bin/env python3
"""Exact finite Boolean substitution operators for Cauchy packets.

The module uses integer residues and exact rational kernel evaluation.
Its dense routines are finite controls, not the polynomial-space inverse
entry algorithm proved using uniform determinant circuits in the note.
"""
from __future__ import annotations

from fractions import Fraction

from cauchy_packet_polynomials import cauchy_polynomials, evaluate


def boundary_rows(length, center):
    return [x for x in range(length) if not 0 <= x + center < length]


def residue(value, modulus):
    value = Fraction(value)
    if not value.denominator & 1:
        raise ValueError('Expected a two-adically integral rational.')
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def kernel_matrices(length, center, point, precision):
    """Return W(t+2L*point), and its boundary rows divided by 2*point.

    At point zero, evaluate the removable quotient by differentiating
    the globally constructed polynomial. No numerical limit is taken.
    """
    modulus = 1 << precision
    kernels = cauchy_polynomials(length)
    boundary = set(boundary_rows(length, center))
    original = []
    normalized = []
    for x in range(length):
        row = []
        scaled = []
        for y in range(length):
            polynomial = kernels[y - x]
            value = evaluate(polynomial, center + 2 * length * point)
            if x not in boundary:
                quotient = value
            elif point:
                quotient = value / (2 * point)
            else:
                derivative = [degree * coefficient for degree, coefficient in enumerate(polynomial)][1:]
                quotient = length * evaluate(derivative, center)
            row.append(residue(value, modulus))
            scaled.append(residue(quotient, modulus))
        original.append(row)
        normalized.append(scaled)
    return original, normalized


def monomial_operator(matrix, modulus):
    """Rows are input squarefree monomials; columns are output monomials."""
    length = len(matrix)
    count = 1 << length
    images = [[sum(matrix[x][y] for y in range(length) if state >> y & 1) % modulus
               for state in range(count)] for x in range(length)]
    values = [[1] * count]
    for subset in range(1, count):
        bit = (subset & -subset).bit_length() - 1
        previous = values[subset & (subset - 1)]
        values.append([a * b % modulus for a, b in zip(previous, images[bit])])
    # Finite differences turn values on the Boolean cube into the unique
    # multilinear coefficients. This includes all repeated-column terms.
    for bit in range(length):
        for state in range(count):
            if state >> bit & 1:
                previous = state ^ (1 << bit)
                for row in values:
                    row[state] = (row[state] - row[previous]) % modulus
    return values


def value_operator(matrix, modulus):
    """Literal Bernstein/Boolean-value operator, constructed independently."""
    length = len(matrix)
    count = 1 << length
    result = [[0] * count for _ in range(count)]
    for column in range(count):
        images = [sum(matrix[x][y] for y in range(length) if column >> y & 1) % modulus
                  for x in range(length)]
        for row in range(count):
            value = 1
            for x, image in enumerate(images):
                value = value * (image if row >> x & 1 else 1 - image) % modulus
            result[row][column] = value
    return result


def monomial_to_value_operator(matrix, modulus):
    """Conjugate by the integral Boolean zeta matrix."""
    result = [row[:] for row in matrix]
    count = len(result)
    length = count.bit_length() - 1
    for bit in range(length):
        for column in range(count):
            if column >> bit & 1:
                for row in result:
                    row[column] = (row[column] + row[column ^ (1 << bit)]) % modulus
    for bit in range(length):
        for row in range(count):
            if not row >> bit & 1:
                other = result[row | (1 << bit)]
                result[row] = [(a - b) % modulus for a, b in zip(result[row], other)]
    return result


def unit_inverse(matrix, modulus):
    """Dense exact inversion when the determinant is odd."""
    count = len(matrix)
    left = [row[:] for row in matrix]
    right = [[int(i == j) for j in range(count)] for i in range(count)]
    for column in range(count):
        pivot = next((row for row in range(column, count) if left[row][column] & 1), None)
        if pivot is None:
            raise ArithmeticError('The matrix is not invertible modulo two.')
        left[column], left[pivot] = left[pivot], left[column]
        right[column], right[pivot] = right[pivot], right[column]
        inverse = pow(left[column][column], -1, modulus)
        left[column] = [value * inverse % modulus for value in left[column]]
        right[column] = [value * inverse % modulus for value in right[column]]
        for row in range(count):
            if row == column:
                continue
            factor = left[row][column]
            if factor:
                left[row] = [(a - factor * b) % modulus for a, b in zip(left[row], left[column])]
                right[row] = [(a - factor * b) % modulus for a, b in zip(right[row], right[column])]
    assert left == [[int(i == j) for j in range(count)] for i in range(count)]
    return right


def multiply_matrices(left, right, modulus):
    count = len(left)
    result = [[0] * count for _ in range(count)]
    for row in range(count):
        for inner, value in enumerate(left[row]):
            if value:
                for column, factor in enumerate(right[inner]):
                    result[row][column] += value * factor
        result[row] = [value % modulus for value in result[row]]
    return result


def smith_orders(matrix, precision):
    """Independent minimum-valuation elimination for small finite controls."""
    modulus = 1 << precision
    work = [row[:] for row in matrix]
    count = len(work)
    result = []
    for index in range(count):
        best = None
        for row in range(index, count):
            for column in range(index, count):
                value = work[row][column]
                if value:
                    order = (value & -value).bit_length() - 1
                    if best is None or order < best[0]:
                        best = order, row, column
        if best is None:
            return result + [None] * (count - index)
        order, row, column = best
        work[index], work[row] = work[row], work[index]
        for values in work:
            values[index], values[column] = values[column], values[index]
        inverse = pow(work[index][index] >> order, -1, modulus)
        work[index] = [value * inverse % modulus for value in work[index]]
        for row in range(index + 1, count):
            factor = work[row][index] >> order
            if factor:
                work[row][index:] = [(a - factor * b) % modulus
                                     for a, b in zip(work[row][index:], work[index][index:])]
        work[index][index + 1:] = [0] * (count - index - 1)
        result.append(order)
    return result


def rotate_subset(subset, length, center):
    return sum(1 << ((position + center) % length)
               for position in range(length) if subset >> position & 1)


def normalized_inverse_mod_two(length, center):
    """Value-basis rows of the normalized inverse, stored as bit vectors."""
    count = 1 << length
    boundary = sum(1 << x for x in boundary_rows(length, center))
    rows = [0] * count
    for subset in range(count):
        if subset & boundary == boundary:
            rows[rotate_subset(subset, length, center)] = 1 << subset
    # Right Boolean zeta transform, then left upper Mobius transform.
    for bit in range(length):
        step = 1 << bit
        lower_mask = sum(1 << state for state in range(count) if not state >> bit & 1)
        rows = [row ^ ((row & lower_mask) << step) for row in rows]
    for bit in range(length):
        for state in range(count):
            if not state >> bit & 1:
                rows[state] ^= rows[state | (1 << bit)]
    return rows


def binary_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def inverse_group_rank(length, center):
    absolute = abs(center)
    if not 0 < absolute < length:
        raise ValueError('A retained coordinate requires 0<|t|<L.')
    coordinate = (absolute & -absolute).bit_length() - 1
    first = [position for position in range(length) if not position >> coordinate & 1]
    second = [position for position in range(length) if position >> coordinate & 1]
    assert len(first) == len(second) == length // 2
    count = 1 << length
    first_labels = [sum(((state >> position) & 1) << i for i, position in enumerate(first))
                    for state in range(count)]
    second_labels = [sum(((state >> position) & 1) << i for i, position in enumerate(second))
                     for state in range(count)]
    flattened = [0] * count
    for row, values in enumerate(normalized_inverse_mod_two(length, center)):
        while values:
            bit = values & -values
            column = bit.bit_length() - 1
            values ^= bit
            first_index = (first_labels[row] << len(first)) | first_labels[column]
            second_index = (second_labels[row] << len(second)) | second_labels[column]
            flattened[first_index] ^= 1 << second_index
    return binary_rank(flattened)
