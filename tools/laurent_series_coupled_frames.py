#!/usr/bin/env python3
"""Exact complete tensors for square-zero Laurent-series physical frames.

For N(w)=sum_j w^-j G_j with N(w)^2=0, apply I+N on x,u and
I+N^T on y,v.  No intermediate polynomiality is assumed.  One common
shift clears all four legs, then the complete negative part is tested.
"""
from __future__ import annotations

from reflection_laurent_physical_frames import (
    add, bell_normalize, bell_restore, degree, scale, valuation,
)
from square_zero_coupled_reflection_frames import apply_matrix, matrix_terms, product, transpose


def validate_coefficients(coefficients, size):
    cleaned = {}
    for pole, matrix in coefficients.items():
        if not isinstance(pole, int) or pole < 1:
            raise ValueError('Laurent coefficients require positive integer pole labels.')
        if len(matrix) != size or any(len(row) != size or any(c not in range(4) for c in row) for row in matrix):
            raise ValueError('Every coefficient is a square F4 physical-letter matrix.')
        if any(c for row in matrix for c in row):
            cleaned[pole] = [list(row) for row in matrix]
    if not cleaned:
        raise ValueError('Use a nonzero Laurent-series generator.')
    squares = {}
    for first, a in cleaned.items():
        for second, b in cleaned.items():
            value = product(a, b)
            previous = squares.setdefault(first + second, [[0] * size for _ in range(size)])
            squares[first + second] = [[x ^ y for x, y in zip(u, v)] for u, v in zip(previous, value)]
    if any(c for matrix in squares.values() for row in matrix for c in row):
        raise ValueError('The complete matrix Laurent polynomial N must satisfy N squared = zero.')
    return cleaned


def finish(seed, numerator, offset):
    mask = (1 << offset) - 1
    if any((a | b) & mask for a, b in numerator):
        raise ArithmeticError('The complete four-leg Laurent-series frame retains a pole.')
    values = [(a >> offset, b >> offset) for a, b in numerator]
    q = seed['q_x']
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        if (a.bit_count() & 1, b.bit_count() & 1) != (int(x == v and y == u), 0):
            raise ArithmeticError('The Laurent-series action changed the strict endpoint.')
    actual = max(degree(value) for value in values)
    if actual != seed['degree']:
        raise ArithmeticError('A strictly negative Laurent frame changed the original leading degree.')
    return dict(seed, values=values, degree=actual,
                monomial_content=min(valuation(value) for value in values if value != (0, 0)))


def construct(seed, coefficients, coordinates='bell', restore=True, bell_data=None):
    q = seed['q_x']
    if seed['q_u'] != q:
        raise ValueError('The four physical alphabets must have the same size.')
    original = validate_coefficients(coefficients, q)
    bell, factors = bell_normalize(seed) if bell_data is None else bell_data
    if coordinates == 'bell':
        matrices = original
    elif coordinates == 'original':
        matrices = {pole: product(product(factors['inverse_factor'], matrix), factors['first_factor'])
                    for pole, matrix in original.items()}
    else:
        raise ValueError('Choose bell or original physical coordinates.')
    maximum = max(matrices)
    current = {i: value for i, value in enumerate(bell['values']) if value != (0, 0)}
    for leg in range(4):
        updated = {index: (a << maximum, b << maximum) for index, (a, b) in current.items()}
        for pole, matrix in matrices.items():
            action = transpose(matrix) if leg % 2 else matrix
            changed = apply_matrix(current, q, leg, action)
            for index, (a, b) in changed.items():
                value = add(updated.get(index, (0, 0)), (a << (maximum - pole), b << (maximum - pole)))
                if value != (0, 0):
                    updated[index] = value
                elif index in updated:
                    del updated[index]
        current = updated
    numerator = [current.get(index, (0, 0)) for index in range(q ** 4)]
    result = finish(bell, numerator, 4 * maximum)
    content = result['monomial_content']
    if restore:
        result = bell_restore(result, factors)
        result['monomial_content'] = content
    result['laurent_series_frame'] = {'coefficients': original, 'coefficient_coordinates': coordinates,
                                      'bell_coefficients': matrices, 'maximum_pole': maximum,
                                      'restored_coordinates': restore}
    return result


def polynomial_product(first, second):
    result = (0, 0)
    a, b = second
    union = a | b
    while union:
        bit = union & -union
        exponent = bit.bit_length() - 1
        coefficient = ((a >> exponent) & 1) + 2 * ((b >> exponent) & 1)
        x, y = scale(first, coefficient)
        result = result[0] ^ (x << exponent), result[1] ^ (y << exponent)
        union ^= bit
    return result


def from_scalar_terms(seed, terms, scalar_coefficients):
    """Reuse the exact coefficients of (I+fG) on all four legs."""
    coefficients = {pole: c for pole, c in scalar_coefficients.items() if c}
    if not coefficients or any(not isinstance(pole, int) or pole < 1 or c not in range(4)
                               for pole, c in coefficients.items()):
        raise ValueError('Use nonzero F4 coefficients at positive pole orders.')
    maximum = max(coefficients)
    polynomial = (sum((c & 1) << (maximum - pole) for pole, c in coefficients.items()),
                  sum((c >> 1) << (maximum - pole) for pole, c in coefficients.items()))
    powers = [(1, 0)]
    for _ in range(4):
        powers.append(polynomial_product(powers[-1], polynomial))
    numerator = [(0, 0)] * len(seed['values'])
    for order, entries in enumerate(terms):
        shift = (4 - order) * maximum
        for index, value in entries.items():
            a, b = polynomial_product(value, powers[order])
            numerator[index] = add(numerator[index], (a << shift, b << shift))
    return finish(seed, numerator, 4 * maximum)
