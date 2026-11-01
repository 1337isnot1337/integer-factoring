#!/usr/bin/env python3
"""Polynomial paired boundary maps, using exact F4 polynomial bit planes."""

from __future__ import annotations

import numpy as np

ONE = (1, 0)
ZERO = (0, 0)


def add(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def binary_product(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def multiply(a, b):
    ac = binary_product(a[0], b[0])
    bd = binary_product(a[1], b[1])
    return ac ^ bd, binary_product(a[0] ^ a[1], b[0] ^ b[1]) ^ ac


def evaluate_at_one(value):
    return (value[0].bit_count() & 1) | ((value[1].bit_count() & 1) << 1)


def paired_maps(a, b):
    """Return L=(1,a), P=(1+ab,b)^T, with the polynomial identity LP=1."""
    left, right = (ONE, a), (add(ONE, multiply(a, b)), b)
    assert add(multiply(left[0], right[0]), multiply(left[1], right[1])) == ONE
    return left, right


def multiply_arrays(planes, scalar):
    result = [np.zeros_like(planes[0]), np.zeros_like(planes[1])]
    for plane, polynomial in enumerate(scalar):
        while polynomial:
            bit = polynomial & -polynomial
            degree = bit.bit_length() - 1
            x, y = planes[0] << degree, planes[1] << degree
            if plane:
                result[0] ^= y
                result[1] ^= x ^ y
            else:
                result[0] ^= x
                result[1] ^= y
            polynomial ^= bit
    return result


def compress_split(seed, low, left, right):
    """Apply a polynomial split to the two final coordinates of all legs."""
    ambient = low + 2
    if seed['dimensions'] != [ambient] * 4:
        raise ValueError('Supply a square seed with exactly two final boundary letters.')
    if len(left) != 2 or len(right) != 2 or add(multiply(left[0], right[0]), multiply(left[1], right[1])) != ONE:
        raise ValueError('Supply two polynomial vectors whose pairing is exactly one.')
    dimensions = [ambient] * 4
    planes = [np.asarray([v[p] for v in seed['values']], dtype=object)
              .reshape(dimensions, order='F') for p in (0, 1)]
    for axis, row in enumerate((left, left, right, right)):
        moved = [np.moveaxis(plane, axis, 0) for plane in planes]
        selected = [np.zeros_like(moved[0][0]), np.zeros_like(moved[1][0])]
        for index, scalar in enumerate(row):
            term = multiply_arrays([moved[0][low + index], moved[1][low + index]], scalar)
            selected[0] ^= term[0]
            selected[1] ^= term[1]
        planes = [np.moveaxis(np.stack(list(moved[p][:low]) + [selected[p]]), 0, axis)
                  for p in (0, 1)]
        dimensions[axis] = low + 1
    values = list(zip(*(map(int, plane.ravel(order='F')) for plane in planes)))
    q = low + 1
    for index, value in enumerate(values):
        x, y, u, v = ((index // q**axis) % q for axis in range(4))
        if evaluate_at_one(value) != int(x == v and y == u):
            raise ArithmeticError('A complete polynomial crossing endpoint failed.')
    return {'values': values, 'dimensions': dimensions, 'q_x': q, 'q_u': q,
            'degree': max(part.bit_length() - 1 for value in values for part in value),
            'strict_endpoint': True}


def compress(seed, low, a, b):
    """Use L=(1,a), P=(1+ab,b)^T in the seed's original variable."""
    return compress_split(seed, low, *paired_maps(a, b))


def reorder(seed, order):
    """Relabel all four physical alphabets by the same permutation."""
    q = seed['q_x']
    if seed['dimensions'] != [q] * 4 or sorted(order) != list(range(q)):
        raise ValueError('Supply a permutation of one square physical alphabet.')
    selector = np.ix_(order, order, order, order)
    planes = [np.asarray([value[p] for value in seed['values']], dtype=object)
              .reshape([q] * 4, order='F')[selector] for p in (0, 1)]
    return {**seed, 'values': list(zip(*(map(int, plane.ravel(order='F')) for plane in planes)))}
