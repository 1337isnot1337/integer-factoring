#!/usr/bin/env python3
"""Exact paired-dual polynomial transvections on full physical alphabets.

F4 coefficients use (low, high) binary polynomial masks, with omega²=omega+1.
No Gaussian/parity hypothesis or implicit degree reduction is used.
"""
from __future__ import annotations

from check_schur_gaussian_reduction import polynomial_multiply
from physical_split_crossings import F4_PRODUCT


def multiply(a, b):
    aa = polynomial_multiply(a[0], b[0])
    bb = polynomial_multiply(a[1], b[1])
    ab = polynomial_multiply(a[0], b[1])
    ba = polynomial_multiply(a[1], b[0])
    return aa ^ bb, ab ^ ba ^ bb


def polynomial(terms):
    low = high = 0
    for degree, scalar in terms:
        if not isinstance(degree, int) or degree < 0 or scalar not in range(4):
            raise ValueError('Use nonnegative polynomial degrees and F4 coefficients.')
        if scalar & 1:
            low ^= 1 << degree
        if scalar & 2:
            high ^= 1 << degree
    return low, high


def leg_shear(values, dimensions, leg, target, source, f):
    q = dimensions[leg]
    if target == source or min(target, source) < 0 or max(target, source) >= q:
        raise ValueError('A transvection needs two distinct physical symbols.')
    stride = 1
    for size in dimensions[:leg]:
        stride *= size
    result = list(values)
    for index, original in enumerate(values):
        if index // stride % q == target:
            other = index + (source - target) * stride
            changed = multiply(f, values[other])
            result[index] = original[0] ^ changed[0], original[1] ^ changed[1]
    return result


def endpoint(seed):
    qx, qy = seed['q_x'], seed['q_u']
    for index, (lo, hi) in enumerate(seed['values']):
        x, y = index % qx, index // qx % qy
        u, v = index // (qx * qy) % qy, index // (qx * qy * qy)
        if (lo.bit_count() % 2, hi.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('The complete crossing endpoint failed.')
    return len(seed['values'])


def paired_shear(seed, pair, target, source, f):
    """Apply I+fE_ij and its inverse transpose to x/v or y/u."""
    if pair not in (0, 1):
        raise ValueError('Pair0 is x/v; pair1 is y/u.')
    first, second = ((0, 3), (1, 2))[pair]
    dimensions = list(seed['dimensions'])
    values = leg_shear(seed['values'], dimensions, first, target, source, f)
    values = leg_shear(values, dimensions, second, source, target, f)
    result = dict(seed, values=values,
                  degree=max(part.bit_length() - 1 for value in values for part in value))
    endpoint(result)
    return result


def apply_word(seed, word):
    result = seed
    for gate in word:
        result = paired_shear(result, gate['pair'], gate['target'], gate['source'],
                              polynomial(gate['terms']))
    return result
