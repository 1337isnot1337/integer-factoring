#!/usr/bin/env python3
"""Literal constant F4 basis changes on paired crossing boundaries.

For S=I+a E_ij, apply S on X and S^{-T} on V, or on Y and U.
Since i!=j and the characteristic is two, S^{-1}=S. Each paired
operation preserves the strict endpoint and the polynomial degree.
"""
from __future__ import annotations


def scale_polynomial(pair, coefficient):
    low, high = pair
    if coefficient == 0:
        return 0, 0
    if coefficient == 1:
        return low, high
    if coefficient == 2:
        return high, low ^ high
    if coefficient == 3:
        return low ^ high, low
    raise ValueError('Use a coefficient in F4 encoded by 0,1,2,3.')


def leg_transvection(values, q, leg, target, source, coefficient):
    if not (len(values) == q**4 and 0 <= leg < 4 and 0 <= target < q
            and 0 <= source < q and target != source and coefficient in (1, 2, 3)):
        raise ValueError('Invalid physical-leg transvection.')
    result = list(map(tuple, values))
    stride = q**leg
    for index in range(q**4):
        if (index//stride) % q == target:
            other = index+(source-target)*stride
            first, second = scale_polynomial(values[other], coefficient)
            low, high = result[index]
            result[index] = low ^ first, high ^ second
    return tuple(result)


def paired_transvection(values, q, pair, target, source, coefficient):
    if pair not in (0, 1):
        raise ValueError('Pair zero means X/V; pair one means Y/U.')
    result = leg_transvection(values, q, pair, target, source, coefficient)
    return leg_transvection(result, q, 3-pair, source, target, coefficient)


def transform(values, q, operations):
    result = tuple(map(tuple, values))
    for operation in operations:
        result = paired_transvection(result, q, *operation)
    return result

