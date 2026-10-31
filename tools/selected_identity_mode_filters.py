#!/usr/bin/env python3
"""Canonical independent identity pairs for any indices at prime-power length."""
from __future__ import annotations

from functools import lru_cache
from math import factorial

from coupled_identity_mode_filters import algebra_multiply
from prime_state_filter_profiles import validate


@lru_cache(maxsize=64)
def geometry_data(prime, n, selected):
    validate(prime, n, 1, (0,)*prime)
    selected = tuple(selected)
    remainder = n
    while remainder > 1 and remainder % prime == 0:
        remainder //= prime
    if remainder != 1:
        raise ValueError('Arbitrary identity selections require a prime-power interval length.')
    if not selected or selected != tuple(sorted(set(selected))) or min(selected) < 0 or max(selected) >= n:
        raise ValueError('Use distinct increasing selected indices in the interval0..n−1.')
    modes = len(selected)
    identity = tuple(tuple({0: 1} if i == j else {} for j in range(modes)) for i in range(modes))
    # (1-tau)^(n-1)=(1-tau^n)/(1-tau), by Frobenius.
    unit = {degree: 1 for degree in range(n)}
    reflected = tuple(tuple(dict(unit) if i == j else {} for j in range(modes)) for i in range(modes))
    return {'prime': prime, 'n': n, 'modes': modes,
            'geometry': 'selected:'+','.join(map(str, selected)),
            'selected': selected, 'retained': tuple(n-1-j for j in selected),
            'active_x': identity, 'active_v': reflected, 'gaussian_v_unit': {0: 1}}


def active_gaussian(prime, n, ramification, modes, geometry):
    selected = tuple(map(int, geometry.removeprefix('selected:').split(',')))
    data = geometry_data(prime, n, selected)
    if modes != len(selected):
        raise ValueError('The selected-index count must equal the number of modes.')
    zero = (0,)*(4*modes)
    result = {zero: {0: 1}}
    for j, index in enumerate(selected):
        for left, right, order in ((modes+j, 2*modes+j, 2*ramification*index),
                                   (j, 3*modes+j, 2*ramification*(n-index))):
            factor = {}
            for digit in range(prime):
                powers = list(zero)
                powers[left] = powers[right] = digit
                factor[tuple(powers)] = {order*digit: pow(factorial(digit), -1, prime)}
            result = algebra_multiply(result, factor, prime)
    return result
