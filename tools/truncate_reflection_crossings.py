#!/usr/bin/env python3
"""Exact low-jet-preserving degree truncations of reflection crossings.

For an even source and an even cutoff d, retain coefficients below d and
replace the coefficient of w^d by the unique endpoint correction.  The
general double-endpoint variant retains coefficients below d-1 and uses
the last two coefficients to enforce value C and derivative zero at one.
"""
from __future__ import annotations

from reflection_laurent_physical_frames import degree, digest


def endpoint(index, q):
    x, y, u, v = ((index // q ** leg) % q for leg in range(4))
    return int(x == v and y == u)


def validate(seed, bound):
    q = seed['q_x']
    if seed['q_u'] != q or len(seed['values']) != q ** 4:
        raise ValueError('Use a complete square crossing tensor.')
    if not isinstance(bound, int) or bound < 2 or bound % 2:
        raise ValueError('Use a positive even degree bound of at least two.')
    for index, pair in enumerate(seed['values']):
        if tuple(value.bit_count() % 2 for value in pair) != (endpoint(index, q), 0):
            raise ArithmeticError('The source fails its strict crossing endpoint.')
    return q


def finish(seed, values, bound, retained, method):
    q = seed['q_x']
    mask = (1 << (retained + 1)) - 1
    odd = sum(1 << power for power in range(1, bound + 1, 2))
    for index, (value, original) in enumerate(zip(values, seed['values'])):
        assert tuple(part & mask for part in value) == tuple(part & mask for part in original)
        assert tuple(part.bit_count() % 2 for part in value) == (endpoint(index, q), 0)
        assert all((part & odd).bit_count() % 2 == 0 for part in value)
    actual = max(map(degree, values))
    assert actual <= bound
    result = dict(seed, values=values, degree=actual, strict_endpoint=True,
                  construction=method, source_degree=seed['degree'],
                  source_tensor_sha256=digest(seed['values']), requested_degree_bound=bound,
                  retained_source_coefficients_through=retained)
    # The truncated tensor need not be the original scalar times a Gaussian.
    result.pop('scalar_w_power', None)
    return result


def truncate_even_crossing(seed, bound=None):
    if bound is None:
        bound = seed['interval'] ** 2
    q = validate(seed, bound)
    source_degree = max(map(degree, seed['values']))
    odd = sum(1 << power for power in range(1, source_degree + 1, 2))
    if any((a | b) & odd for a, b in seed['values']):
        raise ValueError('The even-source formula requires zero odd coefficients.')
    mask = (1 << bound) - 1
    values = []
    for index, pair in enumerate(seed['values']):
        low = tuple(part & mask for part in pair)
        correction = (endpoint(index, q) ^ (low[0].bit_count() % 2), low[1].bit_count() % 2)
        values.append(tuple(part ^ (bit << bound) for part, bit in zip(low, correction)))
    return finish(seed, values, bound, bound - 1, 'even-prefix plus endpoint coefficient')


def truncate_double_endpoint(seed, bound=None):
    if bound is None:
        bound = seed['interval'] ** 2
    q = validate(seed, bound)
    source_degree = max(map(degree, seed['values']))
    all_odd = sum(1 << power for power in range(1, source_degree + 1, 2))
    if any((part & all_odd).bit_count() % 2 for pair in seed['values'] for part in pair):
        raise ValueError('The source must have derivative zero at the endpoint.')
    mask = (1 << (bound - 1)) - 1
    odd = all_odd & mask
    values = []
    for index, pair in enumerate(seed['values']):
        low = tuple(part & mask for part in pair)
        derivative = tuple((part & odd).bit_count() % 2 for part in low)
        correction = (endpoint(index, q) ^ (low[0].bit_count() % 2) ^ derivative[0],
                      (low[1].bit_count() % 2) ^ derivative[1])
        values.append(tuple(part ^ (a << (bound - 1)) ^ (b << bound)
                            for part, a, b in zip(low, derivative, correction)))
    return finish(seed, values, bound, bound - 2, 'prefix plus two endpoint jets')
