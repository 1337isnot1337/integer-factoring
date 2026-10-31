#!/usr/bin/env python3
"""Polynomial crossings from one free height-dressed extreme-core row.

The constructor uses direct Gaussian minors and exact F4 polynomials. It
is intended for small controls; its four bulk alphabets grow exponentially.
The all-size proof is in research/polynomial-frozen-core-crossings.md.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_rectangular_f4_repair import (
    ZERO, HX, HU, add, scale, binary_mul, binary_power, one_particle,
)


def row_components(indices=(0, 1, 2, 3), ambient=4, beta=1):
    """Return the twelve polynomial source tensors before row combination.

    All four retained physical cores are low. Inputs therefore do not
    depend on the mixed input parameter alpha. The source output y-top
    covector is top + beta*both, as in the rectangular construction.
    """
    if beta not in (0, 1, 2, 3):
        raise ValueError('The output parameter must lie in F4.')
    indices = tuple(indices)
    matrix = one_particle(indices, ambient)
    n = len(indices)
    bulk, top = 1 << (n-2), 1 << (n-1)

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i, remaining, answer = bit.bit_length()-1, columns, 0
        while remaining:
            column = remaining & -remaining
            remaining ^= column
            j = column.bit_length()-1
            if matrix[i][j]:
                answer ^= binary_mul(matrix[i][j], minor(rows ^ bit, columns ^ column))
        return answer

    components = [[ZERO]*(bulk**4) for _ in range(12)]
    for source in range(12):
        vcore, ycore = source % 4, source // 4
        vmask = (1, 0, top, top | 1)[vcore]
        ymasks = (((1, 1),), ((0, 1),), ((top, 1), (top | 1, beta)))[ycore]
        height = HX[vcore]+HU[ycore]
        gaussian_factor = binary_power(3, height)
        for x, y, u, v in product(range(bulk), repeat=4):
            xx, uu, vv = 1 | (x << 1), 1 | (u << 1), vmask | (v << 1)
            answer = ZERO
            for ymask, coefficient in ymasks:
                if not coefficient:
                    continue
                yy = ymask | (y << 1)
                polynomial = binary_mul(gaussian_factor, minor(vv | (yy << n), uu | (xx << n)))
                if vv & top and yy & top:
                    # A quartic source always has height two, so its
                    # dressed factor c^2 cancels the original c^-2 pole.
                    assert height == 2
                    polynomial ^= minor((vv ^ top) | ((yy ^ top) << n),
                                        (uu ^ 1) | ((xx ^ 1) << n))
                answer = add(answer, scale((polynomial, 0), coefficient))
            components[source][x+bulk*(y+bulk*(u+bulk*v))] = answer
    heights = [sum(indices[j+1] for j in range(n-2) if mask >> j & 1)
               for mask in range(bulk)]
    return {'indices': list(indices), 'ambient': ambient, 'beta': beta,
            'bulk_size': bulk, 'logical_heights': heights, 'components': components}


def combine_row(data, coefficients, require_crossing=True):
    """Combine source tensors, requiring g[0]+g[10]=1 by default."""
    if len(coefficients) != 12 or any(value not in (0, 1, 2, 3) for value in coefficients):
        raise ValueError('Supply twelve coefficients in F4.')
    endpoint_scalar = coefficients[0] ^ coefficients[10]
    if require_crossing and endpoint_scalar != 1:
        raise ValueError('A crossing row must satisfy g[0]+g[10]=1 in F4.')
    bulk = data['bulk_size']
    values = [ZERO]*(bulk**4)
    for coefficient, component in zip(coefficients, data['components']):
        if coefficient:
            values = [add(current, scale(source, coefficient))
                      for current, source in zip(values, component)]
    return {'field': 4, 'indices': data['indices'], 'ambient': data['ambient'],
            'beta': data['beta'], 'row_coefficients': list(coefficients),
            'bulk_size': bulk, 'dimensions': [bulk]*4,
            'denominator_degree': 0, 'values': values,
            'logical_heights': data['logical_heights'],
            'ramified_endpoint_scalar': endpoint_scalar,
            'strict_endpoint': False,
            'typed_endpoint': 'independent matching X/V and U/Y full-letter types'}


def build_frozen_row(coefficients, indices=(0, 1, 2, 3), ambient=4, beta=1):
    return combine_row(row_components(indices, ambient, beta), coefficients)
