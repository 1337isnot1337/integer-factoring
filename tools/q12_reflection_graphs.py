#!/usr/bin/env python3
"""Full sixteen-letter reflection source and arbitrary twelve-letter F4 graphs.

Reuse the original literal tensor and the generic physical split contraction.
No occupation-closure or Gaussian profile substitution is assumed for a graph.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from physical_split_crossings import F4_PRODUCT, check_split, retract
from reflection_dual_rail_gaussian import tensor as reflection_tensor

HEIGHTS = tuple(sum((7, 5, 3, 1)[i] for i in range(4) if word >> i & 1) for word in range(16))
WORDS = tuple(sorted(range(16), key=lambda word: (HEIGHTS[word], word)))
Q, COMPLEMENT = 12, 4
PARAMETERS = 4*Q*COMPLEMENT
FROBENIUS = (0, 1, 3, 2)


def deflate(value):
    result = 0
    while value:
        bit = value & -value
        power = bit.bit_length()-1
        if power % 2:
            raise ArithmeticError('The full source is not even in w.')
        result |= 1 << (power//2)
        value ^= bit
    return result


@lru_cache(None)
def source():
    original = reflection_tensor(8)
    binary = [deflate(low) for low, high in original['values'] if not high]
    if len(binary) != 16**4:
        raise ArithmeticError('The reference source unexpectedly has nonbinary entries.')
    values = []
    for index in range(16**4):
        x, y, u, v = (WORDS[index//16**leg % 16] for leg in range(4))
        value = binary[x+16*(y+16*(u+16*v))]
        if value.bit_count() % 2 != int(x == v and y == u):
            raise ArithmeticError('A full source endpoint failed.')
        values.append((value, 0))
    return {'q_x': 16, 'q_u': 16, 'dimensions': [16]*4, 'values': values,
            'strict_endpoint': True, 'denominator_degree': 0,
            'degree': max(low.bit_length()-1 for low, _ in values)}


def graph_split(coefficients):
    if len(coefficients) != 2*Q*COMPLEMENT or any(c not in range(4) for c in coefficients):
        raise ValueError('A graph requires 96 F4 coordinates.')
    a = [list(coefficients[i*Q:(i+1)*Q]) for i in range(COMPLEMENT)]
    begin = Q*COMPLEMENT
    b = [list(coefficients[begin+i*COMPLEMENT:begin+(i+1)*COMPLEMENT]) for i in range(Q)]
    right = [[int(i == j) for i in range(Q)]+[a[k][j] for k in range(COMPLEMENT)]
             for j in range(Q)]
    left = []
    for i in range(Q):
        row = []
        for j in range(Q):
            value = int(i == j)
            for k in range(COMPLEMENT):
                value ^= F4_PRODUCT[b[i][k]][a[k][j]]
            row.append(value)
        left.append(row+b[i])
    check_split(left, right)
    return left, right


def tensor(coefficients):
    coefficients = tuple(coefficients)
    if len(coefficients) != PARAMETERS:
        raise ValueError('Supply 192 F4 coefficients for the independent paired graphs.')
    lx, px = graph_split(coefficients[:96])
    ly, py = graph_split(coefficients[96:])
    return retract(source(), lx, px, ly, py)


def canonical(coefficients):
    c = tuple(coefficients)
    return min(c, tuple(FROBENIUS[v] for v in c))


def natural():
    h = tuple(float(HEIGHTS[w]) for w in WORDS[1:Q])
    zero = (0.,)*(Q-1)
    return np.asarray(zero+zero+h+h+h+h+zero+zero+(1.,))
