#!/usr/bin/env python3
"""Closed scalar blocks of an endpoint-preserving two-state polynomial kernel.

At every physical bit tuple a, K_a=P_a I+(1+z)[[0,B_a],[C_a,0]].
The scalar e0^T product K_a e0 has the exact crossing endpoint. A fixed
block may be restricted to any common physical word code before taking
ordinary tensor powers. All polynomials below are exact binary bitmasks.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product


def multiply(a, b):
    answer = 0
    while b:
        bit = b & -b
        answer ^= a << (bit.bit_length()-1)
        b ^= bit
    return answer


def local_entry(x, y, u, v, mask_b, mask_c):
    index = x + 2*y + 4*u + 8*v
    identity = int(x == u and y == v)
    crossing = int(x == v and y == u)
    p = identity | ((identity ^ crossing) << 1)
    return ((p, 3*((mask_b >> index) & 1)),
            (3*((mask_c >> index) & 1), p))


def tensor(length, mask_b, mask_c, labels=None):
    if length < 1 or not 0 <= mask_b < 65536 or not 0 <= mask_c < 65536:
        raise ValueError('positive block length and two 16-bit masks required')
    labels = tuple(range(1 << length)) if labels is None else tuple(labels)
    if len(labels) < 2 or len(set(labels)) != len(labels) or any(
            not 0 <= a < 1 << length for a in labels):
        raise ValueError('physical code must have distinct valid binary words')
    q = len(labels)
    local = [local_entry(i & 1, (i >> 1) & 1, (i >> 2) & 1,
                         (i >> 3) & 1, mask_b, mask_c) for i in range(16)]
    words = [tuple((x >> j) & 1 for j in range(length)) for x in labels]
    values = []
    for v, u, y, x in product(range(q), repeat=4):
        left, right = 1, 0
        for j in range(length):
            index = words[x][j] + 2*words[y][j] + 4*words[u][j] + 8*words[v][j]
            ((p, b), (c, _)) = local[index]
            left, right = multiply(left, p) ^ multiply(right, c), \
                          multiply(left, b) ^ multiply(right, p)
        assert left.bit_length() <= length+1
        assert (left.bit_count() & 1) == int(x == v and y == u)
        values.append((left, 0))
    return {'dimensions': [q]*4, 'values': values, 'strict_endpoint': True,
            'degree': length, 'block_length': length, 'physical_words': labels,
            'mask_b': mask_b, 'mask_c': mask_c}


def product_masks():
    """All nonzero fully separable bit tensors, plus the pairing span."""
    result = set()
    for factors in product((1, 2, 3), repeat=4):
        mask = 0
        for i in range(16):
            coefficient = all(factors[leg] >> ((i >> leg) & 1) & 1 for leg in range(4))
            mask |= int(coefficient) << i
        result.add(mask)
    pairings = []
    for kind in range(3):
        mask = 0
        for i in range(16):
            x, y, u, v = ((i >> k) & 1 for k in range(4))
            condition = ((x == u and y == v), (x == v and y == u),
                         (x == y and u == v))[kind]
            mask |= int(condition) << i
        pairings.append(mask)
    for bits in range(1, 8):
        value = 0
        for i, p in enumerate(pairings):
            if bits >> i & 1:
                value ^= p
        if value:
            result.add(value)
    return sorted(result)


def flatten(values, q, cut):
    out = [[0]*(q*q) for _ in range(q*q)]
    for i, (a, b) in enumerate(values):
        assert b == 0
        x, y, u, v = ((i//q**j) % q for j in range(4))
        row, col = (x+q*y, u+q*v) if cut == 0 else (x+q*u, y+q*v)
        out[row][col] = a
    return out


@lru_cache(None)
def binary_tuples(length):
    return tuple(product(range(2), repeat=length))
