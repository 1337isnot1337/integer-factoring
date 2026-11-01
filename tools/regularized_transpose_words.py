#!/usr/bin/env python3
"""Literal non-diagonal transpose-channel words around a rank-one bottleneck.

F_A,k(z)=I+(1+z^k) S_A, S_A(M)=A M^T A^(-T), over F4. Each F_A,k(1)=I.
Polynomial bit planes fit uint64 for degree at most 62. Both field planes
and all four indices are retained; no Gaussian or compound profile is used.
"""
from __future__ import annotations

from hashlib import sha256
import json

import numpy as np

from physical_split_crossings import F4_PRODUCT

INVERSE = (0, 1, 3, 2)


def field_rank(matrix):
    a = [list(map(int, row)) for row in matrix]
    rank = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        c = INVERSE[a[rank][j]]
        a[rank] = [F4_PRODUCT[c][v] for v in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][j]:
                c = a[i][j]
                a[i] = [v ^ F4_PRODUCT[c][w] for v, w in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def rank_one_factors(matrix):
    n = len(matrix)
    position = next(((i, j) for i in range(n) for j in range(n) if matrix[i][j]), None)
    if position is None:
        return None
    i, j = position
    left = [int(matrix[k][j]) for k in range(n)]
    right = [F4_PRODUCT[INVERSE[int(matrix[i][j])]][int(v)] for v in matrix[i]]
    if any(int(matrix[a][b]) != F4_PRODUCT[left[a]][right[b]] for a in range(n) for b in range(n)):
        raise ArithmeticError('A promised rank-at-most-one coefficient has higher rank.')
    return left, right


def scalar(low, high, c):
    if c == 0:
        return np.zeros_like(low), np.zeros_like(high)
    if c == 1:
        return low, high
    if c == 2:
        return high, low ^ high
    if c == 3:
        return low ^ high, low
    raise ValueError('An F4 scalar is required.')


def add_slice(planes, axis, destination, source, coefficient):
    a = [slice(None)]*4
    b = [slice(None)]*4
    a[axis], b[axis] = destination, source
    a, b = tuple(a), tuple(b)
    low, high = scalar(planes[0][b], planes[1][b], coefficient)
    # Some returned field expressions are views, so preserve both before mutation.
    low, high = low.copy(), high.copy()
    planes[0][a] ^= low
    planes[1][a] ^= high


def apply_shear_channel(planes, shears, left=False):
    q = planes[0].shape[0]
    if any(not (0 <= i < q and 0 <= j < q and i != j and c in (1, 2, 3)) for i, j, c in shears):
        raise ValueError('Use nonzero elementary shears with distinct indices.')
    axes = (1, 0, 2, 3) if left else (0, 1, 3, 2)
    transformed = [plane.transpose(axes).copy() for plane in planes]
    # A=E_1...E_k. Left matrix multiplication applies E_k first;
    # right matrix multiplication applies E_1 first.
    for i, j, c in reversed(shears):
        add_slice(transformed, 0 if left else 2, i if left else j, j if left else i, c)
    for i, j, c in shears:
        add_slice(transformed, 1 if left else 3, i if left else j, j if left else i, c)
    return transformed


def apply_diagonal_channel(planes, diagonal):
    q = planes[0].shape[0]
    transformed = [plane.transpose(0, 1, 3, 2).copy() for plane in planes]
    for v in range(q):
        for u in range(q):
            c = F4_PRODUCT[diagonal[u]][INVERSE[diagonal[v]]]
            low, high = scalar(transformed[0][:, :, v, u], transformed[1][:, :, v, u], c)
            transformed[0][:, :, v, u], transformed[1][:, :, v, u] = low.copy(), high.copy()
    return transformed


def regularize(planes, twisted, power=1):
    return [old ^ value ^ (value << np.uint64(power)) for old, value in zip(planes, twisted)]


def bottleneck_diagonals(q):
    if not 2 <= q <= 64:
        raise ValueError('Use 2 through 64 physical symbols.')
    # The code classes are 00 once, 11 once, and 01 repeated q-2 times.
    # Only the first two labels differ in both filters. Distinct codewords
    # are unnecessary: repeated other labels are killed in the same filter.
    codes, bits = [0, 3]+[2]*(q-2), [0, 1]
    diagonals = [[1]*q, [1, 3]+[1]*(q-2), [1, 2]+[3]*(q-2)]
    return diagonals, codes, bits


def construct(q, prepend=(), append=(), base_powers=(1, 1, 1)):
    diagonals, codes, bits = bottleneck_diagonals(q)
    shear_list = lambda gate: gate['shears'] if isinstance(gate, dict) else gate
    gate_power = lambda gate: int(gate.get('power', 1)) if isinstance(gate, dict) else 1
    if len(base_powers) != 3 or any(not isinstance(k, int) or k < 1 for k in base_powers):
        raise ValueError('Supply three positive integer base powers.')
    powers = list(base_powers)+[gate_power(g) for g in (*prepend, *append)]
    degree = sum(powers)
    if degree > 62 or len(powers) % 2 == 0 or any(k < 1 for k in powers):
        raise ValueError('Use an odd number of factors, positive powers and total degree at most 62.')
    low = np.zeros((q, q, q, q), dtype=np.uint64)
    for x in range(q):
        for y in range(q):
            low[x, y, x, y] = 1
    planes = [low, np.zeros_like(low)]
    for diagonal, power in zip(diagonals, base_powers):
        planes = regularize(planes, apply_diagonal_channel(planes, diagonal), power)
    for gate in append:
        planes = regularize(planes, apply_shear_channel(planes, shear_list(gate)), gate_power(gate))
    for gate in prepend:
        planes = regularize(planes, apply_shear_channel(planes, shear_list(gate), left=True), gate_power(gate))
    actual_degree = max(int(v).bit_length()-1 for plane in planes for v in plane.flat)
    if actual_degree != degree:
        raise ArithmeticError('The full leading operator unexpectedly cancelled.')
    # R coefficient at zero, and opposite C coefficient at the leading degree.
    zero = (planes[0] & 1) | ((planes[1] & 1) << np.uint64(1))
    top = ((planes[0] >> np.uint64(degree)) & 1) | (((planes[1] >> np.uint64(degree)) & 1) << np.uint64(1))
    zero_factors = rank_one_factors(zero.reshape(q*q, q*q))
    top_factors = rank_one_factors(top.transpose(0, 3, 2, 1).reshape(q*q, q*q))
    if top_factors is None:
        raise ArithmeticError('The leading transpose channel vanished.')
    factor_ranks = lambda factors: [field_rank(np.asarray(v).reshape(q, q)) for v in factors] if factors else [0, 0]
    values = list(zip(*(list(map(int, plane.transpose(2, 3, 1, 0).flat)) for plane in planes)))
    for index, (a, b) in enumerate(values):
        x, y, u, v = (index//q**leg % q for leg in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('A complete word endpoint failed.')
    seed = {'q_x': q, 'q_u': q, 'dimensions': [q]*4, 'values': values,
            'degree': degree, 'strict_endpoint': True, 'denominator_degree': 0}
    metadata = {'q': q, 'degree': degree, 'diagonals': diagonals, 'codes': codes,
                'filter_bits': bits, 'prepend': prepend, 'append': append,
                'base_powers': base_powers, 'factor_count': len(powers),
                'zero_R_rank': int(zero_factors is not None), 'top_C_rank': 1,
                'zero_R_vector_matrix_ranks': factor_ranks(zero_factors),
                'top_C_vector_matrix_ranks': factor_ranks(top_factors),
                'values_sha256': sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()}
    return seed, metadata
