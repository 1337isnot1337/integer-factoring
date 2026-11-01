#!/usr/bin/env python3
"""Literal F4 qutrit projector/Clifford words with a crossing endpoint.

For an order-three Weyl operator U, P=I+U+U^2 and
F_U,c(z)=z I+c(1+z)P. Each F_U,c(1)=I. A three-gate cell gives a
Bell rank-one zero coefficient and a full transpose leading channel.
Degree-two loops z^2 I+d(1+z)P preserve the leading channel. All
polynomial degrees and field planes remain.
"""
from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
import json

import numpy as np

from native_matrix_smith_f4 import MatrixOracle
from regularized_transpose_words import field_rank


@lru_cache(None)
def digits(length):
    values = np.arange(3**length, dtype=np.int64)
    return (values[:, None] // (3**np.arange(length-1, -1, -1, dtype=np.int64))) % 3


@lru_cache(None)
def monomial(a, b):
    states = digits(len(a))
    a, b = np.asarray(a), np.asarray(b)
    powers = 3**np.arange(len(a)-1, -1, -1, dtype=np.int64)
    permutation = ((states+a) % 3) @ powers
    exponent = (2*(a@b)+states@b) % 3
    return permutation, np.asarray((1, 2, 3), dtype=np.uint64)[exponent]


def right_pauli(planes, a, b):
    permutation, phases = monomial(tuple(a), tuple(b))
    low, high = (plane[:, permutation].copy() for plane in planes)
    for coefficient in (2, 3):
        mask = phases == coefficient
        a0, a1 = low[:, mask].copy(), high[:, mask].copy()
        if coefficient == 2:
            low[:, mask], high[:, mask] = a1, a0 ^ a1
        else:
            low[:, mask], high[:, mask] = a0 ^ a1, a0
    return low, high


def factor(planes, a, b, coefficient=3, degree=1):
    first = right_pauli(planes, a, b)
    second = right_pauli(planes, tuple(2*x % 3 for x in a), tuple(2*x % 3 for x in b))
    pl = planes[0] ^ first[0] ^ second[0]
    ph = planes[1] ^ first[1] ^ second[1]
    if coefficient == 3:
        cl, ch = pl ^ ph, pl
    elif coefficient == 2:
        cl, ch = ph, pl ^ ph
    elif coefficient == 1:
        cl, ch = pl, ph
    else:
        raise ValueError('Use omega or omega^2 as the projector coefficient.')
    return [(planes[0] << np.uint64(degree)) ^ cl ^ (cl << np.uint64(1)),
            (planes[1] << np.uint64(degree)) ^ ch ^ (ch << np.uint64(1))]


def nullspace_mod3(matrix, columns):
    rows = [[int(v) % 3 for v in row] for row in matrix]
    pivots, rank = [], 0
    for j in range(columns):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = rows[rank][j]
        rows[rank] = [inverse*v % 3 for v in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][j]:
                c = rows[i][j]
                rows[i] = [(v-c*w) % 3 for v, w in zip(rows[i], rows[rank])]
        pivots.append(j)
        rank += 1
    basis = []
    for free in range(columns):
        if free in pivots:
            continue
        value = [0]*columns
        value[free] = 1
        for i, j in enumerate(pivots):
            value[j] = -rows[i][free] % 3
        basis.append(value)
    return basis


def word(k, loops=()):
    if k not in (1, 2, 3):
        raise ValueError('Use one through three qutrits on each physical leg.')
    prefix, suffix = [], []
    for i in range(k):
        sumx, diffz, sumz = [0]*(2*k), [0]*(2*k), [0]*(2*k)
        sumx[i] = sumx[k+i] = 1
        diffz[i], diffz[k+i] = 1, 2
        sumz[i] = sumz[k+i] = 1
        zero = [0]*(2*k)
        prefix.extend([(tuple(sumx), tuple(zero), 3, 1), (tuple(zero), tuple(diffz), 2, 1)])
        suffix.append((tuple(sumx), tuple(sumz), 3, 1))
    middle = []
    for position in range(len(suffix)+1):
        for loop in loops:
            if loop['position'] == position:
                a, b = tuple(loop['a']), tuple(loop['b'])
                if len(a) != 2*k or len(b) != 2*k or any(v not in (0, 1, 2) for v in (*a, *b)) or not any((*a, *b)):
                    raise ValueError('Supply a nonzero ternary Pauli vector of length 4k.')
                coefficient = loop.get('coefficient', 3)
                if coefficient not in (1, 2, 3):
                    raise ValueError('A loop coefficient must be nonzero in F4.')
                middle.append((a, b, coefficient, 2))
        if position < len(suffix):
            middle.append(suffix[position])
    if any(not 0 <= loop['position'] <= len(suffix) for loop in loops):
        raise ValueError('A loop position is outside the completion word.')
    result = prefix+middle
    if sum(gate[3] for gate in result) > 62:
        raise ValueError('The literal uint64 polynomial degree may not exceed 62.')
    return result


def parity(plane):
    result = plane.copy()
    for shift in (32, 16, 8, 4, 2, 1):
        result ^= result >> np.uint64(shift)
    return result & 1


def construct(k, loops=()):
    raw_word = word(k, loops)
    gates = raw_word
    q, n = 3**k, 3**(2*k)
    center_rows = nullspace_mod3([[(loop['a'][i]-loop['a'][k+i]) % 3 for i in range(k)]
                                 for loop in loops], k)
    for a, _, _, _ in gates:
        assert all(sum(c[i]*(a[i]-a[k+i]) for i in range(k)) % 3 == 0 for c in center_rows)
    planes = [np.eye(n, dtype=np.uint64), np.zeros((n, n), dtype=np.uint64)]
    for a, b, coefficient, power in gates:
        planes = factor(planes, a, b, coefficient, power)
    degree = sum(gate[3] for gate in gates)
    if not np.array_equal(parity(planes[0]), np.eye(n, dtype=np.uint64)) or np.any(parity(planes[1])):
        raise ArithmeticError('A complete Pauli-word endpoint failed.')
    top = ((planes[0] >> np.uint64(degree)) & 1) | (((planes[1] >> np.uint64(degree)) & 1) << np.uint64(1))
    zero = (planes[0] & 1) | ((planes[1] & 1) << np.uint64(1))
    table = np.asarray(((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2)), dtype=np.uint8)
    def rank_one_info(matrix):
        where = np.flatnonzero(matrix)
        if not len(where):
            return {'rank': 0, 'vector_matrix_ranks': [0, 0]}
        i, j = divmod(int(where[0]), n)
        left = matrix[:, j].astype(np.uint8)
        right = table[(0,1,3,2)[int(matrix[i,j])], matrix[i].astype(np.uint8)]
        if not np.array_equal(matrix, table[left[:, None], right[None, :]]):
            raise ArithmeticError('A promised rank-at-most-one coefficient has higher rank.')
        return {'rank': 1, 'vector_matrix_ranks': [field_rank(v.reshape(q,q)) for v in (left,right)]}
    zero_info = rank_one_info(zero)
    top_info = rank_one_info(top.reshape(q,q,q,q).transpose(0,3,2,1).reshape(n,n))
    if top_info['rank'] != 1 or top_info['vector_matrix_ranks'] != [q,q]:
        raise ArithmeticError('The leading opposite coefficient is not a full transpose channel.')
    metadata = {'k': k, 'q': q, 'degree': degree, 'loops': loops,
                'word': gates, 'central_rows': center_rows, 'R_charge_sign': -1,
                'zero_R': zero_info, 'top_C': top_info,
                'zero_R_determinant_order': 2*degree*n//3,
                'endpoint_entries': n*n,
                'planes_sha256': sha256(b''.join(plane.astype('<u8',copy=False).tobytes() for plane in planes)).hexdigest()}
    return planes, metadata


def charge_labels(k, center_rows, sign):
    q = 3**k
    a = np.repeat(digits(k), q, axis=0)
    b = np.tile(digits(k), (q, 1))
    if not center_rows:
        return np.zeros(q*q, dtype=np.int64)
    charges = ((a+sign*b) @ np.asarray(center_rows).T) % 3
    return charges @ (3**np.arange(len(center_rows)-1, -1, -1, dtype=np.int64))


def profiles(planes, metadata, precision=512):
    q, n, k = metadata['q'], metadata['q']**2, metadata['k']
    result, records = [], []
    for infinity in (True, False):
        matrix_planes = ([p.reshape(q, q, q, q).transpose(0, 3, 2, 1).reshape(n, n)
                          for p in planes] if infinity else planes)
        labels = charge_labels(k, metadata['central_rows'], 1 if infinity else -1)
        same = labels[:, None] == labels[None, :]
        if any(np.any(p[~same]) for p in matrix_planes):
            raise ArithmeticError('An off-charge polynomial entry is nonzero.')
        cut_orders, blocks = [], []
        for label in sorted(set(map(int, labels))):
            indices = np.flatnonzero(labels == label)
            matrix = [[tuple(int(p[i, j]) for p in matrix_planes) for j in indices] for i in indices]
            oracle = MatrixOracle(matrix, len(indices), precision=precision)
            orders = oracle.profile(infinity=infinity)['orders']
            if not all(v == int(v) for v in orders):
                raise ArithmeticError('A unit-height order was not integral.')
            cut_orders.extend(map(int, orders))
            blocks.append({'charge': label, 'size': len(indices), 'orders': list(map(int, orders))})
        if len(cut_orders) != n:
            raise ArithmeticError('A complete cut was not resolved.')
        cut_orders.sort()
        result.append(cut_orders)
        records.append({'infinity': infinity, 'blocks': blocks})
    if sum(result[1]) != metadata['zero_R_determinant_order']:
        raise ArithmeticError('The exact determinant order disagrees with the full zero profile.')
    return result, records


def literal_seed(planes, metadata):
    q = metadata['q']
    values = list(zip(*(list(map(int, p.reshape(q, q, q, q).transpose(2, 3, 1, 0).flat)) for p in planes)))
    return {'q_x': q, 'q_u': q, 'dimensions': [q]*4, 'degree': metadata['degree'],
            'values': values, 'strict_endpoint': True, 'denominator_degree': 0}
