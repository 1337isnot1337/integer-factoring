#!/usr/bin/env python3
"""Frobenius-linear consequences of the rank-one minor equations over F4."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from quantum_braid_path import MULTIPLY

FIELD = np.array(MULTIPLY, dtype=np.uint8)


def pack(row):
    return tuple(
        int.from_bytes(np.packbits(row >> p & 1, bitorder="little").tobytes(), "little")
        for p in range(2)
    )


def add_row(basis, low, high):
    """Highest-pivot elimination in the basis 1, alpha (alpha^2=alpha+1)."""
    while low | high:
        pivot = (low | high).bit_length() - 1
        value = (low >> pivot & 1) + 2 * (high >> pivot & 1)
        if pivot not in basis:
            if value == 2:
                low, high = low ^ high, low
            elif value == 3:
                low, high = high, low ^ high
            basis[pivot] = low, high
            return
        a, b = basis[pivot]
        if value == 1:
            low, high = low ^ a, high ^ b
        elif value == 2:
            low, high = low ^ b, high ^ a ^ b
        else:
            low, high = low ^ a ^ b, high ^ a


def nullspace(rows, dimension):
    a = np.array(rows, dtype=np.uint8).reshape(-1, dimension)
    pivots, r = [], 0
    for j in range(dimension):
        nonzero = np.flatnonzero(a[r:, j])
        if not len(nonzero):
            continue
        p = r + nonzero[0]
        a[[r, p]] = a[[p, r]]
        a[r] = FIELD[(0, 1, 3, 2)[a[r, j]], a[r]]
        factors = a[:, j].copy()
        factors[r] = 0
        a ^= FIELD[factors[:, None], a[r, None, :]]
        pivots.append(j)
        r += 1
        if r == len(a):
            break
    result = []
    for j in range(dimension):
        if j in pivots:
            continue
        v = np.zeros(dimension, dtype=np.uint8)
        v[j] = 1
        for i, p in enumerate(pivots):
            v[p] = a[i, j]
        result.append(v)
    return np.array(result, dtype=np.uint8).reshape(-1, dimension)


def minor_equations(basis, row, other, columns):
    """Square monomials first, then c_i*c_j for i<j; no field equations."""
    k = len(basis)
    pairs = list(combinations(range(k), 2))
    ii = np.array([i for i, j in pairs], dtype=int)
    jj = np.array([j for i, j in pairs], dtype=int)
    c = np.array([i for i, j in columns], dtype=int)
    e = np.array([j for i, j in columns], dtype=int)
    a, b = basis[:, row, c].T, basis[:, other, e].T
    x, y = basis[:, row, e].T, basis[:, other, c].T
    square = FIELD[a, b] ^ FIELD[x, y]
    mixed = (
        FIELD[a[:, ii], b[:, jj]]
        ^ FIELD[a[:, jj], b[:, ii]]
        ^ FIELD[x[:, ii], y[:, jj]]
        ^ FIELD[x[:, jj], y[:, ii]]
    )
    return np.concatenate((square, mixed), axis=1)


def square_roots(quadratics, dimension):
    # sqrt(a+b*alpha)=(a+b)+b*alpha. Highest pivot <k means no mixed term.
    return [
        [((a ^ b) >> j & 1) + 2 * (b >> j & 1) for j in range(dimension)]
        for p, (a, b) in sorted(quadratics.items())
        if p < dimension
    ]


def coefficient_basis(tensor, column_cut=False):
    if column_cut:
        tensor = tensor.transpose(0, 1, 3, 2, 4)
    size = tensor.shape[1] ** 2
    matrices = tensor.reshape(len(tensor), size, size)
    chosen, indices, packed = [], [], {}
    for i, a in enumerate(matrices):
        old = len(packed)
        add_row(packed, *pack(a.reshape(-1)))
        if len(packed) > old:
            chosen.append(a)
            indices.append(i)
    return np.array(chosen), indices


def rank_one_envelope(basis, stop_dimension):
    """Necessary linear subspace over EVERY extension, with early termination.

    stop_dimension is a target for the linear envelope, not an assumption that
    it is the actual rank-one locus. Callers must check the terminal space.
    """
    k, size, _ = basis.shape
    columns = list(combinations(range(size), 2))
    quadratics, examined = {}, 0
    # Begin with disjoint adjacent row pairs, then complete the list if needed.
    first = [(i, i + 1) for i in range(0, size - 1, 2)]
    row_pairs = first + [p for p in columns if p not in first]
    for row, other in row_pairs:
        for equation in minor_equations(basis, row, other, columns):
            add_row(quadratics, *pack(equation))
        examined += len(columns)
        linears = square_roots(quadratics, k)
        if len(linears) >= k - stop_dimension:
            break
    linears = square_roots(quadratics, k)
    kernel = nullspace(linears, k)
    assert len(kernel) == k - len(linears)
    return kernel, {
        "dimension": k,
        "minor_equations_used": examined,
        "quadratic_row_rank": len(quadratics),
        "pure_square_constraints": len(linears),
        "linear_envelope_dimension": len(kernel),
        "square_root_constraints": linears,
    }


def combine(vector, basis):
    return np.bitwise_xor.reduce(FIELD[vector[:, None, None], basis], axis=0)
