#!/usr/bin/env python3
"""Ordered two-dimensional physical gauges of the three-letter Steiner cubic."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_steiner_cubic_crossings import curve, table
from weighted_smith_oracle import polynomial_multiply

IDENTITY = 9
GL2 = (6, 7, 9, 11, 13, 14)
FIXED_FIRST = (9, 11, 13)
BOUNDARIES = tuple((e, f) for e in range(1, 4) for f in range(1, 4)
                   if (e & f).bit_count() % 2)


def multiply(a, b):
    result = 0
    for i, j in product(range(2), repeat=2):
        bit = 0
        for k in range(2):
            bit ^= ((a >> (2*i+k)) & 1)*((b >> (2*k+j)) & 1)
        result |= bit << (2*i+j)
    return result


MULTIPLICATION = tuple(tuple(multiply(a, b) for b in range(16)) for a in range(16))
INVERSE = {a: next(b for b in GL2 if MULTIPLICATION[a][b] == IDENTITY) for a in GL2}
scalar_product = lru_cache(None)(polynomial_multiply)


def coefficient_rows(first, second, e, f, length):
    """Return both complete coefficient cuts of f^T product(K_i) e.

    K_i(t;x,y,u,v)=Steiner(t;x,y,u,v) F_x G_y G_u^-1 F_v^-1.
    Matrices use row-major four-bit encoding; vectors use two-bit encoding.
    """
    if length < 1 or len(first) != 3 or len(second) != 3:
        raise ValueError('Use positive length and two three-letter gauge tables.')
    if any(a not in GL2 for a in (*first, *second)) or (e, f) not in BOUNDARIES:
        raise ValueError('Use invertible binary matrices and a boundary with f^T e=1.')
    scalar = curve(table(3, 'affine'))
    local = []
    for x, y in product(range(3), repeat=2):
        row = []
        for u, v in product(range(3), repeat=2):
            polynomial = scalar[3*x+y][3*u+v]
            if not polynomial:
                continue
            auxiliary = IDENTITY
            for a in (first[x], second[y], INVERSE[second[u]], INVERSE[first[v]]):
                auxiliary = MULTIPLICATION[auxiliary][a]
            row.append((3*u+v, polynomial, auxiliary))
        local.append(row)
    rows = [[(0, 1, IDENTITY)]]
    for _ in range(length):
        rows = [[(9*c+d, scalar_product(p, s), MULTIPLICATION[a][b])
                 for (c, p, a), (d, s, b) in product(left, right)]
                for left, right in product(rows, local)]
    size = len(rows)
    first_cut = [[0]*size for _ in range(3*length+1)]
    second_cut = [[0]*size for _ in first_cut]
    digits = [[(a//9**(length-1-i)) % 9 for i in range(length)] for a in range(size)]
    contraction = []
    for a in range(16):
        image = sum((sum(((a >> (2*i+j)) & 1)*((e >> j) & 1) for j in range(2)) % 2) << i
                    for i in range(2))
        contraction.append((image & f).bit_count() % 2)
    endpoint_entries = 0
    for r, row in enumerate(rows):
        for c, polynomial, auxiliary in row:
            if not contraction[auxiliary]:
                continue
            new_row = new_column = swap = 0
            for a, b in zip(digits[r], digits[c]):
                x, y = divmod(a, 3)
                u, v = divmod(b, 3)
                new_row = 9*new_row+3*x+u
                new_column = 9*new_column+3*y+v
                swap = 9*swap+3*y+x
            assert polynomial.bit_count() % 2 == int(c == swap)
            endpoint_entries += polynomial.bit_count() % 2
            while polynomial:
                term = polynomial & -polynomial
                degree = term.bit_length()-1
                first_cut[degree][r] ^= 1 << c
                second_cut[degree][new_row] ^= 1 << new_column
                polynomial ^= term
    assert endpoint_entries == size
    return first_cut, second_cut
