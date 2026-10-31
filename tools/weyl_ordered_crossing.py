#!/usr/bin/env python3
"""Qutrit Weyl word crossings over F4, with packed exact coefficient cuts."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_controlled_permutation_smith import binary_rank
from check_steiner_cubic_crossings import curve, table
from weighted_smith_oracle import polynomial_multiply

MULTIPLICATION = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
ROOT_POWERS = (1, 2, 3)
FREE_CHARGES = (1, 2, 3, 4, 6, 8)
scalar_product = lru_cache(None)(polynomial_multiply)


def canonical_weights(code):
    if not 0 <= code < 4096:
        raise ValueError('Use a twelve-bit code for the six free F4 weights.')
    weights = [0]*9
    weights[0], weights[5], weights[7] = 1, 3, 3
    for i, charge in enumerate(FREE_CHARGES):
        weights[charge] = (code >> (2*i)) & 3
    return weights


def coefficient_cuts(length, weights):
    if length < 1 or len(weights) != 9 or weights[0] != 1 or any(v not in range(4) for v in weights):
        raise ValueError('Use positive length and nine F4 weights with w(0,0)=1.')
    scalar = curve(table(3, 'affine'))
    local = []
    for x, y in product(range(3), repeat=2):
        row = []
        for u, v in product(range(3), repeat=2):
            polynomial = scalar[3*x+y][3*u+v]
            if polynomial:
                row.append((3*u+v, polynomial, ((x-v) % 3, (y-u) % 3, -v*(y-u) % 3)))
        local.append(row)
    rows = [[(0, 1, (0, 0, 0))]]
    for _ in range(length):
        rows = [[(9*c+d, scalar_product(p, s),
                  ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3, (a[2]+b[2]+a[1]*b[0]) % 3))
                 for (c, p, a), (d, s, b) in product(left, right)]
                for left, right in product(rows, local)]
    size = len(rows)
    cuts = [[[[0]*size for _ in range(2)] for _ in range(3*length+1)] for _ in range(2)]
    digits = [[(r//9**(length-1-i)) % 9 for i in range(length)] for r in range(size)]
    endpoint_entries = 0
    for r, row in enumerate(rows):
        for c, polynomial, (a, b, phase) in row:
            value = MULTIPLICATION[weights[3*a+b]][ROOT_POWERS[phase]]
            if not value:
                continue
            new_row = new_column = swap = 0
            for first, second in zip(digits[r], digits[c]):
                x, y = divmod(first, 3)
                u, v = divmod(second, 3)
                new_row = 9*new_row+3*x+u
                new_column = 9*new_column+3*y+v
                swap = 9*swap+3*y+x
            endpoint = polynomial.bit_count() % 2
            assert (value if endpoint else 0) == int(c == swap)
            endpoint_entries += endpoint
            while polynomial:
                term = polynomial & -polynomial
                degree = term.bit_length()-1
                for plane in range(2):
                    if value >> plane & 1:
                        cuts[0][degree][plane][r] ^= 1 << c
                        cuts[1][degree][plane][new_row] ^= 1 << new_column
                polynomial ^= term
    assert endpoint_entries == size
    return cuts


def scale(matrix, value):
    if value == 1:
        return matrix
    if value == 2:
        return [matrix[1], [a ^ b for a, b in zip(*matrix)]]
    if value == 3:
        return [[a ^ b for a, b in zip(*matrix)], matrix[0]]
    return [[0]*len(matrix[0]) for _ in range(2)]


def add(first, second):
    return [[a ^ b for a, b in zip(left, right)] for left, right in zip(first, second)]


def boundary_jets():
    # Build affine directions by subtracting the code-zero tensor. The
    # endpoint remains valid in every constructor call.
    baseline = coefficient_cuts(2, canonical_weights(0))
    pick = lambda cuts: [cuts[0][6]]+cuts[1][:4]
    current = pick(baseline)
    changes = []
    for bit in range(12):
        moved = pick(coefficient_cuts(2, canonical_weights(1 << bit)))
        difference = [add(a, b) for a, b in zip(moved, current)]
        assert not any(difference[1][0]) and not any(difference[1][1])
        changes.append(difference)
    previous = 0
    for index in range(4096):
        code = index ^ (index >> 1)
        if index:
            changed = (code ^ previous).bit_length()-1
            current = [add(a, b) for a, b in zip(current, changes[changed])]
        previous = code
        yield code, current


@lru_cache(None)
def spread_masks(width):
    size = 1 << (width-1).bit_length()
    result = []
    step = size//2
    while step:
        mask = sum(((1 << step)-1) << (2*step*i) for i in range(size//step))
        result.append((step, mask))
        step //= 2
    return tuple(result)


def scalar_restriction(matrix):
    masks = spread_masks(len(matrix[0]))
    def spread(value):
        for step, mask in masks:
            value = (value | (value << step)) & mask
        return value
    rows = []
    for low, high in zip(*matrix):
        a, b = spread(low), spread(high)
        rows.extend((a | (b << 1), b | ((a ^ b) << 1)))
    return rows


def field_rank(matrix, low=False):
    result = binary_rank(scalar_restriction(matrix), low=low)
    assert result % 2 == 0
    return result//2


def toeplitz_rank(coefficients, precision, low=False):
    size = len(coefficients[0])
    rows = []
    for i in range(precision):
        for r in range(size):
            rows.append(sum(coefficients[i-j][r] << (size*j)
                            for j in range(max(0, i-len(coefficients)+1), i+1)))
    result = binary_rank(rows, low=low)
    assert result % 2 == 0
    return result//2
