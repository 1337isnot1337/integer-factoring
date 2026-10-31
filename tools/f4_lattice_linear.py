#!/usr/bin/env python3
"""Packed homogeneous linear algebra over F4, with optional F16 endpoint lifts.

An equation is a pair of integer bit planes.  Bit j of the first/second
plane is the coefficient of 1/omega in column j, where omega^2=omega+1.
The same row space is valid after a fixed coefficient-field extension.
"""
from __future__ import annotations

from itertools import product


MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
INV = (0, 1, 3, 2)


def add(first, second):
    return first[0] ^ second[0], first[1] ^ second[1]


def scale(row, coefficient):
    a, b = row
    if coefficient == 0:
        return 0, 0
    if coefficient == 1:
        return row
    if coefficient == 2:
        return b, a ^ b
    if coefficient == 3:
        return a ^ b, a
    raise ValueError('Use a coefficient in F4.')


def coefficient(row, column):
    return ((row[0] >> column) & 1) | (((row[1] >> column) & 1) << 1)


def encode(values):
    return (sum((value & 1) << j for j, value in enumerate(values)),
            sum((value >> 1) << j for j, value in enumerate(values)))


def dot(first, second):
    a, b = first
    c, d = second
    low = ((a & c).bit_count() ^ (b & d).bit_count()) & 1
    high = ((a & d).bit_count() ^ (b & c).bit_count() ^ (b & d).bit_count()) & 1
    return low | (high << 1)


def eliminate(rows):
    pivots = {}
    for row in rows:
        while row[0] | row[1]:
            pivot = (row[0] | row[1]).bit_length() - 1
            value = coefficient(row, pivot)
            if pivot in pivots:
                row = add(row, scale(pivots[pivot], value))
            else:
                pivots[pivot] = scale(row, INV[value])
                break
    return pivots


def reduce_row(row, pivots):
    for pivot in sorted(pivots, reverse=True):
        value = coefficient(row, pivot)
        if value:
            row = add(row, scale(pivots[pivot], value))
    return row


def endpoint_remainders(endpoint_masks, pivots):
    return tuple(reduce_row((mask, 0), pivots) for mask in endpoint_masks)


def image_basis(remainders):
    """Return (endpoint vector, free-variable assignment) pairs over F4."""
    columns = 0
    for low, high in remainders:
        columns |= low | high
    basis = {}
    while columns:
        bit = columns & -columns
        free = bit.bit_length() - 1
        columns ^= bit
        vector = tuple(coefficient(row, free) for row in remainders)
        lift = bit, 0
        while any(vector):
            pivot = max(j for j, value in enumerate(vector) if value)
            value = vector[pivot]
            if pivot in basis:
                old, old_lift = basis[pivot]
                vector = tuple(a ^ MUL[value][b] for a, b in zip(vector, old))
                lift = add(lift, scale(old_lift, value))
            else:
                inverse = INV[value]
                basis[pivot] = (tuple(MUL[inverse][entry] for entry in vector),
                                scale(lift, inverse))
                break
    return tuple(basis.values())


def f16_multiply(first, second):
    """F16=F4[u]/(u^2+u+omega), with two F4 digits."""
    a, b = first & 3, first >> 2
    c, d = second & 3, second >> 2
    low = MUL[a][c] ^ MUL[2][MUL[b][d]]
    high = MUL[a][d] ^ MUL[b][c] ^ MUL[b][d]
    return low | (high << 2)


def f16_scale_f4(value, scalar):
    return MUL[value & 3][scalar] | (MUL[value >> 2][scalar] << 2)


def full_support_lift(remainders, pivots, field_size=16):
    """Construct a nullspace vector whose endpoint coordinates are nonzero.

    For at most field_size coordinates, every coordinate projection being
    nonzero suffices.  The greedy construction below needs strictly fewer
    coordinates; at equality or above it uses a complete finite image
    search.  The present callers use at most eight coordinates.
    """
    if field_size not in (4, 16):
        raise ValueError('The supported coefficient fields are F4 and F16.')
    if any(not (a | b) for a, b in remainders):
        return None
    basis = image_basis(remainders)
    size = len(remainders)
    endpoint = [0] * size
    free = {}
    if size < field_size:
        for vector, lift in basis:
            union = [j for j in range(size) if endpoint[j] or vector[j]]
            if all(endpoint[j] for j in union):
                continue
            scalar = next(value for value in range(field_size)
                          if all(endpoint[j] ^ f16_scale_f4(value, vector[j]) for j in union))
            endpoint = [a ^ f16_scale_f4(scalar, b) for a, b in zip(endpoint, vector)]
            bits = lift[0] | lift[1]
            while bits:
                bit = bits & -bits
                column = bit.bit_length() - 1
                bits ^= bit
                free[column] = free.get(column, 0) ^ f16_scale_f4(scalar, coefficient(lift, column))
            if all(endpoint):
                break
        assert all(endpoint)
    else:
        # Normalize the first nonzero image coefficient to 1.  This removes
        # scalar duplicates without losing any full-support image vector.
        found = None
        for first in range(len(basis)):
            for tail in product(range(field_size), repeat=len(basis) - first - 1):
                scalars = (0,) * first + (1,) + tail
                values = [0] * size
                for scalar, (vector, _) in zip(scalars, basis):
                    values = [a ^ f16_scale_f4(scalar, b) for a, b in zip(values, vector)]
                if all(values):
                    found = scalars, values
                    break
            if found:
                break
        if found is None:
            return None
        scalars, endpoint = found
        for scalar, (_, lift) in zip(scalars, basis):
            bits = lift[0] | lift[1]
            while bits:
                bit = bits & -bits
                column = bit.bit_length() - 1
                bits ^= bit
                free[column] = free.get(column, 0) ^ f16_scale_f4(scalar, coefficient(lift, column))
    assignment = {column: value for column, value in free.items() if value}
    for pivot, equation in sorted(pivots.items()):
        value = 0
        bits = (equation[0] | equation[1]) ^ (1 << pivot)
        while bits:
            bit = bits & -bits
            column = bit.bit_length() - 1
            bits ^= bit
            value ^= f16_scale_f4(assignment.get(column, 0), coefficient(equation, column))
        if value:
            assignment[pivot] = value
    return tuple(endpoint), assignment, len(basis)


def verify_assignment(equations, assignment):
    for equation in equations:
        value = 0
        bits = equation[0] | equation[1]
        while bits:
            bit = bits & -bits
            column = bit.bit_length() - 1
            bits ^= bit
            value ^= f16_scale_f4(assignment.get(column, 0), coefficient(equation, column))
        if value:
            return False
    return True
