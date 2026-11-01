#!/usr/bin/env python3
"""Independent affine-block counts for Gray-code F4 scalar normalization.

Counts use binary linear systems, not traversal of the individual points.
The rank-product calculations remain the native scanner's separate work.
"""
from __future__ import annotations


def add_equation(pivots, row, rhs):
    while row:
        pivot = row.bit_length()-1
        if pivot not in pivots:
            pivots[pivot] = row, rhs
            return True
        old, value = pivots[pivot]
        row, rhs = row ^ old, rhs ^ value
    return rhs == 0


def affine_zero_counts(start, stop, directions, gray_mask):
    length = stop-start
    if length <= 0 or length & (length-1) or start % length:
        raise ValueError('Use one aligned binary interval.')
    n, free = len(directions), length.bit_length()-1
    if start < 0 or stop > 1 << (2*n) or gray_mask < 0 or gray_mask >> (2*n):
        raise ValueError('The interval and mask must fit the declared coordinate space.')
    local_mask = (1 << free)-1
    width = max((v.bit_length() for v in directions), default=0)
    forms = []
    for coordinate in range(width):
        pair = []
        for plane in range(2):
            full, constant = 0, 0
            for j, direction in enumerate(directions):
                if not direction >> coordinate & 1:
                    continue
                bit = j+plane*n
                full ^= 1 << bit
                if bit+1 < 2*n:
                    full ^= 1 << (bit+1)
                constant ^= (gray_mask >> bit) & 1
            pair.append((full & local_mask, constant ^ ((full & start).bit_count() % 2)))
        forms.append(pair)
    # Higher original coordinates vanish; the leading coordinate equals
    # the F4 digit one exactly when its low/high bits are 1,0.
    pivots, normalized, consistent = {}, 0, True
    for pair in reversed(forms):
        trial = dict(pivots)
        if (add_equation(trial, pair[0][0], pair[0][1] ^ 1)
                and add_equation(trial, pair[1][0], pair[1][1])):
            normalized += 1 << (free-len(trial))
        if (not add_equation(pivots, pair[0][0], pair[0][1])
                or not add_equation(pivots, pair[1][0], pair[1][1])):
            consistent = False
            break
    zeros = (1 << (free-len(pivots))) if consistent else 0
    return {'visited': length, 'normalized_threshold_tests': normalized,
            'zeros': zeros, 'scalar_orbit_skips': length-normalized-zeros}
