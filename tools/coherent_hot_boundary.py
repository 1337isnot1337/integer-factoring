#!/usr/bin/env python3
"""Exact paired compression of a final equal-energy physical block.

NumPy object arrays accelerate only exact bit-plane XOR and F4 products.
The first low letters are retained literally. Boundary maps are stored
as left rows and right columns, as in physical_split_crossings.retract.
"""
from __future__ import annotations
import numpy as np
from physical_split_crossings import check_split


def compress(seed, low, left, right, *, swap_y=False):
    ambient = seed['q_x']
    if seed['q_u'] != ambient or not isinstance(low, int) or not 0 <= low < ambient:
        raise ValueError('Use equal ambient alphabets and a nonempty final boundary block.')
    if check_split(left, right) != ambient - low:
        raise ValueError('The paired split must act on the final boundary block.')
    width = len(left)
    dimensions = [ambient] * 4
    planes = [np.asarray([value[plane] for value in seed['values']], dtype=object)
              .reshape(dimensions, order='F') for plane in (0, 1)]
    maps = [left, right if swap_y else left, left if swap_y else right, right]
    for axis, rows in enumerate(maps):
        moved = [np.moveaxis(plane, axis, 0) for plane in planes]
        outputs = [list(plane[:low]) for plane in moved]
        for row in rows:
            a = np.zeros(moved[0].shape[1:], dtype=object)
            b = np.zeros(moved[0].shape[1:], dtype=object)
            for index, scalar in enumerate(row):
                if scalar:
                    x, y = moved[0][low + index], moved[1][low + index]
                    if scalar == 1:
                        a ^= x; b ^= y
                    elif scalar == 2:
                        a ^= y; b ^= x ^ y
                    else:
                        a ^= x ^ y; b ^= x
            outputs[0].append(a); outputs[1].append(b)
        planes = [np.moveaxis(np.stack(output), 0, axis) for output in outputs]
        dimensions[axis] = low + width
    values = list(zip(*(map(int, plane.ravel(order='F')) for plane in planes)))
    q = low + width
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // q ** axis) % q for axis in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('The compressed tensor lost its strict crossing endpoint.')
    return {'q_x': q, 'q_u': q, 'dimensions': dimensions, 'values': values,
            'degree': max(max(a.bit_length(), b.bit_length()) - 1 for a, b in values),
            'strict_endpoint': True}
