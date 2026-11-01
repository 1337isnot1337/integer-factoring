#!/usr/bin/env python3
"""Polynomial codimension-one physical frames with exact endpoint identity."""

from __future__ import annotations

import numpy as np

from moving_boundary_crossings import evaluate_at_one, multiply_arrays


def compress_frame(seed, a, b, a_shifts=None, b_shifts=None):
    """L=[I,a], P=[I+a*b^T;b^T], so LP=I in characteristic two.

    Rank-one updates avoid constructing or applying dense polynomial maps.
    Neither a represented lower coordinate core nor homogeneity is assumed.
    """
    q = len(a)
    if len(b) != q or seed['dimensions'] != [q + 1] * 4:
        raise ValueError('Supply two vectors and an ambient alphabet one dimension larger.')
    a_shifts = [0] * q if a_shifts is None else list(a_shifts)
    b_shifts = [0] * q if b_shifts is None else list(b_shifts)
    if len(a_shifts) != q or len(b_shifts) != q or any(type(v) is not int for v in a_shifts + b_shifts):
        raise ValueError('Supply one integer Laurent shift per vector coordinate.')
    left_shift = max([0] + [-shift for value, shift in zip(a, a_shifts) if value != (0, 0)])
    right_shift = max([0] + [left_shift - shift for value, shift in zip(b, b_shifts) if value != (0, 0)])
    dimensions = [q + 1] * 4
    planes = [np.asarray([value[p] for value in seed['values']], dtype=object)
              .reshape(dimensions, order='F') for p in (0, 1)]
    for axis in range(4):
        moved = [np.moveaxis(plane, axis, 0) for plane in planes]
        common = [plane[-1] << left_shift for plane in moved]
        if axis >= 2:
            for i, scalar in enumerate(a):
                if scalar == (0, 0):
                    continue
                shift = left_shift + a_shifts[i]
                term = multiply_arrays([moved[0][i], moved[1][i]], (scalar[0] << shift, scalar[1] << shift))
                common[0] ^= term[0]
                common[1] ^= term[1]
        outputs = [[], []]
        for i, scalar in enumerate(a if axis < 2 else b):
            base_shift = left_shift if axis < 2 else right_shift
            if scalar == (0, 0):
                outputs[0].append(moved[0][i] << base_shift)
                outputs[1].append(moved[1][i] << base_shift)
                continue
            if axis < 2:
                shift = a_shifts[i]
                # common already carries left_shift; use an unshifted last coordinate if needed.
                source = [plane[-1] for plane in moved]
                shift += left_shift
            else:
                source = common
                shift = right_shift + b_shifts[i] - left_shift
            term = multiply_arrays(source, (scalar[0] << shift, scalar[1] << shift))
            outputs[0].append((moved[0][i] << base_shift) ^ term[0])
            outputs[1].append((moved[1][i] << base_shift) ^ term[1])
        planes = [np.moveaxis(np.stack(output), 0, axis) for output in outputs]
        dimensions[axis] = q
    values = list(zip(*(map(int, plane.ravel(order='F')) for plane in planes)))
    for index, value in enumerate(values):
        x, y, u, v = ((index // q**axis) % q for axis in range(4))
        if evaluate_at_one(value) != int(x == v and y == u):
            raise ArithmeticError('A polynomial frame crossing endpoint failed.')
    return {'values': values, 'dimensions': dimensions, 'q_x': q, 'q_u': q,
            'degree': max(part.bit_length() - 1 for value in values for part in value),
            'strict_endpoint': True, 'laurent_clearing_power': 2 * (left_shift + right_shift)}
