#!/usr/bin/env python3
"""Exact physical-code restrictions of the uniform implicit reflection tensor.

Only requested Pfaffian coefficients are constructed. The ambient tensor
need not have a feasible complete physical alphabet. This module makes
no rank or pressure claim for the resulting restricted crossing.
"""
from __future__ import annotations

from functools import lru_cache

from check_schur_gaussian_reduction import polynomial_multiply
from reflection_dual_rail_gaussian import gauged_kernel, logical_heights, validate


def restricted_tensor(interval, xletters, yletters=None):
    validate(interval)
    modes = interval // 2
    xletters = tuple(xletters)
    yletters = xletters if yletters is None else tuple(yletters)
    for letters in (xletters, yletters):
        if (not letters or len(set(letters)) != len(letters)
                or any(not isinstance(value, int) or not 0 <= value < 1 << modes
                       for value in letters)):
            raise ValueError('Use nonempty sets of distinct ambient physical labels.')
    kernel = gauged_kernel(interval, 'zero')
    heights = logical_heights(interval, 'zero')
    scalar = interval * interval // 2

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        row = first.bit_length() - 1
        rest = mask ^ first
        remaining, value = rest, 0
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            column = chosen.bit_length() - 1
            if kernel[row][column]:
                value ^= polynomial_multiply(kernel[row][column], pfaffian(rest ^ chosen))
        return value

    physical = (xletters, yletters, yletters, xletters)
    shifts = [[sum(heights[leg * modes + bit] for bit in range(modes)
                   if symbol >> bit & 1) for symbol in letters]
              for leg, letters in enumerate(physical)]
    values = []
    for vi, v in enumerate(xletters):
        for ui, u in enumerate(yletters):
            for yi, y in enumerate(yletters):
                for xi, x in enumerate(xletters):
                    mask = x | (y << modes) | (u << (2 * modes)) | (v << (3 * modes))
                    gauged = pfaffian(mask) << scalar
                    shift = shifts[0][xi] + shifts[1][yi] + shifts[2][ui] + shifts[3][vi]
                    if gauged & ((1 << shift) - 1):
                        raise ArithmeticError('The restricted reflection coefficient has a pole.')
                    value = gauged >> shift
                    if value.bit_count() % 2 != int(xi == vi and yi == ui):
                        raise ArithmeticError('The literal restricted crossing endpoint failed.')
                    values.append((value, 0))
    qx, qy = len(xletters), len(yletters)
    info = pfaffian.cache_info()
    return {'values': values, 'dimensions': [qx, qy, qy, qx], 'q_x': qx, 'q_u': qy,
            'denominator_degree': 0, 'strict_endpoint': True, 'interval': interval,
            'ambient_physical_modes': modes, 'ambient_x_labels': xletters,
            'ambient_y_labels': yletters, 'scalar_w_power': scalar,
            'degree': max(value.bit_length() - 1 for value, _ in values),
            'construction_counts': {'requested_entries': len(values),
                                    'cached_pfaffians': info.currsize,
                                    'pfaffian_cache_hits': info.hits}}
