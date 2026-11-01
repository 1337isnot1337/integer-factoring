#!/usr/bin/env python3
"""Strict crossings from the selected compounds of a hot one-particle matrix.

The variable is t=w**2. The ordinary row profile is the scalar-cancelled
natural profile of the corresponding deletion-downset reflection code.
The opposite cut is a separate computation; no pressure is inherited.
"""
from __future__ import annotations

from hot_code_charge_profiles import natural_profiles


def compound_crossing(xletters, yletters=None, *, modes=None):
    data = natural_profiles(xletters, yletters, modes=modes, keep_blocks=True)
    xs, ys = data['x_letters'], data['y_letters']
    qx, qy = len(xs), len(ys)
    values = [(0, 0)] * (qx * qx * qy * qy)
    pairs = [(ix, iy, x.bit_count() + y.bit_count())
             for ix, x in enumerate(xs) for iy, y in enumerate(ys)]
    for block in data['blocks']:
        chosen = [(ix, iy) for ix, iy, charge in pairs if charge == block['charge']]
        if len(chosen) != block['size']:
            raise ArithmeticError('The complete occupation block has the wrong size.')
        for row, (x, y) in enumerate(chosen):
            for col, (v, u) in enumerate(chosen):
                values[x + qx * (y + qy * (u + qy * v))] = (block['matrix'][row][col], 0)
    for index, (a, b) in enumerate(values):
        x = index % qx
        y = index // qx % qy
        u = index // (qx * qy) % qy
        v = index // (qx * qy * qy) % qx
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('The compound construction changed the strict crossing endpoint.')
    vacuum = xs.index(0) + qx * (ys.index(0) + qy * (ys.index(0) + qy * xs.index(0)))
    if values[vacuum] != (1, 0):
        raise ArithmeticError('The vacuum coefficient must be one, proving scalar content one.')
    return {'q_x': qx, 'q_u': qy, 'dimensions': [qx, qy, qy, qx],
            'degree': max(value[0].bit_length() - 1 for value in values),
            'values': values, 'x_letters': xs, 'y_letters': ys,
            'modes': data['modes'], 'variable': 't=w^2',
            'ordinary_row_orders_t': [order // 2 for order in data['orders_w']],
            'scope': 'Strict compound crossing; only its ordinary row profile is inherited.'}
