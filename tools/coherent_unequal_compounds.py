#!/usr/bin/env python3
"""Complete binary compound profiles with different paired physical alphabets."""
from itertools import product
import time

from coherent_hot_compounds import compound_block, one_particle, support, validate_maps


def profiles(embedding_x, retraction_x, embedding_y, retraction_y, modes, keep_blocks=False):
    px, lx, hx, cx, checks_x = validate_maps(embedding_x, retraction_x, modes)
    py, ly, hy, cy, checks_y = validate_maps(embedding_y, retraction_y, modes)
    px, lx, py, ly = (tuple(map(support, values)) for values in (px, lx, py, ly))
    matrix, started, cuts = one_particle(modes), time.monotonic(), []
    for cut in ('zero', 'infinity'):
        orders, blocks = [], []
        for charge in sorted({a+b for a in cx for b in cy}):
            pairs = [(x, y) for x in range(len(px)) for y in range(len(py)) if cx[x]+cy[y] == charge]
            rows = [tuple(a | (b << modes) for a, b in product(lx[x], ly[y] if cut == 'zero' else py[y]))
                    for x, y in pairs]
            columns = [tuple(a | (b << modes) for a, b in product(py[y] if cut == 'zero' else ly[y], px[x]))
                       for x, y in pairs]
            block = compound_block(matrix, rows, columns, coefficients=keep_blocks)
            block['pairs'] = pairs
            blocks.append(block)
            orders.extend(block['orders_t'])
        assert len(orders) == len(px)*len(py)
        cuts.append(dict(cut=cut, orders_t=sorted(orders), blocks=blocks))
    return dict(modes=modes, q_x=len(px), q_y=len(py), heights_x=list(hx), heights_y=list(hy),
                charges_x=list(cx), charges_y=list(cy), deletion_checks=checks_x+checks_y,
                zero_orders_t=cuts[0]['orders_t'], infinity_orders_t=cuts[1]['orders_t'],
                cuts=cuts, seconds=time.monotonic()-started)
