#!/usr/bin/env python3
"""Polynomial crossings from an asymmetric three-letter involutive core.

Direct outputs have an exact crossing after common physical word-type
restriction. Ramified repair and temporal outputs have a strict local
crossing endpoint. No numerical pressure bound is asserted by this tool.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from check_schur_downset_crossing import sparse_kernel
from weighted_smith_oracle import polynomial_exact_divide, polynomial_multiply


# Extreme states: vacuum, low, top, both. Physical legs: x,y,u,v.
CORE_MAPS = (((0,), (1,), (3,)), ((0,), (2,), (3,)),
             ((0,), (2,), (1, 3)), ((0,), (1,), (2, 3)))


def polynomial_power(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_multiply(result, value)
        exponent >>= 1
        if exponent:
            value = polynomial_multiply(value, value)
    return result


def minor_function(matrix):
    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        row = bit.bit_length() - 1
        remaining = columns
        result = 0
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            column = cb.bit_length() - 1
            if matrix[row][column]:
                result ^= polynomial_multiply(matrix[row][column], minor(rows ^ bit, columns ^ cb))
        return result
    return minor


def choices(label, leg, n):
    bulk = 1 << (n - 2)
    body = (label % bulk) << 1
    core = label // bulk
    return [(body | (state & 1) | ((state >> 1) << (n - 1)),
             int(leg == 2 and core == 2 and state == 1),
             int(leg == 3 and core == 2 and state == 3))
            for state in CORE_MAPS[leg][core]]


def logical_heights(indices):
    n, top = len(indices), indices[-1]
    bulk = 1 << (n - 2)
    return [[sum(indices[i + 1] for i in range(n - 2) if label >> i & 1)
             + ((0, 0, top) if leg in (0, 3) else (0, top, top))[label // bulk]
             for label in range(3 * bulk)] for leg in range(4)]


def coded_layer(indices, ambient, input_power=0, output_power=0, quartic_power=0,
                ramified=False):
    indices = tuple(indices)
    assert len(indices) >= 2 and indices == tuple(sorted(set(indices))) and indices[0] == 0
    assert ambient > indices[-1] and ambient & (ambient - 1) == 0
    n, top = len(indices), indices[-1]
    q = 3 << (n - 2)
    if ramified:
        assert input_power == output_power == quartic_power == 0
        t = 3  # t=1+z, so the final polynomial endpoint is z=1.
        w = 1 ^ polynomial_power(t, top)
        matrix = [[0] * (2 * n) for _ in range(2 * n)]
        for i, a in enumerate(indices):
            for j, b in enumerate(indices):
                if a & ~b == 0:
                    matrix[i][j] = polynomial_power(t, top * (2 * a + 1) + a - b)
                if b & ~a == 0:
                    matrix[i][n + j] = polynomial_multiply(
                        polynomial_power(w, 2 * ambient), polynomial_power(t, a - b))
                matrix[n + i][j] = int(i == j)
                matrix[n + i][n + j] = polynomial_multiply(
                    polynomial_power(t, top + a - b),
                    polynomial_power(w, 2 * (ambient - 1 - (a | b))))
        clearing = 0
    else:
        kernel = sparse_kernel(indices, ambient)
        matrix = [row[:2 * n] for row in kernel[2 * n:]]
        clearing = -min(0, input_power, output_power, input_power + output_power,
                        quartic_power, quartic_power + input_power,
                        quartic_power + output_power, quartic_power + input_power + output_power)
    minor = minor_function(matrix)
    values = [0] * q ** 4
    for x, y, u, v in product(range(q), repeat=4):
        numerator = 0
        for xc, yc, uc, vc in product(choices(x, 0, n), choices(y, 1, n),
                                       choices(u, 2, n), choices(v, 3, n)):
            xx, yy, uu, vv = xc[0], yc[0], uc[0], vc[0]
            rows, columns = vv | (yy << n), uu | (xx << n)
            factor = 1 if ramified or not uc[1] else 3
            shift = 0 if ramified else input_power * uc[1] + output_power * vc[2] + clearing
            gaussian = polynomial_multiply(factor, minor(rows, columns))
            numerator ^= (gaussian if ramified else polynomial_multiply(5, gaussian)) << shift
            if xx & 1 and uu & 1 and yy >> (n - 1) & 1 and vv >> (n - 1) & 1:
                quartic = polynomial_multiply(
                    factor, minor(rows ^ (1 << (n - 1)) ^ (1 << (2 * n - 1)),
                                  columns ^ 1 ^ (1 << n)))
                numerator ^= quartic << (shift + (0 if ramified else quartic_power))
        index = x + q * y + q * q * u + q ** 3 * v
        values[index] = numerator if ramified else polynomial_exact_divide(numerator, 5)
    return values, clearing


def repair_outputs(values, q, bulk, power=0, ramified=False):
    clearing = max(0, -power)
    result = [value << clearing for value in values]
    factor = 1 if ramified else 3
    for x, u, by, bv in product(range(q), range(q), range(bulk), range(bulk)):
        # Explicit pairs are (target_y,target_v,source_y,source_v).
        for ty, tv, sy, sv in ((0, 1, 1, 0), (2, 1, 1, 2)):
            target = x + q * (ty * bulk + by) + q * q * u + q ** 3 * (tv * bulk + bv)
            source = x + q * (sy * bulk + by) + q * q * u + q ** 3 * (sv * bulk + bv)
            result[target] ^= polynomial_multiply(factor, values[source]) << (clearing + power)
    return result, clearing


def compose_core(values, bulk):
    q = 3 * bulk
    grouped_bulk = bulk * bulk
    grouped_q = 3 * grouped_bulk
    result = [0] * grouped_q ** 4

    def index(x, y, u, v):
        return x + q * y + q * q * u + q ** 3 * v

    for x, y, u, v in product(range(grouped_q), repeat=4):
        cores = [label // grouped_bulk for label in (x, y, u, v)]
        first_bulk = [label % bulk for label in (x, y, u, v)]
        second_bulk = [(label // bulk) % bulk for label in (x, y, u, v)]
        value = 0
        for i, j in product(range(3), repeat=2):
            first = values[index(cores[0] * bulk + first_bulk[0], i * bulk + first_bulk[1],
                                 cores[2] * bulk + first_bulk[2], j * bulk + first_bulk[3])]
            second = values[index(j * bulk + second_bulk[0], cores[1] * bulk + second_bulk[1],
                                  i * bulk + second_bulk[2], cores[3] * bulk + second_bulk[3])]
            if first and second:
                value ^= polynomial_multiply(first, second)
        result[x + grouped_q * y + grouped_q ** 2 * u + grouped_q ** 3 * v] = value
    return result, grouped_q


def check_endpoint(values, q, heights=None):
    positive = 0
    for x, y, u, v in product(range(q), repeat=4):
        actual = values[x + q * y + q * q * u + q ** 3 * v].bit_count() & 1
        expected = int(x == v and y == u)
        if heights is None:
            assert actual == expected, (x, y, u, v)
        else:
            delta = heights[1][y] + heights[3][v] - heights[0][x] - heights[2][u]
            if actual:
                assert delta >= 0, (x, y, u, v, delta)
            if delta == 0:
                assert actual == expected, (x, y, u, v)
            positive += actual and delta > 0
    return positive


def build(indices=(0, 1, 2), ambient=4, input_power=0, output_power=0, quartic_power=0,
          repair_power=None, mode='direct'):
    indices = tuple(indices)
    bulk = 1 << (len(indices) - 2)
    q = 3 * bulk
    ramified = mode != 'direct'
    values, clearing = coded_layer(indices, ambient, input_power, output_power, quartic_power, ramified)
    if mode == 'temporal':
        values, q = compose_core(values, bulk)
    elif mode == 'ramified-repair':
        values, _ = repair_outputs(values, q, bulk, ramified=True)
    elif repair_power is not None:
        values, extra = repair_outputs(values, q, bulk, repair_power)
        clearing += extra
    heights = None if ramified else logical_heights(indices)
    positive = check_endpoint(values, q, heights)
    return {'q': q, 'field': 2, 'indices': indices, 'ambient': ambient,
            'mode': mode, 'values': values, 'degree': max(p.bit_length() - 1 for p in values),
            'index_order': 'x+q*y+q^2*u+q^3*v', 'endpoint': 1, 'strict_endpoint': ramified,
            'input_power': input_power, 'output_power': output_power, 'quartic_power': quartic_power,
            'repair_power': repair_power, 'common_w_clearing_power': clearing,
            'logical_heights_by_leg': heights, 'positive_endpoint_entries': positive,
            'endpoint_entries_checked': q ** 4}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--indices', default='0,1,2')
    parser.add_argument('--ambient', type=int, default=4)
    parser.add_argument('--input-power', type=int, default=0)
    parser.add_argument('--output-power', type=int, default=0)
    parser.add_argument('--quartic-power', type=int, default=0)
    parser.add_argument('--repair-power', type=int)
    parser.add_argument('--mode', choices=('direct', 'ramified-repair', 'temporal'), default='direct')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = build(tuple(map(int, args.indices.split(','))), args.ambient, args.input_power,
                 args.output_power, args.quartic_power, args.repair_power, args.mode)
    args.output.write_text(json.dumps(data, separators=(',', ':')) + '\n')
    print(json.dumps({key: value for key, value in data.items() if key != 'values'}, indent=2))


if __name__ == '__main__':
    main()
