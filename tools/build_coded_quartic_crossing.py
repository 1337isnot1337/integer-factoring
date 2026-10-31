#!/usr/bin/env python3
"""Build a typed polynomial crossing from arbitrary modes and a residual code.

Physical polynomial masks use x+q*y+q^2*u+q^3*v indexing.  A coefficient
argument POLYNOMIAL:SHIFT represents w^SHIFT times the binary polynomial
mask POLYNOMIAL.  Every supplied coefficient must evaluate to one at w=1.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_schur_downset_crossing import sparse_kernel
from weighted_smith_oracle import polynomial_exact_divide, polynomial_multiply


def coefficient(text):
    fields = text.split(':')
    if len(fields) not in (1, 2):
        raise argparse.ArgumentTypeError('Use a polynomial mask, optionally followed by :SHIFT.')
    try:
        polynomial = int(fields[0], 0)
        shift = int(fields[1]) if len(fields) == 2 else 0
    except ValueError as error:
        raise argparse.ArgumentTypeError('Polynomial masks and shifts must be integers.') from error
    if polynomial < 1 or polynomial.bit_count() % 2 != 1:
        raise argparse.ArgumentTypeError('The polynomial must evaluate to one at w=1.')
    return polynomial, shift


def build(indices, ambient, coefficients=((1, 0),) * 4, quartic=(1, 0)):
    indices = tuple(indices)
    assert len(indices) >= 2 and indices == tuple(sorted(set(indices)))
    assert indices[0] == 0 and ambient > indices[-1]
    assert ambient & (ambient - 1) == 0
    assert all(p > 0 and p.bit_count() % 2 == 1 for p, _ in (*coefficients, quartic))
    n = len(indices)
    q = 1 << (n - 1)
    bulk = 1 << (n - 2)
    kernel = sparse_kernel(indices, ambient)
    matrix = [row[:2 * n] for row in kernel[2 * n:]]

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        row = bit.bit_length() - 1
        remaining = columns
        value = 0
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            column = cb.bit_length() - 1
            if matrix[row][column]:
                value ^= polynomial_multiply(matrix[row][column], minor(rows ^ bit, columns ^ cb))
        return value

    def input_choices(label):
        body = (label & (bulk - 1)) << 1
        if label < bulk:
            return [(body, 0)]
        return [(body | 1, 1), (body | 1 | (1 << (n - 1)), 0)]

    def output_choices(label):
        body = (label & (bulk - 1)) << 1
        if label < bulk:
            return [(body, 0)]
        return [(body | (1 << (n - 1)), 0), (body | 1 | (1 << (n - 1)), 1)]

    factors = [(1, 0)] * 16
    for mask in range(1, 16):
        bit = mask & -mask
        previous = factors[mask ^ bit]
        current = coefficients[bit.bit_length() - 1]
        factors[mask] = (polynomial_multiply(previous[0], current[0]), previous[1] + current[1])
    clearing = -min(0, min(shift for _, shift in factors),
                    min(shift + quartic[1] for _, shift in factors))
    values = [0] * q ** 4
    heights = [sum(indices[i + 1] for i in range(n - 2) if label >> i & 1)
               + indices[-1] * (label >> (n - 2)) for label in range(q)]
    nonzero_quartic = 0
    diagonal_quartic = []
    positive_endpoint = 0
    for x, y, u, v in product(range(q), repeat=4):
        gaussian = perturbation = 0
        for (xx, bx), (yy, gy), (uu, bu), (vv, gv) in product(
                input_choices(x), output_choices(y), input_choices(u), output_choices(v)):
            mask = bx | (bu << 1) | (gy << 2) | (gv << 3)
            p, shift = factors[mask]
            c_factor = (1, 3, 5)[bx + bu]
            factor = polynomial_multiply(c_factor, p)
            rows, columns = vv | (yy << n), uu | (xx << n)
            gaussian ^= polynomial_multiply(factor, minor(rows, columns)) << (shift + clearing)
            if xx & 1 and uu & 1 and yy >> (n - 1) & 1 and vv >> (n - 1) & 1:
                numerator = polynomial_multiply(
                    factor, minor(rows ^ (1 << (n - 1)) ^ (1 << (2 * n - 1)),
                                  columns ^ 1 ^ (1 << n)))
                quotient = polynomial_exact_divide(numerator, 5)
                perturbation ^= polynomial_multiply(quotient, quartic[0]) << (shift + clearing + quartic[1])
        value = gaussian ^ perturbation
        index = x + q * y + q * q * u + q ** 3 * v
        values[index] = value
        delta = heights[y] + heights[v] - heights[x] - heights[u]
        endpoint = value.bit_count() & 1
        if endpoint:
            assert delta >= 0, (x, y, u, v, delta)
        if delta == 0:
            assert endpoint == int(x == v and y == u), (x, y, u, v)
        positive_endpoint += endpoint and delta > 0
        nonzero_quartic += bool(perturbation)
        if x == y == u == v and perturbation:
            diagonal_quartic.append(x)
    return {'q': q, 'field': 2, 'indices': indices, 'ambient': ambient,
            'degree': max(value.bit_length() - 1 for value in values),
            'values': values, 'index_order': 'x+q*y+q^2*u+q^3*v',
            'endpoint': 1, 'strict_endpoint': False,
            'endpoint_kind': 'nonnegative logical height; equality is crossing',
            'logical_energy': heights,
            'coefficients_hx_hu_gy_gv': coefficients, 'quartic_numerator': quartic,
            'common_w_clearing_power': clearing,
            'local_endpoint_entries_checked': q ** 4,
            'strictly_positive_endpoint_entries': positive_endpoint,
            'nonzero_quartic_entries': nonzero_quartic,
            'nonzero_quartic_diagonal_labels': diagonal_quartic}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--indices', default='0,1,2,3')
    parser.add_argument('--ambient', type=int, default=4)
    for name in ('hx', 'hu', 'gy', 'gv', 'quartic'):
        parser.add_argument('--' + name, type=coefficient, default=(1, 0))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = build(tuple(map(int, args.indices.split(','))), args.ambient,
                 tuple(getattr(args, name) for name in ('hx', 'hu', 'gy', 'gv')), args.quartic)
    rendered = json.dumps(data, separators=(',', ':')) + '\n'
    args.output.write_text(rendered)
    summary = {key: value for key, value in data.items() if key != 'values'}
    summary['output_sha256'] = sha256(rendered.encode()).hexdigest()
    summary['output'] = str(args.output)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
