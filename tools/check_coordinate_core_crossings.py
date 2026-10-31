#!/usr/bin/env python3
"""Exact rational coordinate-core constructors and typed polynomial endpoints.

This certifies an admissible family, not an improved factoring exponent.
Polynomial coefficients are bit masks over F2.  Physical entries have index
x + q*y + q**2*u + q**3*v.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path

from check_schur_gaussian_reduction import polynomial_multiply as mul


ROOT = Path(__file__).resolve().parent.parent


def power(base, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, base)
        base = mul(base, base)
        exponent >>= 1
    return result


def exact_divide(numerator, denominator):
    quotient = 0
    while numerator.bit_length() >= denominator.bit_length():
        shift = numerator.bit_length() - denominator.bit_length()
        numerator ^= denominator << shift
        quotient ^= 1 << shift
    assert not numerator
    return quotient


def at_one(polynomial):
    return polynomial.bit_count() & 1


def c_order(polynomial):
    if not polynomial:
        return None
    order = 0
    while not at_one(polynomial):
        polynomial = exact_divide(polynomial, 3)
        order += 1
    return order


def one_particle(indices, ambient):
    """Direct formula, rows (v,y) and columns (u,x)."""
    indices = tuple(indices)
    n = len(indices)
    assert n >= 2 and indices == tuple(sorted(set(indices))) and indices[0] == 0
    assert ambient > indices[-1] and ambient & (ambient - 1) == 0
    matrix = [[0] * (2 * n) for _ in range(2 * n)]
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            matrix[i][j] = mul(3, power(5, a)) if a & b == a else 0
            matrix[i][n + j] = (1 << (2 * ambient)) if a & b == b else 0
            matrix[n + i][j] = int(i == j)
            matrix[n + i][n + j] = 3 << (2 * (ambient - 1 - (a | b)))
    return matrix


def core_maps(n):
    bulk_size = 1 << (n - 2)
    top = 1 << (n - 1)
    # x/v use (low,vacuum,both); y/u use (low,vacuum,top).
    first = [((label % bulk_size) << 1) | (1, 0, top | 1)[label // bulk_size]
             for label in range(3 * bulk_size)]
    second = [((label % bulk_size) << 1) | (1, 0, top)[label // bulk_size]
              for label in range(3 * bulk_size)]
    return first, second, second, first


def build(indices=(0, 1, 2), ambient=4):
    indices = tuple(indices)
    n = len(indices)
    matrix = one_particle(indices, ambient)

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        value = 0
        remaining = columns
        while remaining:
            column = remaining & -remaining
            remaining ^= column
            j = column.bit_length() - 1
            if matrix[i][j]:
                value ^= mul(matrix[i][j], minor(rows ^ bit, columns ^ column))
        return value

    maps = core_maps(n)
    q = len(maps[0])
    values = [0] * q**4
    corrections = [0] * q**4
    top = 1 << (n - 1)
    for x, y, u, v in product(range(q), repeat=4):
        xx, yy, uu, vv = (maps[leg][label]
                          for leg, label in enumerate((x, y, u, v)))
        base = minor(vv | (yy << n), uu | (xx << n))
        correction = 0
        if xx & 1 and uu & 1 and yy & top and vv & top:
            correction = minor((vv ^ top) | ((yy ^ top) << n),
                               (uu ^ 1) | ((xx ^ 1) << n))
        position = x + q*y + q*q*u + q*q*q*v
        values[position] = mul(5, base) ^ correction
        corrections[position] = correction
    height = [sum(indices[j + 1] for j in range(n - 2) if label >> j & 1)
              + indices[-1] * int(label // (1 << (n - 2)) == 2)
              for label in range(q)]
    return {'q': q, 'indices': list(indices), 'ambient': ambient,
            'degree': max(p.bit_length() - 1 for p in values),
            'denominator_degree': 2, 'strict_endpoint': False,
            'typed_endpoint': 'one common full physical word type on all four legs',
            'core_size': 3, 'bulk_size': 1 << (n - 2),
            'core_maps': maps, 'logical_heights': height, 'values': values}, corrections


def independent_exterior_build(indices, ambient):
    """Multiply literal exterior factors, independently of all-minors recursion."""
    n = len(indices)
    matrix = one_particle(indices, ambient)
    support = {0: 1}
    for i, row in enumerate(matrix):
        output = (3*n + i) if i < n else i
        for j, coefficient in enumerate(row):
            if not coefficient:
                continue
            input_index = 2*n + j if j < n else j - n
            monomial = (1 << output) | (1 << input_index)
            for mask, value in list(support.items()):
                if mask & monomial:
                    continue
                new = mask | monomial
                support[new] = support.get(new, 0) ^ mul(coefficient, value)
                if not support[new]:
                    del support[new]
    quartic = 1 | (1 << (2*n - 1)) | (1 << (2*n)) | (1 << (4*n - 1))
    numerator = {mask: mul(5, value) for mask, value in support.items()}
    for mask, value in support.items():
        if not mask & quartic:
            target = mask | quartic
            numerator[target] = numerator.get(target, 0) ^ value
    maps = core_maps(n)
    q = len(maps[0])
    values = [0] * q**4
    for x, y, u, v in product(range(q), repeat=4):
        mask = sum(maps[leg][label] << (n*leg)
                   for leg, label in enumerate((x, y, u, v)))
        values[x + q*y + q*q*u + q*q*q*v] = numerator.get(mask, 0)
    return values


def residual_core():
    # Literal constant exterior factors; variable order x0,xT,y0,yT,u0,uT,v0,vT.
    factors = [(1 << 0) | (1 << 6), (1 << 1) | (1 << 7),
               (1 << 4) | (1 << 2), (1 << 5) | (1 << 3),
               (1 << 5) | (1 << 6), (1 << 1) | (1 << 2),
               (1 << 0) | (1 << 3) | (1 << 4) | (1 << 7)]
    support = {0}
    for factor in factors:
        support ^= {mask | factor for mask in support if not mask & factor}
    first, second = (1, 0, 3), (1, 0, 2)
    rows = [0] * 9
    for x, y, u, v in product(range(3), repeat=4):
        mask = first[x] | (second[y] << 2) | (second[u] << 4) | (first[v] << 6)
        if mask in support:
            rows[v + 3*y] |= 1 << (x + 3*u)
    expected = [1 << i for i in range(9)]
    expected[0] ^= 1 << 5
    expected[3] ^= 1 << 7
    assert rows == expected
    off = [rows[i] ^ (1 << i) for i in range(9)]
    for row in off:
        image = 0
        for j in range(9):
            if row >> j & 1:
                image ^= off[j]
        assert not image
    return rows


def substitute_ramified(polynomial, top):
    """w=1+t**top, without evaluating at any numerical point."""
    result = 0
    for exponent in range(polynomial.bit_length()):
        if polynomial >> exponent & 1:
            sub = exponent
            while True:
                result ^= 1 << (top * sub)
                if not sub:
                    break
                sub = (sub - 1) & exponent
    return result


def check_ramified_endpoint(data, residual):
    q, top = data['q'], data['indices'][-1]
    bulk = data['bulk_size']
    heights = data['logical_heights']
    checks = 0
    for x, y, u, v in product(range(q), repeat=4):
        position = x + q*y + q*q*u + q*q*q*v
        value = substitute_ramified(data['values'][position], top)
        shift = heights[y] + heights[v] - heights[x] - heights[u] - 2*top
        if shift < 0:
            assert not value & ((1 << (-shift)) - 1)
            value >>= -shift
        else:
            value <<= shift
        core_row = (v // bulk) + 3*(y // bulk)
        core_column = (x // bulk) + 3*(u // bulk)
        expected = ((residual[core_row] >> core_column) & 1)
        expected *= int(x % bulk == v % bulk and y % bulk == u % bulk)
        assert value & 1 == expected
        checks += 1
    return checks


def check_two_copy_types(data):
    q = data['q']
    values = data['values']
    checks = 0
    denominator = power(3, 4)
    for pair in combinations_with_replacement(range(q), 2):
        words = sorted(set(permutations(pair)))
        for x, y, u, v in product(words, repeat=4):
            numerator = 1
            for site in range(2):
                position = x[site] + q*y[site] + q*q*u[site] + q*q*q*v[site]
                numerator = mul(numerator, values[position])
            polynomial = exact_divide(numerator, denominator)
            assert at_one(polynomial) == int(x == v and y == u)
            checks += 1
    return checks


def digest(values):
    return sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--seed-directory', type=Path)
    parser.add_argument('--compare-seed', action='append', type=Path, default=[])
    args = parser.parse_args()
    residual = residual_core()
    cases = [((0, 1), 2), ((0, 1, 2), 4), ((0, 1, 3), 4),
             ((0, 2, 3), 4), ((0, 1, 2, 3), 4), ((0, 2, 5), 8)]
    records = []
    cached = {}
    for indices, ambient in cases:
        data, corrections = build(indices, ambient)
        assert data['values'] == independent_exterior_build(indices, ambient)
        q = data['q']
        orders = [c_order(p) for p in data['values'] if p]
        assert min(orders) >= 1
        if len(indices) == 2:
            assert not any(corrections)
            assert min(orders) >= 2
        for x, y, u, v in product(range(q), repeat=4):
            position = x + q*y + q*q*u + q*q*q*v
            correction = corrections[position]
            if not correction:
                continue
            assert c_order(correction) >= 1
            if x // data['bulk_size'] == 2:
                assert c_order(correction) >= 2
        ramified = check_ramified_endpoint(data, residual)
        typed = check_two_copy_types(data)
        pole_entries = sum(order == 1 for order in orders)
        if len(indices) > 2:
            assert pole_entries > 0  # A genuine negative control for local polynomiality.
        records.append({'indices': indices, 'ambient': ambient, 'q': q,
                        'degree': data['degree'], 'denominator_degree': 2,
                        'nonzero_entries': sum(bool(p) for p in data['values']),
                        'minimum_common_c_order': min(orders),
                        'strictly_polar_entries': pole_entries,
                        'numerator_c_order_histogram': dict(sorted(Counter(orders).items())),
                        'ramified_endpoint_entries': ramified,
                        'two_copy_common_type_entries': typed,
                        'independent_exterior_entries': q**4,
                        'coefficient_sha256': digest(data['values'])})
        cached[(indices, ambient)] = data
        if args.seed_directory:
            args.seed_directory.mkdir(parents=True, exist_ok=True)
            tag = '-'.join(map(str, indices))
            (args.seed_directory / f'coordinate-core-{tag}-M{ambient}.json').write_text(
                json.dumps(data, separators=(',', ':')) + '\n')
    comparisons = []
    for path in args.compare_seed:
        supplied = json.loads(path.read_text())
        key = tuple(supplied['indices']), supplied['ambient']
        data = cached[key]
        assert supplied['values'] == data['values']
        assert supplied['denominator_degree'] == 2
        comparisons.append({'name': path.name, 'entries': len(data['values']),
                            'coefficient_sha256': digest(data['values'])})
    result = {'status': 'exact algebra and endpoint checks passed; no exponent improvement claimed',
              'residual_rows_pair_index_v_plus_3y': residual,
              'cases': records, 'independent_supplied_seed_comparisons': comparisons}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
