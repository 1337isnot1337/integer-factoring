#!/usr/bin/env python3
"""Exact p-group swap endpoints and binary first-jet exclusions."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_steiner_cubic_crossings import coefficient_rows, toeplitz_rank


def cyclic(order):
    return order, lambda a, b: (a+b) % order


def dihedral(half_order):
    def multiply(a, b):
        i, j = divmod(a, 2)
        k, l = divmod(b, 2)
        return 2*((i+(-1)**j*k) % half_order)+(j ^ l)
    return 2*half_order, multiply


def quaternion():
    units = [[(0, 0), (1, 0), (2, 0), (3, 0)],
             [(1, 0), (0, 1), (3, 0), (2, 1)],
             [(2, 0), (3, 1), (0, 1), (1, 0)],
             [(3, 0), (2, 0), (1, 1), (0, 1)]]
    def multiply(a, b):
        unit, sign = units[a % 4][b % 4]
        return unit+4*(sign ^ (a//4) ^ (b//4))
    return 8, multiply


def curve(q, multiply):
    inverse = [next(b for b in range(q) if multiply(a, b) == multiply(b, a) == 0) for a in range(q)]
    pairs, n = list(product(range(q), repeat=2)), q*q
    first = [q*multiply(x, y)+y for x, y in pairs]
    second = [q*x+multiply(y, inverse[x]) for x, y in pairs]
    third = [q*multiply(y, x)+y for x, y in pairs]
    inversion = [q*x+inverse[y] for x, y in pairs]
    swap = [q*y+x for x, y in pairs]
    assert [inversion[third[second[first[j]]]] for j in range(n)] == swap
    for permutation in (first, second, third):
        powered = list(range(n))
        for _ in range(q):
            powered = [permutation[i] for i in powered]
        assert powered == list(range(n))
    matrix = [[0]*n for _ in range(n)]
    for j in range(n):
        column = {j: 1}
        for permutation in (first, second, third):
            updated = column.copy()
            for k, value in column.items():
                updated[k] = updated.get(k, 0) ^ (value << 1)
                target = permutation[k]
                updated[target] = updated.get(target, 0) ^ (value << 1)
            column = {k: value for k, value in updated.items() if value}
        for k, value in column.items():
            matrix[inversion[k]][j] = value
    assert all(matrix[i][j].bit_count() % 2 == int(i == swap[j]) for i, j in product(range(n), repeat=2))
    return matrix


def check():
    records = []
    for name, (q, multiply) in (('C8', cyclic(8)), ('D8', dihedral(4)), ('Q8', quaternion()),
                               ('C16', cyclic(16)), ('D16', dihedral(8)), ('D32', dihedral(16))):
        assert all(multiply(multiply(a, b), c) == multiply(a, multiply(b, c))
                   for a, b, c in product(range(q), repeat=3))
        matrix = curve(q, multiply)
        n = q*q
        column = [[matrix[q*x+y][q*u+v] for y, v in product(range(q), repeat=2)]
                  for x, u in product(range(q), repeat=2)]
        leading = [sum(((value >> 3) & 1) << j for j, value in enumerate(row)) for row in matrix]
        r = binary_rank(leading)
        assert r == binary_rank(leading, low=True)
        coefficients = coefficient_rows(column)
        ranks = [toeplitz_rank(coefficients, k) for k in (1, 2, 3)]
        zero, one, two = ranks[0], ranks[1]-2*ranks[0], ranks[2]-2*ranks[1]+ranks[0]
        assert zero == 1 and one >= 0 and two >= 0
        floor = r*two
        assert floor > q**3
        records.append(dict(group=name, alphabet=q, row_leading_rank=r,
                            column_toeplitz_ranks=ranks,
                            column_zero_order_counts={'0': zero, '1': one, '2': two},
                            exact_all_scale_bare_floor=floor,
                            strict_bare_rho_lower_bound='3/4'))
    return dict(records=records, scope='Displayed six groups, characteristic two and its extensions; bare cubic criterion only.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
