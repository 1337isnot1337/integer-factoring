#!/usr/bin/env python3
"""Independent exact controls and a bounded regular-cyclic-algebra search."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import exp, log
from pathlib import Path

from check_controlled_permutation_smith import binary_rank, toeplitz_rank
from cyclic_algebra_crossing import coefficient_rows, construct, default_cycle


def matrix_multiply(first, second):
    result = []
    for row in first:
        value = 0
        while row:
            bit = row & -row
            value ^= second[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return result


def direct_operator(q, g):
    """Evaluate g(I+U0) by literal binary matrix multiplication."""
    n = q*q
    identity = [1 << i for i in range(n)]
    u0 = [0] * n
    for source, target in enumerate(default_cycle(q)):
        u0[target] |= 1 << source
    s0 = [a ^ b for a, b in zip(identity, u0)]
    power, s = identity, [0] * n
    for degree in range(n):
        if g >> degree & 1:
            s = [a ^ b for a, b in zip(s, power)]
        power = matrix_multiply(power, s0)
    assert not any(power)
    u = [a ^ b for a, b in zip(identity, s)]
    # g(S0) is regular: its last nonzero power has rank one.
    power = identity
    for degree in range(n):
        if degree in {1, n//2, n-1}:
            assert binary_rank(power) == n-degree
        power = matrix_multiply(power, s)
    assert not any(power)
    return u, s


def phi(operator, q):
    return [sum(((operator[q*x+u] >> (q*v+y)) & 1) << (q*u+v)
                for u in range(q) for v in range(q))
            for x in range(q) for y in range(q)]


def jet_orders(matrix, precision):
    coefficients = [coefficient_rows(matrix, j) for j in range(precision)]
    ranks = [0] + [toeplitz_rank(coefficients, k) for k in range(1, precision+1)]
    increments = [b-a for a, b in zip(ranks, ranks[1:])]
    assert increments[-1] == len(matrix), 'increase the jet precision'
    orders, previous = [], 0
    for order, count in enumerate(increments):
        orders.extend([order] * (count-previous))
        previous = count
    return orders


def first_jet_floor(matrix):
    coefficients = [coefficient_rows(matrix, j) for j in (0, 1)]
    r = binary_rank(coefficients[0])
    k = toeplitz_rank(coefficients, 2) - 2*r
    assert k >= 0
    n = len(matrix)
    return n*r + (n-1)*k, r, k


def controls():
    records = []
    for q, masks in ((2, (2, 6, 10, 14)), (4, (2, 6, 14, 130, 65534)),
                     (8, (2, 6, 14, 3314))):
        for g in masks:
            n = q*q
            row, column, infinity = construct(q, g)
            for x in range(q):
                for y in range(q):
                    for u in range(q):
                        for v in range(q):
                            value = row[q*x+y][q*u+v]
                            assert value.bit_count() % 2 == int(x == v and y == u)
                            assert column[q*x+u][q*y+v] == value
            operator, s = direct_operator(q, g)
            power = operator
            # Independently compare every polynomial coefficient.
            for j in range(n):
                assert coefficient_rows(infinity, j) == phi(power, q)
                power = matrix_multiply(power, s)
            assert not any(power)
            orders = jet_orders(infinity, 10)
            column_orders = None
            if q <= 4:
                column_orders = jet_orders(column, n+1)
                assert column_orders == [0] + [n]*(n-1)
            floor, r, k = first_jet_floor(infinity)
            assert (r, k) == (orders.count(0), orders.count(1))
            records.append(dict(q=q, g=g, row_infinity_orders=orders,
                                column_zero_orders=column_orders,
                                first_jet_floor=floor))
    # Invalid transforms must not silently use the regular-nilpotent theorem.
    for g in (0, 1, 3, 4, 5, 16):
        try:
            construct(2, g)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid transform accepted')
    return records


def numerical_minimum(orders, q):
    """Exploratory convex minimization; the exact first-jet floor certifies scope."""
    n = q*q
    degrees = [n-1-v for v in orders]
    def value(t):
        peak = max(degrees)*t
        weights = [exp(d*t-peak) for d in degrees]
        tail = (n-1)*exp(-n*t)
        slope = sum(d*w for d, w in zip(degrees, weights))/sum(weights)-n*tail/(1+tail)
        cost = peak+log(sum(weights))+log(1+tail)
        return cost/(4*log(q)), slope
    low, high = 0., 1.
    while value(high)[1] < 0:
        high *= 2
    for _ in range(70):
        middle = (low+high)/2
        if value(middle)[1] < 0:
            low = middle
        else:
            high = middle
    return value((low+high)/2)[0]


def search(bits):
    from native_matrix_smith_f4 import MatrixOracle
    q, n, hist, minimum = 8, 64, Counter(), None
    best = None
    for mask in range(1 << bits):
        g = 2 + (mask << 2)
        _, _, matrix = construct(q, g)
        orders = MatrixOracle([[(v, 0) for v in row] for row in matrix],
                              n, precision=256).profile()['orders']
        assert all(v == round(v) for v in orders)
        orders = tuple(map(round, orders))
        if orders not in hist:
            assert jet_orders(matrix, max(orders)+1) == list(orders)
        floor, r, k = first_jet_floor(matrix)
        assert (r, k) == (orders.count(0), orders.count(1))
        # Exact certificate rho >= 31/32, well above the half-rate threshold.
        assert floor**8 >= q**31
        if minimum is None or (floor, g) < minimum:
            minimum = (floor, g)
        rho = numerical_minimum(orders, q)
        if best is None or (rho, g) < best:
            best = (rho, g)
        hist[orders] += 1
    return dict(q=q, transforms=1 << bits, masks=f'g=2+4k, 0 <= k < 2^{bits}',
                profiles=[dict(orders=list(p), count=c) for p, c in sorted(hist.items())],
                minimum_exact_first_jet_floor=minimum[0], minimum_floor_mask=minimum[1],
                certified_rho_lower_bound='31/32', numerical_best_rho=best[0],
                numerical_best_mask=best[1], native_precision=256,
                independent_full_profile_representatives=len(hist),
                scope='Only this base cycle, mask schedule, and bare two-place criterion.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search-bits', type=int, choices=range(1, 15))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(controls=controls(), factoring_exponent='1/2; unchanged')
    if args.search_bits is not None:
        result['search'] = search(args.search_bits)
    result['source_sha256'] = {path.name: sha256(path.read_bytes()).hexdigest()
                              for path in (Path(__file__), Path(__file__).with_name('cyclic_algebra_crossing.py'))}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(controls=len(result['controls']), search=result.get('search', {}).get('transforms'),
                          factoring_exponent=result['factoring_exponent'])))


if __name__ == '__main__':
    main()
