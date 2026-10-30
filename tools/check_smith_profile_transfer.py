#!/usr/bin/env python3
"""Independent finite checks of the local Smith crossing upper bound."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
from pathlib import Path
import random


def multiply(a, b):
    out = 0
    while b:
        bit = b & -b
        out ^= a << (bit.bit_length()-1)
        b ^= bit
    return out


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            k = row.bit_length()-1
            if k not in pivots:
                pivots[k] = row
                break
            row ^= pivots[k]
    return len(pivots)


def valuation(a):
    assert a
    return (a & -a).bit_length()-1


def smith_orders(matrix, degree, precision=None):
    a = [row[:] for row in matrix]
    m, n = len(a), len(a[0])
    precision = n*degree+1 if precision is None else precision
    mask = (1 << precision)-1
    orders = []
    for k in range(n):
        choices = [(valuation(a[i][j]), i, j)
                   for i in range(k, m) for j in range(k, n) if a[i][j]]
        assert choices, 'generic rank deficiency, not an infinite Smith order certificate'
        order, pi, pj = min(choices)
        orders.append(order)
        a[k], a[pi] = a[pi], a[k]
        for row in a:
            row[k], row[pj] = row[pj], row[k]
        unit = a[k][k] >> order
        inverse = 1
        for j in range(1, precision):
            if multiply(unit, inverse) >> j & 1:
                inverse ^= 1 << j
        assert multiply(unit, inverse) & mask == 1
        for j in range(k, n):
            a[k][j] = multiply(a[k][j], inverse) & mask
        assert a[k][k] == 1 << order
        for i in range(k+1, m):
            factor = a[i][k] >> order
            for j in range(k+1, n):
                a[i][j] ^= multiply(factor, a[k][j]) & mask
            a[i][k] = 0
    return orders


def toeplitz_rank(matrix, precision):
    m, n = len(matrix), len(matrix[0])
    rows = []
    for output_degree in range(precision):
        for i in range(m):
            row = 0
            for input_degree in range(output_degree+1):
                coefficient = output_degree-input_degree
                for j in range(n):
                    if matrix[i][j] >> coefficient & 1:
                        row |= 1 << (input_degree*n+j)
            rows.append(row)
    return rank(rows)


def elementary(auxiliary=1):
    local = {}
    for x, y, u, v in itertools.product(range(2), repeat=4):
        identity = int(x == u and y == v)
        crossing = int(x == v and y == u)
        value = identity | ((identity ^ crossing) << 1)
        local[x, y, u, v] = [[value if a == b else 0
                              for b in range(auxiliary)] for a in range(auxiliary)]
    return local, 1


def loop(first, second):
    assert first and second and not (first & second).bit_count() % 2
    local, _ = elementary(2)
    for x, y, u, v in local:
        for a, b in itertools.product(range(2), repeat=2):
            row, column, swapped = 2*(2*x+y)+a, 2*(2*u+v)+b, 2*(2*v+u)+b
            if first >> row & 1:
                if second >> column & 1:
                    local[x, y, u, v][a][b] ^= 6
                if ((second >> column) ^ (second >> swapped)) & 1:
                    local[x, y, u, v][a][b] ^= 12
    return local, 3


def pairing(matrices, grading):
    powers = (0, 0, 1) if grading == 0 else (0, 1, 2)
    local = {}
    for x, y, u, v in itertools.product(range(2), repeat=4):
        R = int(x == y and u == v)
        C = int(x == u and y == v)
        E = int(x == v and y == u)
        H = R ^ C ^ E
        value = [[0]*2 for _ in range(2)]
        for selected, power, matrix in zip((R, H, C), powers, matrices):
            if not selected:
                continue
            for a, b in itertools.product(range(2), repeat=2):
                polynomial = int(a == b) ^ (3 if matrix >> (2*a+b) & 1 else 0)
                value[a][b] ^= polynomial << power
        local[x, y, u, v] = value
    return local, max(powers)+1


def flatten(local, cut, side):
    auxiliary = len(next(iter(local.values())))
    matrix = [[0]*4 for _ in range(4*auxiliary*auxiliary)]
    for (x, y, u, v), value in local.items():
        p, q = (2*x+y, 2*u+v) if cut == 0 else (2*x+u, 2*y+v)
        for a, b in itertools.product(range(auxiliary), repeat=2):
            row, column = (p, q) if side == 0 else (q, p)
            matrix[(row*auxiliary+a)*auxiliary+b][column] = value[a][b]
    return matrix


def spectra(local, degree):
    result = []
    toeplitz_checks = 0
    for cut, side in itertools.product(range(2), repeat=2):
        matrix = flatten(local, cut, side)
        orders = smith_orders(matrix, degree)
        for h in range(1, 4*degree+2):
            assert toeplitz_rank(matrix, h) == sum(max(0, h-v) for v in orders)
            toeplitz_checks += 1
        result.append(orders)
        reverse = [[sum(((entry >> j) & 1) << (degree-j) for j in range(degree+1))
                    for entry in row] for row in matrix]
        infinity_orders = smith_orders(reverse, degree)
        for h in range(1, 4*degree+2):
            assert toeplitz_rank(reverse, h) == sum(max(0, h-v) for v in infinity_orders)
            toeplitz_checks += 1
        result.append([degree-v for v in infinity_orders])
    return result, toeplitz_checks


def matrix_product(a, b):
    n = len(a)
    out = [[0]*n for _ in range(n)]
    for i, k, j in itertools.product(range(n), repeat=3):
        out[i][j] ^= multiply(a[i][k], b[k][j])
    return out


def full_coefficients(local, degree, length):
    d = 1 << length
    auxiliary = len(next(iter(local.values())))
    row = [[0]*(d*d) for _ in range(degree*length+1)]
    col = [[0]*(d*d) for _ in range(degree*length+1)]
    endpoint_checks = 0
    for x, y, u, v in itertools.product(range(d), repeat=4):
        matrix = [[int(i == j) for j in range(auxiliary)] for i in range(auxiliary)]
        for position in range(length):
            labels = tuple(w >> position & 1 for w in (x, y, u, v))
            matrix = matrix_product(matrix, local[labels])
        polynomial = matrix[0][0]
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        endpoint_checks += 1
        while polynomial:
            bit = polynomial & -polynomial
            k = bit.bit_length()-1
            polynomial ^= bit
            row[k][x*d+y] |= 1 << (u*d+v)
            col[k][x*d+u] |= 1 << (y*d+v)
    return [list(map(rank, row)), list(map(rank, col))], endpoint_checks


def word_counts(values, length):
    counts = {0: 1}
    for _ in range(length):
        following = Counter()
        for total, multiplicity in counts.items():
            for value in values:
                following[total+value] += multiplicity
        counts = following
    return counts


def check_family(local, degree, lengths):
    spectrum, toeplitz_checks = spectra(local, degree)
    rank_checks = endpoint_checks = chernoff_checks = 0
    for length in lengths:
        ranks, endpoints = full_coefficients(local, degree, length)
        endpoint_checks += endpoints
        for first_cut in range(2):
            for first_side, second_side in itertools.product(range(2), repeat=2):
                degrees = spectrum[4*first_cut+2*first_side+1]
                orders = spectrum[4*(1-first_cut)+2*second_side]
                high = word_counts(degrees, length)
                low = word_counts(orders, length)
                for k in range(degree*length+1):
                    count_high = sum(v for s, v in high.items() if s >= k)
                    count_low = sum(v for s, v in low.items() if s <= k)
                    assert ranks[first_cut][k] <= (degree*length-k+1)*count_high
                    assert ranks[1-first_cut][k] <= (k+1)*count_low
                    rank_checks += 2
                    for s in (Fraction(1), Fraction(3, 2), Fraction(3)):
                        U = sum(s**v for v in degrees)
                        V = sum(s**(-v) for v in orders)
                        assert count_high <= s**(-k)*U**length
                        assert count_low <= s**k*V**length
                        assert ranks[0][k]*ranks[1][k] <= (degree*length+1)**2*(U*V)**length
                        chernoff_checks += 3
    return spectrum, toeplitz_checks, rank_checks, endpoint_checks, chernoff_checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    args = parser.parse_args()
    cases = [('scalar elementary', *elementary()),
             ('auxiliary elementary', *elementary(2)),
             ('ordered pairing', *pairing((1, 6, 9), 0)),
             ('higher degree pairing', *pairing((5, 10, 3), 1)),
             ('nilpotent loop', *loop(1, 64)),
             ('mixed nilpotent loop', *loop(13, 18))]
    totals = [0, 0, 0, 0]
    records = []
    for name, local, degree in cases:
        spectrum, *counts = check_family(local, degree, (1, 2, 3))
        totals = [a+b for a, b in zip(totals, counts)]
        records.append({'name': name, 'degree': degree, 'spectra': spectrum})
    baseline = records[0]['spectra']
    assert baseline[1] == [1, 0, 0, -1]
    assert baseline[4] == [0, 1, 1, 2]
    value = sum(Fraction(3)**d for d in baseline[1])*sum(Fraction(3)**(-v) for v in baseline[4])
    assert value == Fraction(256, 27)
    for numerator in range(1, 100):
        s = Fraction(numerator, 7)
        assert (s+1)**4/s**3-Fraction(256, 27) == (s-3)**2*(27*s*s+14*s+3)/(27*s**3)
    # Same two leading coefficients do not determine the full Smith spectrum.
    first = [[1, 0], [0, 1 << 2]]
    second = [[1, 0], [0, 1 << 3]]
    assert toeplitz_rank(first, 2) == toeplitz_rank(second, 2)
    assert smith_orders(first, 3) == [0, 2]
    assert smith_orders(second, 3) == [0, 3]
    # Exact local-alphabet cutoff for the basic Smith envelope.
    numerator, denominator = 577, 816
    assert 2*numerator*numerator > denominator*denominator
    for q in (2, 3, 4):
        n = q*q
        assert (4*n*n)**(n*denominator) > (n+1)**((n+1)*denominator)*q**(4*numerator*n)
    assert 17**13872 < 2**56704
    # Exercise the arithmetic degree constraints with exact rational bounds.
    rng = random.Random(667)
    degree_budget_checks = 0
    for _ in range(240):
        n, degree = 4, rng.randrange(1, 6)
        degrees = [degree]+[rng.randrange(-degree, degree+1) for _ in range(n-1)]
        if sum(degrees) < 0:
            continue
        orders = [0]+[rng.randrange(0, 2*degree+1) for _ in range(n-1)]
        if sum(orders) > n*degree:
            continue
        for s in (Fraction(1), Fraction(4, 3), Fraction(2), Fraction(5)):
            B = sum(s**d for d in degrees)*sum(s**(-v) for v in orders)
            assert B**n*(n+1)**(n+1) >= (4*n*n)**n
            degree_budget_checks += 1
    report = {
        'status': 'passed',
        'block_toeplitz_smith_checks': totals[0],
        'coefficient_rank_upper_bound_checks': totals[1],
        'ordered_mpo_endpoint_entries': totals[2],
        'exact_rational_count_and_product_checks': totals[3],
        'elementary_minimum': {'numerator': 256, 'denominator': 27, 's': 3},
        'degree_budget_rational_checks': degree_budget_checks,
        'exact_alphabet_exclusions': [2, 3, 4],
        'cases': records,
        'scope': 'verifies the implementations and finite instances; the all-length bound is proved in the accompanying note',
    }
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
