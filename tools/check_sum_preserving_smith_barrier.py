#!/usr/bin/env python3
"""Finite controls for the sum-preserving conjugation Smith-envelope floor."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import random

import check_smith_profile_transfer as basic


def binary_product(left, right):
    out = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            value ^= right[bit.bit_length()-1]
            row ^= bit
        out.append(value)
    return out


def binary_inverse(matrix):
    q = len(matrix)
    rows = [row | (1 << (q+i)) for i, row in enumerate(matrix)]
    for column in range(q):
        pivot = next(i for i in range(column, q) if rows[i] >> column & 1)
        rows[pivot], rows[column] = rows[column], rows[pivot]
        for i in range(q):
            if i != column and rows[i] >> column & 1:
                rows[i] ^= rows[column]
    assert [row & ((1 << q)-1) for row in rows] == [1 << i for i in range(q)]
    return [row >> q for row in rows]


def centralizer(q, g, generator, permutation):
    identity = [1 << i for i in range(q)]
    exchange = [1 << (i ^ g) for i in range(q)]
    if permutation:
        if g == 0:
            mapping = list(range(q))
            generator.shuffle(mapping)
        else:
            representatives = [i for i in range(q) if i < (i ^ g)]
            images = representatives[:]
            generator.shuffle(images)
            mapping = [0]*q
            for first, second in zip(representatives, images):
                second ^= g if generator.randrange(2) else 0
                mapping[first], mapping[first ^ g] = second, second ^ g
        result = [1 << i for i in mapping]
    else:
        result = identity[:]
        representatives = [i for i in range(q) if i <= (i ^ g)]
        for _ in range(2*q+3):
            change = identity[:]
            if g == 0:
                first, second = generator.sample(range(q), 2)
                change[first] ^= 1 << second
            elif len(representatives) >= 2 and generator.randrange(2):
                first, second = generator.sample(representatives, 2)
                change[first] ^= 1 << second
                change[first ^ g] ^= 1 << (second ^ g)
            else:
                first, second = generator.choices(representatives, k=2)
                change[first] ^= (1 << second) | (1 << (second ^ g))
                change[first ^ g] ^= (1 << second) | (1 << (second ^ g))
            assert binary_product(change, change) == identity
            result = binary_product(change, result)
    inverse = binary_inverse(result)
    assert binary_product(result, inverse) == identity
    assert binary_product(result, exchange) == binary_product(exchange, result)
    return result, inverse


def polynomial_conjugate(matrix, left, right):
    q = len(matrix)
    intermediate = [[0]*q for _ in range(q)]
    result = [[0]*q for _ in range(q)]
    for i in range(q):
        for j in range(q):
            if left[i] >> j & 1:
                for k in range(q):
                    intermediate[i][k] ^= matrix[j][k]
    for i in range(q):
        for j in range(q):
            for k in range(q):
                if right[j] >> k & 1:
                    result[i][k] ^= intermediate[i][j]
    return result


def original_block(q, g):
    weight = g.bit_count()
    powers = [1]
    for _ in range(weight):
        powers.append(basic.multiply(powers[-1], 3))
    return [[0 if (x ^ u) & ~g else
             powers[weight-(x ^ u).bit_count()] << (x ^ u).bit_count()
             for u in range(q)] for x in range(q)]


def reverse(polynomial, degree):
    return sum(1 << (degree-i) for i in range(degree+1) if polynomial >> i & 1)


def expected_row_degrees(width):
    out = Counter({0: 1})
    for _ in range(width):
        new = Counter()
        for degree, count in out.items():
            for shift, multiplicity in ((-1, 1), (0, 2), (1, 1)):
                new[degree+shift] += count*multiplicity
        out = new
    return out


def check_case(width, seed, permutation):
    q, n = 1 << width, 1 << (2*width)
    generator = random.Random(seed)
    blocks = []
    row_degrees = []
    endpoint_entries = changed_complement_entries = 0
    toeplitz_checks = 0
    for g in range(q):
        left, right = centralizer(q, g, generator, permutation)
        block = polynomial_conjugate(original_block(q, g), left, right)
        blocks.append(block)
        degree = g.bit_count()
        assert all(entry.bit_length() <= degree+1 for row in block for entry in row)
        for x in range(q):
            for u in range(q):
                assert block[x][u] & 1 == int(x == u)
                assert block[x][u].bit_count() % 2 == int(x ^ u == g)
                endpoint_entries += 2
            complement = block[x][x ^ g]
            assert complement and basic.valuation(complement) <= degree
            if permutation:
                assert complement == 1 << degree
            changed_complement_entries += complement != 1 << degree
        reversed_block = [[reverse(entry, degree) for entry in row] for row in block]
        orders = basic.smith_orders(reversed_block, degree)
        row_degrees.extend(degree-order for order in orders)
        if g in (0, q-1):
            for precision in sorted({1, degree+1, q*degree+1}):
                assert basic.toeplitz_rank(reversed_block, precision) == sum(
                    max(0, precision-order) for order in orders)
                toeplitz_checks += 1
    assert Counter(row_degrees) == expected_row_degrees(width)

    column_orders = []
    reverse_column_degrees = []
    forced_order_checks = 0
    for w in range(q):
        block = [[blocks[x ^ y][x][x ^ w] for y in range(q)] for x in range(q)]
        orders = basic.smith_orders(block, width)
        column_orders.extend(orders)
        assert min(orders) <= w.bit_count()
        assert all(block[x][x ^ w] and basic.valuation(block[x][x ^ w]) <= w.bit_count()
                   for x in range(q))
        forced_order_checks += q+1
        reversed_block = [[reverse(entry, width) for entry in row] for row in block]
        reverse_orders = basic.smith_orders(reversed_block, width)
        reverse_column_degrees.extend(width-order for order in reverse_orders)
        if w in (0, q-1):
            for precision in (1, width+1, q*width+1):
                assert basic.toeplitz_rank(block, precision) == sum(
                    max(0, precision-order) for order in orders)
                toeplitz_checks += 1
    assert sum(column_orders) <= n*width
    assert sum(reverse_column_degrees) >= 0

    rational_checks = 0
    for s in (Fraction(1), Fraction(4, 3), Fraction(3, 2), Fraction(2),
              Fraction(5, 2), Fraction(3), Fraction(7, 2), Fraction(4), Fraction(8)):
        U = sum(s**value for value in row_degrees)
        V = sum(s**(-value) for value in column_orders)
        assert U == (s+2+1/s)**width
        assert V >= (1+1/s)**width
        assert V >= (4/s)**width
        assert U*V >= Fraction(64, 9)**width
        assert sum(s**value for value in reverse_column_degrees) >= n
        rational_checks += 5
    return {
        'width': width, 'labels': q,
        'mode': 'permutation' if permutation else 'general linear centralizer',
        'endpoint_entry_checks': endpoint_entries,
        'forced_entry_and_minimum_order_checks': forced_order_checks,
        'complement_entries_changed_by_conjugation': changed_complement_entries,
        'independent_toeplitz_checks': toeplitz_checks,
        'exact_rational_inequality_checks': rational_checks,
    }


def spike_control():
    width, q = 3, 8
    blocks = [original_block(q, g) for g in range(q)]
    blocks[-1] = [[sum(1 << (2*k) for k in range(width+1) if entry >> k & 1)
                   for entry in row] for row in blocks[-1]]
    w = q-1
    block = [[blocks[x ^ y][x][x ^ w] for y in range(q)] for x in range(q)]
    orders = basic.smith_orders(block, 2*width)
    assert min(orders) == 2*width > w.bit_count()
    return {'width': width, 'spiked_block': w,
            'orders': orders, 'old_degree_cap_fails': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [check_case(width, 1709+41*seed+width, permutation)
             for width in (1, 2, 3) for seed in range(4)
             for permutation in (False, True)]
    assert sum(case['complement_entries_changed_by_conjugation'] for case in cases) > 0
    assert 2**317 > 3**200
    assert 2*283**2 > 400**2
    # Check both elementary branches of the exact minimization on a
    # rational grid. Their monotonicity is proved in the note.
    for numerator in range(16, 257):
        s = Fraction(numerator, 16)
        base = (s+2+1/s)*max(1+1/s, 4/s)
        assert base >= Fraction(64, 9)
    report = {
        'status': 'passed', 'field': 2, 'cases': cases,
        'strict_threshold_certificate': {'2^317 > 3^200': True,
                                         '2*283^2 > 400^2': True},
        'spike_negative_control': spike_control(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'base_checker_sha256': hashlib.sha256(Path(basic.__file__).read_bytes()).hexdigest(),
        'scope': 'Finite controls for an all-width Smith-envelope obstruction; no actual-rank or factoring lower bound.',
    }
    print(json.dumps(report, indent=2))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
