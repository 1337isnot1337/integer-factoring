#!/usr/bin/env python3
"""Exhaustive finite control for order-derived relation saturation.

For a cyclic group C_h and target logs L_i, let K be the relation kernel
{x : sum L_i x_i = 0 mod h}.  After stripping the common divisor of h and all
L_i, the checker builds diagonal order rows and one pair row per lower
2-adic level plus the top-level pair rows.  It verifies that the generated
lattice is contained in K and that K/L has odd order, so the order-derived
rows 2-saturate K in every tested case.

The full-prime extension also builds the anchor pair rows for every prime
dividing h, using the units modulo the relevant prime powers.  It verifies
that these rows generate all of K, not only its 2-primary part.  These are
finite exhaustive controls for the proofs in
research/two-primary-order-saturation.md.
"""
import argparse
import itertools
import json
import sys
from math import gcd

from check_mixed_pair_gap_cycle import row_hnf_index


CASES = ((8, 4), (16, 3), (8, 5), (12, 4),
         (24, 3), (36, 3), (20, 3))


def two_adic(value):
    return (value & -value).bit_length() - 1 if value else None


def build_rows(h, logs):
    count = len(logs)
    common = 0
    for value in logs:
        common = gcd(common, gcd(h, value))
    h //= common
    logs = [value//common for value in logs]
    orders = []
    units = []
    levels = []
    for value in logs:
        divisor = gcd(h, value)
        order = h//divisor
        orders.append(order)
        units.append(value//divisor)
        levels.append(two_adic(order))
    rows = [[0]*count for _ in range(count)]
    for index, order in enumerate(orders):
        rows[index][index] = order
    if h % 2 == 0:
        top = two_adic(h)
        assert max(levels) == top
        maxima = [index for index, level in enumerate(levels)
                  if level == top]
        anchor = maxima[0]
        for index, level in enumerate(levels):
            if not 1 <= level < top:
                continue
            modulus = 1 << level
            left = units[anchor] % modulus
            right = units[index] % modulus
            row = [0]*count
            row[index] = left*(orders[index]//modulus)
            row[anchor] -= right*(orders[anchor]//modulus)
            rows.append(row)
        for index in maxima[1:]:
            modulus = 1 << top
            left = units[index] % modulus
            right = units[anchor] % modulus
            row = [0]*count
            row[anchor] = left*(orders[anchor]//modulus)
            row[index] -= right*(orders[index]//modulus)
            rows.append(row)
    return h, logs, rows


def prime_factors(value):
    factors = []
    divisor = 2
    while divisor*divisor <= value:
        if value % divisor:
            divisor += 1
            continue
        factors.append(divisor)
        while value % divisor == 0:
            value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    return factors


def valuation(value, prime):
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def build_full_rows(h, logs):
    count = len(logs)
    common = 0
    for value in logs:
        common = gcd(common, gcd(h, value))
    h //= common
    logs = [value//common for value in logs]
    orders = []
    units = []
    gcds = []
    for value in logs:
        divisor = gcd(h, value)
        gcds.append(divisor)
        orders.append(h//divisor)
        units.append(value//divisor)
    rows = [[0]*count for _ in range(count)]
    for index, order in enumerate(orders):
        rows[index][index] = order
    for prime in prime_factors(h):
        anchor = next(index for index, value in enumerate(logs)
                      if value % prime)
        levels = [valuation(order, prime) for order in orders]
        assert levels[anchor] == max(levels)
        for index in range(count):
            if index == anchor:
                continue
            level = levels[index]
            prime_power = prime**level
            coefficient = units[anchor]*(orders[index]//prime_power)
            numerator = coefficient*logs[index]
            assert numerator % logs[anchor] == 0
            other = numerator//logs[anchor]
            row = [0]*count
            row[index] = coefficient
            row[anchor] -= other
            rows.append(row)
    return h, logs, rows


def check_full_case(h, logs):
    reduced_h, reduced_logs, rows = build_full_rows(h, logs)
    assert all(sum(value*coefficient
                   for value, coefficient in zip(reduced_logs, row))
               % reduced_h == 0 for row in rows)
    basis, rank, index = row_hnf_index(rows, len(logs))
    assert rank == len(logs)
    assert index == reduced_h, (h, logs, reduced_h, index)
    return index//reduced_h


def check_case(h, logs):
    reduced_h, reduced_logs, rows = build_rows(h, logs)
    basis, rank, index = row_hnf_index(rows, len(logs))
    assert rank == len(logs)
    assert all(sum(value*coefficient
                   for value, coefficient in zip(reduced_logs, row))
               % reduced_h == 0 for row in rows)
    assert index % reduced_h == 0
    quotient = index//reduced_h
    return quotient


def check():
    cases = []
    total = 0
    for h, width in CASES:
        count = 0
        for logs in itertools.product(range(h), repeat=width):
            quotient = check_case(h, logs)
            assert quotient & 1
            assert check_full_case(h, logs) == 1
            count += 1
        cases.append({
            'h': h,
            'width': width,
            'cases_checked': count,
            'all_2_primary_quotients_odd': True,
            'all_full_lattices_equal_kernel': True,
        })
        total += count
    return {
        'cases': cases,
        'cases_checked': total,
        'full_lattice_cases_checked': total,
        'failures': 0,
        'scope': ('Exhaustive small cyclic-group control for the order-derived '
                  'two-primary and full-prime saturation theorems; no uniform '
                  'cost or factoring claim.'),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    arguments = parser.parse_args()
    receipt = check()
    text = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if arguments.output:
        with open(arguments.output, 'w') as handle:
            handle.write(text)
    print(text)


if __name__ == '__main__':
    main()
