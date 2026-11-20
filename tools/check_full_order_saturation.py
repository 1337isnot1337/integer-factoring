#!/usr/bin/env python3
"""Exhaustive control for full order-saturation of cyclic relation kernels.

For a cyclic group C_h and target logs L_i, let K be the relation kernel.
After normalizing gcd(h, L_i)=1, the checker builds diagonal order rows
o_i e_i and, for every prime l dividing h, one unit-aligned pair row per
non-anchor target.  It verifies that every row lies in K and that the
generated lattice equals K, not just its 2-primary part.  This is a finite
control for the proof in research/full-order-saturation.md; no uniform cost
or factoring claim is made.
"""
import argparse
import itertools
import json
from math import gcd

from check_mixed_pair_gap_cycle import row_hnf_index


CASES = ((8, 4), (16, 3), (8, 5), (12, 4),
         (24, 3), (36, 3), (20, 3),
         (9, 3), (15, 3), (25, 3), (27, 3), (45, 3))


def prime_factors(value):
    factors = []
    divisor = 2
    while divisor*divisor <= value:
        if value % divisor == 0:
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


def build_rows(h, logs):
    width = len(logs)
    common = 0
    for value in logs:
        common = gcd(common, gcd(h, value))
    h //= common
    logs = [value//common for value in logs]
    orders = []
    units = []
    for value in logs:
        divisor = gcd(h, value)
        orders.append(h//divisor)
        units.append(value//divisor)
    rows = [[0]*width for _ in range(width)]
    for position, order in enumerate(orders):
        rows[position][position] = order
    for prime in prime_factors(h):
        top = valuation(h, prime)
        anchors = [position for position, value in enumerate(logs)
                   if value % prime]
        assert anchors
        anchor = anchors[0]
        anchor_order = orders[anchor]
        anchor_unit = units[anchor]
        for position, value in enumerate(logs):
            if position == anchor:
                continue
            order = orders[position]
            unit = units[position]
            level = valuation(order, prime)
            assert level <= top
            modulus = prime**level
            row = [0]*width
            row[position] = anchor_unit*(order//modulus)
            row[anchor] -= unit*(anchor_order//modulus)
            rows.append(row)
    return h, logs, rows


def check_case(h, logs):
    reduced_h, reduced_logs, rows = build_rows(h, logs)
    basis, rank, index = row_hnf_index(rows, len(logs))
    assert rank == len(logs)
    assert all(
        sum(value*coefficient for value, coefficient
            in zip(reduced_logs, row)) % reduced_h == 0
        for row in rows)
    assert index % reduced_h == 0
    return index//reduced_h


def check():
    cases = []
    total = 0
    for h, width in CASES:
        count = 0
        for logs in itertools.product(range(h), repeat=width):
            quotient = check_case(h, logs)
            assert quotient == 1
            count += 1
        cases.append({
            'h': h,
            'width': width,
            'cases_checked': count,
            'all_quotients_one': True,
        })
        total += count
    return {
        'cases': cases,
        'cases_checked': total,
        'failures': 0,
        'scope': ('Exhaustive small cyclic-group control for full '
                  'order-saturation from all-prime unit-aligned rows; no '
                  'uniform cost or factoring claim.'),
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
