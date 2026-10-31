#!/usr/bin/env python3
"""Exact finite controls for the grouped normalized-inverse theorem."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path
import random

from boolean_cauchy_operator import (
    boundary_rows, kernel_matrices, monomial_operator, rotate_subset,
    smith_orders, unit_inverse,
)


def groups(length, center):
    coordinate = (abs(center) & -abs(center)).bit_length() - 1
    return ([x for x in range(length) if not x >> coordinate & 1],
            [x for x in range(length) if x >> coordinate & 1])


def flatten(matrix, length, center):
    first, second = groups(length, center)
    count = 1 << length
    first_labels = [sum(((state >> x) & 1) << i for i, x in enumerate(first))
                    for state in range(count)]
    second_labels = [sum(((state >> x) & 1) << i for i, x in enumerate(second))
                     for state in range(count)]
    result = [[0] * count for _ in range(count)]
    for row in range(count):
        for column in range(count):
            left = (first_labels[row] << len(first)) | first_labels[column]
            right = (second_labels[row] << len(second)) | second_labels[column]
            result[left][right] = matrix[row][column]
    return result


def normalized_group_inverse(physical, center, scalar, precision):
    length = len(physical)
    modulus = 1 << precision
    operator = monomial_operator(physical, modulus)
    inverse = unit_inverse(operator, modulus)
    boundary = sum(1 << x for x in boundary_rows(length, center))
    b = abs(center)
    normalized = [[entry * pow(scalar, b - (column & boundary).bit_count(), modulus) % modulus
                   for column, entry in enumerate(row)] for row in inverse]
    return flatten(normalized, length, center)


def expected_orders(length, boundary_count, scalar):
    if not scalar:
        units = 1 << (length - boundary_count)
        return [0] * units + [None] * ((1 << length) - units)
    step = (abs(scalar) & -abs(scalar)).bit_length() - 1
    return [j * step for j in range(boundary_count + 1)
            for _ in range(comb(boundary_count, j) * (1 << (length - boundary_count)))]


def valuation(value, precision):
    return (value & -value).bit_length() - 1 if value else precision


def tensor_product(left, right, modulus):
    return [[left_row[i] * right_row[j] % modulus
             for i in range(len(left)) for j in range(len(right))]
            for left_row in left for right_row in right]


def check():
    rng = random.Random(502713)
    precision = 48
    modulus = 1 << precision
    weighted_entries = 0
    records = []

    # These are general near-identity physical matrices, not Cauchy kernels.
    for length in (2, 4):
        count = 1 << length
        for _ in range(8):
            physical = [[int(x == y) + 2 * rng.randrange(-3, 4)
                         for y in range(length)] for x in range(length)]
            operator = monomial_operator(physical, modulus)
            inverse = unit_inverse(operator, modulus)
            for selected in range(count):
                occupations = [(state & selected).bit_count() for state in range(count)]
                for matrix in (operator, inverse):
                    for row in range(count):
                        for column in range(count):
                            assert valuation(matrix[row][column], precision) >= abs(
                                occupations[row] - occupations[column])
                            weighted_entries += 1

    # Small Cauchy profiles cover all retained centers at unit points.
    cases = [(length, center, point) for length in (2, 4)
             for center in range(-length + 1, length) if center
             for point in (1, 3)]
    # Outer centers additionally cover arbitrarily different valuation scales,
    # and the genuinely rank-deficient zero-point limit.
    cases += [(length, sign * (length - 1), point)
              for length in (2, 4) for sign in (-1, 1) for point in (0, 2, 16)]
    cases += [(8, center, point) for center, point in
              [(-7, 0), (-7, 2), (-3, 1), (1, 3), (5, 1), (7, 1), (7, 16)]]
    for length, center, point in cases:
        _, physical = kernel_matrices(length, center, point, precision)
        grouped = normalized_group_inverse(physical, center, 2 * point, precision)
        observed = smith_orders(grouped, precision)
        expected = expected_orders(length, abs(center), 2 * point)
        assert observed == expected, (length, center, point, Counter(observed), Counter(expected))
        records.append({'kind': 'Cauchy', 'L': length, 't': center, 'zeta': point,
                        'orders': {str(key): value for key, value in Counter(observed).items()}})

    # Test the more general physical-matrix statement, including outer points
    # with high valuation. Boundary indices are still defined by the center.
    for length in (2, 4):
        for center in range(1, length):
            for trial in range(4):
                physical = [[int(y == (x + center) % length) + 2 * rng.randrange(-3, 4)
                             for y in range(length)] for x in range(length)]
                scalars = [2] if center != length - 1 else [0, 2, 8]
                for scalar in scalars:
                    grouped = normalized_group_inverse(physical, center, scalar, precision)
                    observed = smith_orders(grouped, precision)
                    assert observed == expected_orders(length, abs(center), scalar)
                    records.append({'kind': 'general_physical', 'L': length, 't': center,
                                    'trial': trial, 'scalar': scalar,
                                    'orders': {str(key): value for key, value in Counter(observed).items()}})

    # A separate check of the outer proof's exact charge pattern and ordering.
    charge_controls = 0
    for length in (2, 4):
        for center in (1 - length, length - 1):
            _, physical = kernel_matrices(length, center, 2, precision)
            inverse = unit_inverse(monomial_operator(physical, modulus), modulus)
            count = 1 << length
            first, second = groups(length, center)
            boundary = set(boundary_rows(length, center))
            unique = next(x for x in range(length) if x not in boundary)
            if unique not in first:
                first, second = second, first
            first_mask = sum(1 << x for x in first)
            # Define E directly from its original inverse entry, avoiding the
            # packed flattening routine used in the Smith checks.
            local = [state for state in range(count) if not state & ~first_mask]
            labels = [(row, column) for row in local for column in local]
            def data(label):
                row, column = label
                d0 = len(first) - 1 - column.bit_count() + int(bool(column >> unique & 1))
                d1 = len(second) - row.bit_count()
                return d0, d1, row.bit_count() - column.bit_count()
            labels.sort(key=lambda label: data(label)[:2])
            for ri, (jr, ir) in enumerate(labels):
                d0r, d1r, qr = data((jr, ir))
                for si, (js, is_) in enumerate(labels):
                    d0s, d1s, qs = data((js, is_))
                    original_row = jr | rotate_subset(is_, length, center)
                    original_column = ir | rotate_subset(js, length, -center)
                    entry = inverse[original_row][original_column]
                    if qr < qs:
                        assert entry == 0
                    if ri < si and entry:
                        assert d1r <= d1s
                    assert d0r - d1r == len(first) - 1 - len(second) + qr + int(bool(ir >> unique & 1))
                    charge_controls += 1

    products = []
    for left_center, right_center in [(1, 1), (-1, 1), (1, -1)]:
        _, left_physical = kernel_matrices(2, left_center, 1, precision)
        _, right_physical = kernel_matrices(2, right_center, 3, precision)
        left = normalized_group_inverse(left_physical, left_center, 2, precision)
        right = normalized_group_inverse(right_physical, right_center, 6, precision)
        product = tensor_product(left, right, modulus)
        expected = expected_orders(4, 2, 2)
        assert smith_orders(product, precision) == expected
        needed_precision = 3
        reduced = [[value % (1 << needed_precision) for value in row] for row in product]
        at_needed_precision = smith_orders(reduced, needed_precision)
        assert at_needed_precision == expected
        assert all(order is not None for order in at_needed_precision)
        products.append({'centers': [left_center, right_center], 'slots': 4,
                         'normalization_valuation': 2, 'tested_precision': needed_precision,
                         'nonzero_factors': 16, 'orders': dict(Counter(expected))})

    root = Path(__file__).resolve().parent.parent
    sources = ['tools/check_grouped_cauchy_inverse.py', 'tools/boolean_cauchy_operator.py',
               'tools/cauchy_packet_polynomials.py', 'research/grouped-cauchy-inverse-profile.md']
    return {'status': 'pass', 'precision': precision,
            'weighted_boolean_entry_inequalities': weighted_entries,
            'grouped_smith_profiles': len(records),
            'outer_charge_and_ordering_controls': charge_controls,
            'records': records, 'tensor_products': products,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in sources},
            'scope': 'Exact finite controls for all-size grouped inverse Smith theorems. '
                     'The tensor-product obstruction concerns full packet operators and '
                     'direct groupwise sums of products, not arbitrary factoring algorithms.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.write_text(json.dumps({'status': 'fail', 'error': repr(error)}, indent=2) + '\n')
        raise
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
