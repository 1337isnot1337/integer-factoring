#!/usr/bin/env python3
"""Exact four-mode constant-flag controls over F2 and F4.

The finite exclusion keeps the displayed zero gauge and physical weights.
It is not a general obstruction to other flags, gauges, or mode dimensions.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

from weighted_smith_oracle import WeightedSmithOracle as BinaryOracle
from weighted_smith_oracle_f4 import (WeightedSmithOracle, polynomial_exact_divide,
                                     polynomial_multiply, scalar_multiply)


ROOT = Path(__file__).resolve().parent.parent


def matrix_inverse(matrix):
    n = len(matrix)
    work = [row[:] + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if work[i][j])
        work[j], work[pivot] = work[pivot], work[j]
        inverse = (0, 1, 3, 2)[work[j][j]]
        work[j] = [scalar_multiply(value, inverse) for value in work[j]]
        for i in range(n):
            if i != j and work[i][j]:
                coefficient = work[i][j]
                work[i] = [a ^ scalar_multiply(b, coefficient) for a, b in zip(work[i], work[j])]
    return [row[n:] for row in work]


def constant_product(first, second):
    n = len(first)
    result = [[0] * n for _ in range(n)]
    for i, j, k in product(range(n), repeat=3):
        result[i][j] ^= scalar_multiply(first[i][k], second[k][j])
    return result


def c_power(exponent):
    result = 0
    subset = exponent
    while True:
        result ^= 1 << (2 * subset)
        if not subset:
            return result
        subset = (subset - 1) & exponent


def flag_matrix(n, code):
    matrix = [[int(i == j) for j in range(n)] for i in range(n)]
    binary = True
    for i in range(n):
        for j in range(i):
            matrix[i][j] = code & 3
            binary &= matrix[i][j] < 2
            code >>= 2
    return matrix, binary


def crossing_matrix(n, flag):
    pascal = [[int(j & ~i == 0) for j in range(n)] for i in range(n)]
    transpose = list(map(list, zip(*pascal)))
    inverse_v = constant_product(flag, transpose)
    v = matrix_inverse(inverse_v)
    first = [[0] * n for _ in range(n)]
    second = [[0] * n for _ in range(n)]
    for i, j in product(range(n), repeat=2):
        first[i][j] = scalar_multiply(c_power(2 * i + 1), v[i][j])
        for k in range(n):
            if pascal[k][j]:
                second[i][j] ^= scalar_multiply(c_power(2 * n - 2 * k - 1), inverse_v[i][k])
    # Verify AB=(1+w)^(2n) P before any Smith calculation.
    for i, j in product(range(n), repeat=2):
        value = 0
        for k in range(n):
            value ^= polynomial_multiply(first[i][k], second[k][j])
        assert value == (c_power(2 * n) if pascal[i][j] else 0)
    matrix = [first[i] + [(1 << (4 * n)) if pascal[i][j] else 0 for j in range(n)]
              for i in range(n)]
    matrix += [[int(i == j) for j in range(n)] + second[i] for i in range(n)]
    return matrix, second


def independent_product(first, second):
    # Direct coefficient convolution, independent of the packed multiplier.
    multiplication = ((0, 0, 0, 0), (0, 1, 2, 3),
                      (0, 2, 3, 1), (0, 3, 1, 2))
    result = 0
    for i in range((first.bit_length() + 1) // 2):
        for j in range((second.bit_length() + 1) // 2):
            result ^= multiplication[(first >> (2 * i)) & 3][(second >> (2 * j)) & 3] << (2 * (i + j))
    return result


def independent_determinant(matrix):
    result = 0
    for order in permutations(range(len(matrix))):
        term = 1
        for i, j in enumerate(order):
            term = independent_product(term, matrix[i][j])
        result ^= term
    return result


def minor_prefix_controls(matrix, row_weights, column_weights, orders):
    checked = 0
    for size in (1, 2, 3):
        minimum = None
        for rows in combinations(range(len(matrix)), size):
            for columns in combinations(range(len(matrix)), size):
                value = independent_determinant([[matrix[i][j] for j in columns] for i in rows])
                if not value:
                    continue
                order = -((value.bit_length() - 1) // 2)
                order += sum(row_weights[i] for i in rows) + sum(column_weights[j] for j in columns)
                minimum = order if minimum is None else min(minimum, order)
                checked += 1
        assert minimum == sum(orders[:size])
    return checked


def binary_polynomial(value):
    result = 0
    for degree in range((value.bit_length() + 1) // 2):
        coefficient = (value >> (2 * degree)) & 3
        assert coefficient in (0, 1)
        result |= coefficient << degree
    return result


def check():
    arithmetic_checks = 0
    for first, second in product(range(256), repeat=2):
        value = polynomial_multiply(first, second)
        assert value == independent_product(first, second)
        if second:
            assert polynomial_exact_divide(value, second) == first
        arithmetic_checks += 1
    n = 4
    rows = [2 * i - n + 1 for i in range(n)] + [-n] * n
    columns = [2 * i - n + 1 for i in range(n)] + [n] * n
    baseline = tuple(range(-2 * n + 1, 2 * n, 2))
    all_profiles, binary_profiles = Counter(), Counter()
    minor_checks = 0
    for code in range(4 ** 6):
        flag, binary = flag_matrix(n, code)
        matrix, second = crossing_matrix(n, flag)
        for i, j in product(range(n), repeat=2):
            if second[i][j]:
                valuation = ((second[i][j] & -second[i][j]).bit_length() - 1) // 2
                assert valuation >= 2 * n - 2 - 2 * i - 2 * j
        orders, prefixes = WeightedSmithOracle(matrix, infinity=True).profile(rows, columns)
        assert prefixes[-1].polynomial == 1
        assert sum(orders) == sum(baseline) == 0
        assert all(sum(orders[:k]) <= sum(baseline[:k]) for k in range(1, 2 * n))
        all_profiles[orders] += 1
        if binary:
            binary_matrix = [[binary_polynomial(value) for value in row] for row in matrix]
            binary_orders, _ = BinaryOracle(binary_matrix, infinity=True).profile(rows, columns)
            assert orders == binary_orders
            binary_profiles[orders] += 1
        if code in (0, 1, 37, 4095):
            minor_checks += minor_prefix_controls(matrix, rows, columns, orders)
    assert sum(binary_profiles.values()) == 64
    assert len(all_profiles) == len(binary_profiles) == 9
    assert all_profiles[baseline] == 1
    return {
        'status': 'exact fixed-gauge finite flag control; no new factoring exponent',
        'modes': n, 'fields': [2, 4], 'binary_flags': 64, 'F4_flags': 4096,
        'baseline_infinity_orders': baseline,
        'zero_orders_each_direction': list(range(1, 2 * n, 2)),
        'physical_mode_weights': [8 * i - 4 * (n - 1) for i in range(n)],
        'every_profile_majorizes_baseline': True,
        'unique_normalized_flag_attaining_baseline': True,
        'independent_coefficient_convolution_checks': arithmetic_checks,
        'independent_nonzero_minor_checks': minor_checks,
        'independent_binary_oracle_comparisons': 64,
        'profiles': [{'orders': orders, 'F4_count': count,
                      'F2_count': binary_profiles[orders]}
                     for orders, count in sorted(all_profiles.items())],
        'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ('tools/check_constant_flag_crossings.py',
                                       'tools/weighted_smith_oracle_f4.py',
                                       'tools/weighted_smith_oracle.py')},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
