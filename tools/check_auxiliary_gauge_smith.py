#!/usr/bin/env python3
"""Exact Laurent-gauge and coefficient-rank checks for ordered crossings.

The finite checks support the rational-gauge theorem. They do not certify
an improving crossing family. Python's standard library is sufficient.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random

import check_smith_profile_transfer as basic


def exponents(polynomial):
    return {i for i in range(polynomial.bit_length()) if polynomial >> i & 1}


def product_series(left, right):
    result = set()
    for a in left:
        for b in right:
            degree = a+b
            if degree in result:
                result.remove(degree)
            else:
                result.add(degree)
    return result


def product_matrix(left, right):
    dimension = len(left)
    result = [[set() for _ in range(dimension)] for _ in range(dimension)]
    for a, b, c in itertools.product(range(dimension), repeat=3):
        result[a][c] ^= product_series(left[a][b], right[b][c])
    return result


def identity_series(dimension):
    return [[{0} if a == b else set() for b in range(dimension)]
            for a in range(dimension)]


def gauge(local, heights, infinity=False):
    direction = -1 if infinity else 1
    return {labels: [[{direction*k+heights[a]-heights[b]
                       for k in exponents(value[a][b])}
                      for b in range(len(heights))]
                     for a in range(len(heights))]
            for labels, value in local.items()}


def flatten(local, cut, side):
    dimension = len(next(iter(local.values())))
    matrix = [[set() for _ in range(4)] for _ in range(4*dimension*dimension)]
    for (x, y, u, v), value in local.items():
        left, right = (2*x+y, 2*u+v) if cut == 0 else (2*x+u, 2*y+v)
        for a, b in itertools.product(range(dimension), repeat=2):
            row, column = (left, right) if side == 0 else (right, left)
            matrix[(row*dimension+a)*dimension+b][column] = value[a][b]
    return matrix


def restored_orders(matrix):
    support = [degree for row in matrix for entry in row for degree in entry]
    assert support
    minimum, maximum = min(support), max(support)
    normalized = [[sum(1 << (degree-minimum) for degree in entry)
                   for entry in row] for row in matrix]
    degree = maximum-minimum
    orders = basic.smith_orders(normalized, degree)
    # These four scales independently check the normalization and the
    # complete precision bound; the base checker tests all small scales.
    precisions = sorted({1, 2, degree+1, len(matrix[0])*degree+1})
    for precision in precisions:
        expected = sum(max(0, precision-order) for order in orders)
        assert basic.toeplitz_rank(normalized, precision) == expected
    return [order+minimum for order in orders], len(precisions)


def verify_telescoping(local, heights, length, infinity):
    dimension = len(heights)
    transformed = gauge(local, heights, infinity)
    direction = -1 if infinity else 1
    all_labels = list(local)
    if length <= 2:
        words = itertools.product(all_labels, repeat=length)
    else:
        generator = random.Random(12937+17*heights[-1]+int(infinity))
        words = [tuple(generator.choice(all_labels) for _ in range(length))
                 for _ in range(128)]
    entries = altered_boundaries = 0
    for word in words:
        ordinary = [[int(a == b) for b in range(dimension)] for a in range(dimension)]
        changed = identity_series(dimension)
        for labels in word:
            ordinary = basic.matrix_product(ordinary, local[labels])
            changed = product_matrix(changed, transformed[labels])
        for a, b in itertools.product(range(dimension), repeat=2):
            expected = {direction*k+heights[a]-heights[b]
                        for k in exponents(ordinary[a][b])}
            assert changed[a][b] == expected
            unadjusted = {direction*k for k in exponents(ordinary[a][b])}
            altered_boundaries += expected != unadjusted
            entries += 1
        # The chosen scalar boundary is the first coordinate on both
        # sides. Its product is unchanged, even for negative heights.
        assert changed[0][0] == {direction*k for k in exponents(ordinary[0][0])}
    return entries, altered_boundaries


def check_case(name, local, degree):
    profiles = {}
    toeplitz_checks = telescoping_entries = altered_boundaries = 0
    for height in (-2, 0, 2):
        heights = (0, height)
        for infinity, cut, side in itertools.product((False, True), range(2), range(2)):
            orders, checks = restored_orders(flatten(gauge(local, heights, infinity), cut, side))
            profile = [-order for order in orders] if infinity else orders
            profiles[infinity, cut, side, height] = profile
            toeplitz_checks += checks
        for infinity, length in itertools.product((False, True), (1, 2, 3)):
            entries, changed = verify_telescoping(local, heights, length, infinity)
            telescoping_entries += entries
            altered_boundaries += changed

    rank_checks = count_checks = endpoint_checks = 0
    best = None
    for length in (1, 2, 3):
        ranks, endpoints = basic.full_coefficients(local, degree, length)
        endpoint_checks += endpoints
        for first_cut in range(2):
            first_options = [(side, height) for side in range(2) for height in (-2, 0, 2)]
            for (first_side, first_height), (second_side, second_height) in itertools.product(first_options, repeat=2):
                degrees = profiles[True, first_cut, first_side, first_height]
                orders = profiles[False, 1-first_cut, second_side, second_height]
                high = basic.word_counts(degrees, length)
                low = basic.word_counts(orders, length)
                for s in (Fraction(1), Fraction(3, 2), Fraction(3)):
                    U = sum(s**value for value in degrees)
                    V = sum(s**(-value) for value in orders)
                    B = U*V
                    if best is None or B < best[0]:
                        best = (B, first_cut, first_side, first_height,
                                second_side, second_height, s, degrees, orders)
                    for k in range(degree*length+1):
                        accepted_high = sum(count for weight, count in high.items() if weight >= k)
                        accepted_low = sum(count for weight, count in low.items() if weight <= k)
                        high_splits = max(0, max(degrees)*length-k+1)
                        low_splits = max(0, k-min(orders)*length+1)
                        assert ranks[first_cut][k] <= high_splits*accepted_high
                        assert ranks[1-first_cut][k] <= low_splits*accepted_low
                        assert accepted_high <= s**(-k)*U**length
                        assert accepted_low <= s**k*V**length
                        assert ranks[0][k]*ranks[1][k] <= high_splits*low_splits*B**length
                        rank_checks += 2
                        count_checks += 3
    assert best is not None
    B, cut, first_side, first_height, second_side, second_height, s, degrees, orders = best
    return {
        'name': name, 'degree': degree,
        'toeplitz_checks': toeplitz_checks,
        'telescoping_matrix_entries': telescoping_entries,
        'entries_requiring_boundary_shift': altered_boundaries,
        'original_crossing_endpoint_entries': endpoint_checks,
        'coefficient_rank_checks': rank_checks,
        'rational_count_and_product_checks': count_checks,
        'best_sampled_bound': {'B': str(B), 's': str(s), 'infinity_cut': cut,
                               'infinity_side': first_side, 'infinity_height': first_height,
                               'zero_side': second_side, 'zero_height': second_height,
                               'degrees': degrees, 'orders': orders},
    }


def negative_controls():
    orders, _ = restored_orders([[{-3}, set()], [set(), {-1}]])
    assert orders == [-3, -1]
    assert [-order for order in orders] == [3, 1]
    # Arbitrary transition exponents need not telescope. The two-cycle
    # has nonzero total exponent and cannot be a difference of heights.
    left = [[set(), {0}], [set(), set()]]
    right = [[set(), set()], [{0}, set()]]
    changed_right = [[set(), set()], [{1}, set()]]
    assert product_matrix(left, right)[0][0] == {0}
    assert product_matrix(left, changed_right)[0][0] == {1}
    return {'restored_negative_orders': orders,
            'unrestored_orders_would_be_wrong': [0, 2],
            'nonpotential_cycle_changes_boundary_entry': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [('scalar pencil with auxiliary identity', *basic.elementary(2)),
             ('ordered pairing', *basic.pairing((1, 6, 9), 0)),
             ('noncommuting nilpotent loop', *basic.loop(13, 18))]
    records = [check_case(name, local, degree) for name, local, degree in cases]
    assert sum(record['entries_requiring_boundary_shift'] for record in records) > 0
    report = {
        'status': 'passed', 'field': 2, 'cases': records,
        'negative_controls': negative_controls(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'base_smith_checker_sha256': hashlib.sha256(Path(basic.__file__).read_bytes()).hexdigest(),
        'scope': 'Exact finite gauge, Smith normalization, and upper-bound checks; no improving crossing certificate.',
    }
    print(json.dumps(report, indent=2))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
