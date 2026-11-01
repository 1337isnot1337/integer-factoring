#!/usr/bin/env python3
"""Independent exact controls for homogeneous generic Smith specialization."""
from __future__ import annotations

import argparse
import ctypes
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from cube_generic_smith128 import (FIELD_MODULUS, WORD_MASK, field_product,
                                  generic_smith_profile, native_library)
from cube_homogeneous_gaussian import (ZERO, ONE, axis_parts, identity, kernel,
                                      matrix_product, monomial, multiply, valuation)


def coupled_control_block(axis):
    """One explicit cyclic two-mode kernel, with physical heights 1 and 3."""
    heights = (1, 3)
    coefficients = [[[ZERO] * 2 for _ in range(2)] for _ in range(8)]
    for selected in range(3):
        for i, height in enumerate(heights):
            powers = [0, 0, 0]
            powers[selected] = 6 * height
            coefficients[1 << selected][i][i] = monomial(*powers)
    for first, second in ((0, 1), (1, 2), (2, 0)):
        block = coefficients[(1 << first) | (1 << second)]
        for i, height in enumerate(heights):
            powers = [0, 0, 0]
            powers[first] = powers[second] = 3 * height
            block[i][i] = monomial(*powers)
        balanced, changed = [0, 0, 0], [0, 0, 0]
        balanced[first] = balanced[second] = 6
        changed[first], changed[second] = 5, 7
        block[0][1] = block[1][0] = monomial(*balanced) ^ monomial(*changed)
    for i, height in enumerate(heights):
        coefficients[7][i][i] = monomial(2 * height, 2 * height, 2 * height)
    change = monomial(4, 4, 4)
    triple = (3, 4, 5)
    for shift in range(3):
        change ^= monomial(*(triple[(selected + shift) % 3] for selected in range(3)))
    coefficients[7][0][1] = coefficients[7][1][0] = change
    for block in coefficients[1:]:
        for i in range(2):
            for j in range(2):
                assert len(block[i][j]) % 2 == int(i == j)
                assert all(sum(term) == 3 * (heights[i] + heights[j])
                           for term in block[i][j])
    others = [selected for selected in range(3) if selected != axis]
    blocks = [coefficients[1 << axis],
              coefficients[(1 << axis) | (1 << others[0])],
              coefficients[(1 << axis) | (1 << others[1])], coefficients[7]]
    result = [[ZERO] * 8 for _ in range(8)]
    for row in range(4):
        for column in range(4):
            if column & ~row == 0:
                for i in range(2):
                    for j in range(2):
                        result[2 * row + i][2 * column + j] = blocks[row ^ column][i][j]
    return result


def complete_minor_profile(matrix, axis):
    """Independent polynomial determinants of every nonempty square minor."""
    order = len(matrix)

    @lru_cache(None)
    def determinant(rows, columns):
        if not rows:
            return ONE
        bit = rows & -rows
        row, rest, scan, value = bit.bit_length() - 1, rows ^ bit, columns, ZERO
        while scan:
            selected = scan & -scan
            column = selected.bit_length() - 1
            scan ^= selected
            value ^= multiply(matrix[row][column], determinant(rest, columns ^ selected))
        return value

    previous, orders, count = 0, [], 0
    for size in range(1, order + 1):
        masks = [mask for mask in range(1 << order) if mask.bit_count() == size]
        minimum = None
        for rows in masks:
            for columns in masks:
                value = determinant(rows, columns)
                count += 1
                if value:
                    candidate = valuation(value, axis)
                    minimum = candidate if minimum is None else min(minimum, candidate)
        assert minimum is not None
        orders.append(minimum - previous)
        previous = minimum
    assert orders == sorted(orders)
    return {'orders': orders, 'nonempty_minors_evaluated': count}


def polynomial_gcd(first, second):
    while second:
        remainder = first
        while remainder and remainder.bit_length() >= second.bit_length():
            remainder ^= second << (remainder.bit_length() - second.bit_length())
        first, second = second, remainder
    return first


def field_controls(rng):
    value = 2
    for step in range(128):
        value = field_product(value, value)
        if step == 63:
            assert polynomial_gcd(value ^ 2, FIELD_MODULUS) == 1
    assert value == 2
    library = native_library()
    for _ in range(1000):
        first, second = rng.getrandbits(128), rng.getrandbits(128)
        a = (ctypes.c_uint64 * 2)(first & WORD_MASK, first >> 64)
        b = (ctypes.c_uint64 * 2)(second & WORD_MASK, second >> 64)
        result = (ctypes.c_uint64 * 2)()
        library.cube_gf128_product(a, b, result)
        assert (result[0] | (result[1] << 64)) == field_product(first, second)
    return {'rabin_squarings': 128, 'rabin_prime_divisor_gcds': 1,
            'independent_field_products': 1000}


def known_profile_matrix(order, rng):
    left, right = identity(order), identity(order)
    for i in range(order):
        for j in range(i + 1, order):
            if rng.randrange(2):
                left[i][j] = monomial(rng.randrange(3), 0, rng.randrange(3))
            if rng.randrange(2):
                right[j][i] = monomial(rng.randrange(3), 0, rng.randrange(3))
    orders = sorted(rng.randrange(6) for _ in range(order))
    diagonal = [[monomial(orders[i], 0, 0) if i == j else ZERO
                 for j in range(order)] for i in range(order)]
    source = matrix_product(matrix_product(left, diagonal), right)
    degree = max(sum(term) for row in source for value in row for term in value)
    homogeneous = [[frozenset((a, degree - a - c, c) for a, _, c in value)
                    for value in row] for row in source]
    assert order * degree < 128
    return homogeneous, degree, orders


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='evidence/cube-generic-smith128.json')
    args = parser.parse_args()
    started = time.monotonic()
    rng = random.Random(2026110129)
    field = field_controls(rng)
    known = []
    for order in range(1, 7):
        for case in range(16):
            matrix, degree, expected = known_profile_matrix(order, rng)
            result = generic_smith_profile(matrix, 0, [degree] * order, [degree] * order)
            assert result['orders'] == expected
            known.append({'order': order, 'sample': case, 'degree': degree,
                          'orders': expected})
    frames = []
    source = kernel()
    for axis in range(3):
        rows, columns = axis_parts(axis)
        matrix = [[source[i][j] for j in columns] for i in rows]
        result = generic_smith_profile(matrix, axis, [6] * 4, [6] * 4)
        assert result['orders'] == [2, 4, 8, 10]
        frames.append(result)
    matrix = coupled_control_block(0)
    independent = complete_minor_profile(matrix, 0)
    coupled = generic_smith_profile(matrix, 0, [6, 18] * 4, [6, 18] * 4)
    assert coupled['orders'] == independent['orders']
    for axis in (1, 2):
        result = generic_smith_profile(coupled_control_block(axis), axis,
                                       [6, 18] * 4, [6, 18] * 4)
        assert result['orders'] == coupled['orders']
    misleading = [[monomial(1, 0, 0) ^ monomial(0, 1, 0) ^ monomial(0, 0, 1)]]
    assert generic_smith_profile(misleading, 0, [1], [1])['orders'] == [0]
    assert generic_smith_profile([[monomial(1, 0, 0)]], 0, [1], [1])['orders'] == [1]
    boundary = generic_smith_profile([[monomial(127, 0, 0)]], 0, [127], [127])
    assert boundary['orders'] == [127] and boundary['absolute_series_precision'] == 128
    rejections = 0
    for matrix, rows, columns in [([[monomial(128, 0, 0)]], [128], [128]),
                                  ([[monomial(1, 0, 0)]], [0], [0])]:
        try:
            generic_smith_profile(matrix, 0, rows, columns)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError('An unsupported generic certificate was accepted.')
    paths = [Path(__file__), Path(__file__).with_name('cube_generic_smith128.py'),
             Path(__file__).with_name('cube_generic_smith128.cpp'),
             Path(__file__).with_name('cube_homogeneous_gaussian.py')]
    report = {'schema': 'cube-generic-smith128-v1', 'field_controls': field,
              'known_profiles': known, 'scalar_three_cut_frames': frames,
              'coupled_generic_profile': coupled,
              'independent_complete_minors': independent['nonempty_minors_evaluated'],
              'specialization_cancellation_control': {'generic_order': 0, 'both_variables_one_order': 1},
              'last_supported_order': boundary, 'unsupported_input_rejections': rejections,
              'new_factoring_score': False,
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in paths},
              'elapsed_seconds': time.monotonic() - started}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'known_profiles': len(known),
                      'field_products': field['independent_field_products'],
                      'independent_minors': report['independent_complete_minors'],
                      'seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()
