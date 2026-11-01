#!/usr/bin/env python3
"""Exact dense controls for the global pole trace's Cauchy determinant slope."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import random

from interval_pole_trace import (factorial_vandermonde_content,
                                 interval_vandermonde_content, square_pole_domains)


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def multiply(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def inverse(x):
    a, b = x
    norm = a * a - a * b + b * b
    assert norm
    return (a - b) / norm, -b / norm


def determinant(matrix):
    matrix = [row[:] for row in matrix]
    result = (F(1), F(0))
    for i in range(len(matrix)):
        pivot = next(j for j in range(i, len(matrix)) if matrix[j][i] != (0, 0))
        if pivot != i:
            matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
            result = (-result[0], -result[1])
        result = multiply(result, matrix[i][i])
        reciprocal = inverse(matrix[i][i])
        for j in range(i + 1, len(matrix)):
            ratio = multiply(matrix[j][i], reciprocal)
            for k in range(i + 1, len(matrix)):
                term = multiply(ratio, matrix[i][k])
                matrix[j][k] = add(matrix[j][k], (-term[0], -term[1]))
            matrix[j][i] = (F(0), F(0))
    return result


def solve(matrix, right):
    matrix = [row[:] for row in matrix]
    right = [row[:] for row in right]
    size = len(matrix)
    for i in range(size):
        pivot = next(j for j in range(i, size) if matrix[j][i] != (0, 0))
        matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
        right[i], right[pivot] = right[pivot], right[i]
        reciprocal = inverse(matrix[i][i])
        matrix[i] = [multiply(entry, reciprocal) for entry in matrix[i]]
        right[i] = [multiply(entry, reciprocal) for entry in right[i]]
        for j in range(size):
            if j == i:
                continue
            factor = matrix[j][i]
            for row, source in ((matrix[j], matrix[i]), (right[j], right[i])):
                for k in range(len(row)):
                    term = multiply(factor, source[k])
                    row[k] = add(row[k], (-term[0], -term[1]))
    return right


def v2(value):
    value = F(value)
    assert value
    numerator, denominator = abs(value.numerator), value.denominator
    return ((numerator & -numerator).bit_length()
            - (denominator & -denominator).bit_length())


def valuation(value):
    return min(v2(coefficient) for coefficient in value if coefficient)


def check():
    rng = random.Random(20261040)
    controls = positive = flow_controls = trace_controls = max_pairs = empty_controls = 0
    block_controls, negative_schur_example = 0, None
    examples = []
    for index in range(48):
        left_size, right_size = rng.randrange(1, 6), rng.randrange(1, 6)
        lower, other = rng.randrange(1, 12), rng.randrange(1, 12)
        original_left = range(lower, lower + left_size)
        original_right = range(other, other + right_size)
        number = 2 * rng.randrange(1, 256) + 1
        if index % 2 == 0:
            # Force a genuine pole, including multiple-pair and square cases.
            odd_left = [a for a in original_left if a % 2]
            odd_right = [b for b in original_right if b % 2]
            if odd_left and odd_right:
                number = rng.choice(odd_left) * rng.choice(odd_right)
        if index == 0:
            number, original_left, original_right = 15, range(3, 6), range(3, 6)
        elif index == 1:
            number, original_left, original_right = 3, range(9, 12), range(7, 10)
        elif index == 2:
            number, original_left, original_right = 9, range(1, 4), range(1, 4)
        elif index == 4:
            number, original_left, original_right = 45, range(5, 10), range(5, 10)
        domains = square_pole_domains(number, (original_left.start, original_left.stop),
                                      (original_right.start, original_right.stop))
        left, right = [[x for interval in domain for x in range(*interval)] for domain in domains]
        size = len(left)
        assert size == len(right)
        pairs = [(a, b) for a in original_left for b in original_right if a * b == number]
        max_pairs = max(max_pairs, len(pairs))
        empty_controls += size == 0
        assert pairs == [(a, b) for a in left for b in right if a * b == number]
        bound = max((abs(number - a * b) for a in left for b in right), default=0).bit_length()
        scale_bits = bound + 3
        orders, schur_orders, previous_matrix = [], [], None
        for k in (scale_bits, scale_bits + 1):
            z = (F(number), F(1 << k))
            matrix = [[inverse((F(number - a * b), F(1 << k))) for b in right]
                      for a in left]
            direct = determinant(matrix)
            # The full Cauchy determinant formula, not just its valuation.
            formula = (F(1), F(0))
            for _ in range(size * (size - 1) // 2):
                formula = multiply(formula, z)
            factor = 1
            for i, a in enumerate(left):
                for b in list(left)[i + 1:]:
                    factor *= b - a
            for i, a in enumerate(right):
                for b in list(right)[i + 1:]:
                    factor *= b - a
            formula = multiply(formula, (F(factor), F(0)))
            for row in matrix:
                for entry in row:
                    formula = multiply(formula, entry)
            assert direct == formula
            nonzero_content = sum(v2(number - a * b) for a in left for b in right
                                  if a * b != number)
            known_content = sum(interval_vandermonde_content(domain) for domain in domains)
            expected = known_content - nonzero_content - k * len(pairs)
            assert valuation(direct) == expected
            orders.append(size * k + valuation(direct))
            scaled = [[multiply(entry, (F(0), F(1 << k))) for entry in row] for row in matrix]
            if previous_matrix is not None:
                for i, row in enumerate(previous_matrix):
                    for j, entry in enumerate(row):
                        flowed = multiply((F(2), F(0)),
                                          multiply(entry, inverse(add((F(1), F(0)), entry))))
                        assert scaled[i][j] == flowed
                        flow_controls += 1
            previous_matrix = scaled
            if size >= 2:
                cut = size - 1 if index == 0 else max(1, size // 2)
                block = [row[:cut] for row in scaled[:cut]]
                upper_right = [row[cut:] for row in scaled[:cut]]
                solved = solve(block, upper_right)
                schur = []
                for i in range(cut, size):
                    row = []
                    for j in range(cut, size):
                        entry = scaled[i][j]
                        for r in range(cut):
                            term = multiply(scaled[i][r], solved[r][j - cut])
                            entry = add(entry, (-term[0], -term[1]))
                        expected_entry = scaled[i][j]
                        for r in range(cut):
                            expected_entry = multiply(expected_entry, z)
                            expected_entry = multiply(expected_entry,
                                (F((left[i] - left[r]) * (right[j] - right[r])), F(0)))
                            expected_entry = multiply(expected_entry, matrix[i][r])
                            expected_entry = multiply(expected_entry, matrix[r][j])
                        assert entry == expected_entry
                        row.append(entry)
                        block_controls += 1
                    schur.append(row)
                # determinant(scaled) is c^size times determinant(unscaled).
                scaled_det = direct
                for _ in range(size):
                    scaled_det = multiply(scaled_det, (F(0), F(1 << k)))
                assert multiply(determinant(block), determinant(schur)) == scaled_det
                schur_orders.append(valuation(determinant(schur)))
            quotient = solve(scaled, [[multiply(entry, entry) for entry in row] for row in scaled])
            scalar_trace = (F(0), F(0))
            total = (F(0), F(0))
            for i, row in enumerate(scaled):
                scalar_trace = add(scalar_trace, quotient[i][i])
                for entry in row:
                    total = add(total, entry)
            correction = multiply((F(0), F(size * (size - 1) // 2 * (1 << k))), inverse(z))
            assert scalar_trace == add(total, (-correction[0], -correction[1]))
            for coefficient, expected_coefficient in zip(scalar_trace, (len(pairs), 0)):
                assert coefficient.denominator % 2
                residue = coefficient.numerator * pow(coefficient.denominator, -1, 8) % 8
                assert residue == expected_coefficient % 8
            trace_controls += 1
            controls += 1
        assert orders[1] - orders[0] == size - len(pairs)
        if index == 0:
            assert schur_orders[1] - schur_orders[0] == -1
            negative_schur_example = {'number': number, 'nodes': left,
                                     'eliminated_prefix_size': size - 1,
                                     'scale_bits': scale_bits, 'schur_orders': schur_orders,
                                     'schur_order_slope': -1}
        positive += bool(pairs)
        if len(examples) < 3 and pairs:
            examples.append({'number': number, 'original_left': [original_left.start, original_left.stop],
                             'original_right': [original_right.start, original_right.stop],
                             'padded_domains': domains, 'pair_count': len(pairs),
                             'scale_bits': scale_bits, 'scaled_determinant_orders': orders})
    for length in range(1, 129):
        literal = sum(v2(j - i) for i in range(length) for j in range(i + 1, length))
        assert factorial_vandermonde_content(length) == literal
    for _ in range(128):
        lower, first_length, gap, second_length = [rng.randrange(1, 20) for _ in range(4)]
        intervals = ((lower, lower + first_length),
                     (lower + first_length + gap, lower + first_length + gap + second_length))
        nodes = [x for interval in intervals for x in range(*interval)]
        literal = sum(v2(b - a) for i, a in enumerate(nodes) for b in nodes[i + 1:])
        assert interval_vandermonde_content(intervals) == literal
    huge = []
    for bits in (80, 200):
        length = (1 << bits) + 12
        content = factorial_vandermonde_content(length)
        assert factorial_vandermonde_content(length + 1) - content == length - length.bit_count()
        huge.append({'length': length, 'vandermonde_content_bits': content.bit_length()})
    return {'scope': 'Exact square Cauchy determinant and valuation slope; no fast global valuation algorithm.',
            'seed': 20261040, 'dense_determinant_controls': controls,
            'slope_controls': controls // 2, 'positive_rectangles': positive,
            'maximum_factor_pair_count': max_pairs, 'empty_padded_domains': empty_controls,
            'literal_entry_flow_controls': flow_controls, 'normalized_trace_controls': trace_controls,
            'literal_schur_entry_controls': block_controls,
            'nonintegral_schur_control': negative_schur_example,
            'vandermonde_content_controls': 128, 'two_interval_content_controls': 128,
            'huge_known_content_controls': huge, 'examples': examples,
            'unresolved': 'Evaluate the determinant-content difference or normalized inverse/Hadamard-square trace without expanding the matrix.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
