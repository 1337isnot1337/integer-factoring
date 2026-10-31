#!/usr/bin/env python3
"""Exact finite lattice certificate for an eight-wire Cauchy packet.

All matrix operations use Python integers modulo a fixed power of two.
The degree bound supplies the exact rational rank; no floating rank is used.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path

from cauchy_packet_polynomials import cauchy_polynomials, shift


LENGTH = 8
DEGREE = LENGTH * (2 * LENGTH - 1)
PRECISION = 32
MODULUS = 1 << PRECISION
MASK = MODULUS - 1


def residue(value):
    value = Fraction(value)
    if not value.denominator & 1:
        raise ValueError('The local coefficient is not two-adically integral.')
    return value.numerator * pow(value.denominator, -1, MODULUS) & MASK


def multiply_truncated(left, right, order=LENGTH):
    result = [0] * order
    for i, a in enumerate(left[:order]):
        if a:
            for j, b in enumerate(right[:order - i]):
                if b:
                    result[i + j] += a * b
    return result


def packet_jet_matrix():
    global_kernels = cauchy_polynomials(LENGTH)
    rows = []
    for center in range(-LENGTH, LENGTH):
        kernels = {difference: [residue(value) for value in shift(polynomial, center, 2 * LENGTH)[:LENGTH]]
                   for difference, polynomial in global_kernels.items()}
        columns = []
        for state in range(1 << LENGTH):
            polynomial = [1] + [0] * (LENGTH - 1)
            for x in range(LENGTH):
                factor = [-sum(kernels[y - x][degree] for y in range(LENGTH) if state >> y & 1)
                          for degree in range(LENGTH)]
                factor[0] += 1
                polynomial = [value & MASK for value in multiply_truncated(polynomial, factor)]
            columns.append(polynomial)
        rows.extend([[columns[state][degree] for state in range(1 << LENGTH)]
                     for degree in range(LENGTH)])
    return rows


def tapered_jet(center, retained_levels):
    h = LENGTH // retained_levels
    polynomial = [1] + [0] * (LENGTH - 1)
    denominator = 1
    degree = 0
    for root in range(-LENGTH, LENGTH):
        if not root:
            continue
        absolute = abs(root)
        valuation = (absolute & -absolute).bit_length() - 1
        multiplicity = h * max(retained_levels - valuation, 0)
        for _ in range(multiplicity):
            polynomial = multiply_truncated(polynomial, [root - center, -2 * LENGTH])
            denominator *= root
            degree += 1
    assert degree <= DEGREE
    return [residue(Fraction(value, denominator)) for value in polynomial], degree


def normalized_newton_jet(center):
    denominator_exponent = sum(DEGREE // (1 << level) for level in range(1, 5))
    polynomial = [1] + [0] * (LENGTH - 1)
    for root in range(DEGREE):
        polynomial = multiply_truncated(polynomial, [center - root, 2 * LENGTH])
    divisor = 1 << denominator_exponent
    assert all(value % divisor == 0 for value in polynomial)
    return [(value // divisor) & MASK for value in polynomial]


def right_sides():
    rows = []
    for center in range(-LENGTH, LENGTH):
        tapered = [tapered_jet(center, levels)[0] for levels in (1, 2, 3)]
        negative = normalized_newton_jet(center)
        for order in range(LENGTH):
            monomials = [(comb(degree, order) * center**(degree - order) * (2 * LENGTH)**order) & MASK
                         if order <= degree else 0 for degree in range(DEGREE + 1)]
            rows.append(monomials + [values[order] for values in tapered]
                        + [negative[order], (negative[order] << 26) & MASK])
    return rows


def smith_membership(matrix, right):
    matrix = [row[:] for row in matrix]
    right = [row[:] for row in right]
    row_count, column_count = len(matrix), len(matrix[0])
    targets = len(right[0])
    first_failure = [None] * targets
    pivots = []
    for index in range(min(row_count, column_count)):
        best = None
        for i in range(index, row_count):
            for j in range(index, column_count):
                value = matrix[i][j]
                if value:
                    order = (value & -value).bit_length() - 1
                    if best is None or order < best[0]:
                        best = order, i, j
            if best is not None and best[0] == 0:
                break
        if best is None:
            break
        order, i, j = best
        matrix[index], matrix[i] = matrix[i], matrix[index]
        right[index], right[i] = right[i], right[index]
        for row in matrix:
            row[index], row[j] = row[j], row[index]
        inverse = pow(matrix[index][index] >> order, -1, MODULUS)
        matrix[index] = [(value * inverse) & MASK for value in matrix[index]]
        right[index] = [(value * inverse) & MASK for value in right[index]]
        assert matrix[index][index] == 1 << order
        divisor_mask = (1 << order) - 1
        for target, value in enumerate(right[index]):
            if value & divisor_mask and first_failure[target] is None:
                first_failure[target] = {'pivot': index, 'pivot_valuation': order,
                                         'right_valuation': (value & -value).bit_length() - 1}
        pivot_row = matrix[index]
        pivot_right = right[index]
        for i in range(index + 1, row_count):
            factor = matrix[i][index] >> order
            if factor:
                row = matrix[i]
                row[index:] = [(value - factor * pivot) & MASK
                               for value, pivot in zip(row[index:], pivot_row[index:])]
                right[i] = [(value - factor * pivot) & MASK
                            for value, pivot in zip(right[i], pivot_right)]
            assert matrix[i][index] == 0
        # The corresponding column operations affect no remaining row,
        # because the pivot column has just been cleared below the pivot.
        matrix[index][index + 1:] = [0] * (column_count - index - 1)
        pivots.append(order)

    rank = len(pivots)
    tails = [sum(bool(right[row][target]) for row in range(rank, row_count))
             for target in range(targets)]
    return pivots, first_failure, tails


def check():
    matrix = packet_jet_matrix()
    right = right_sides()
    pivots, failures, tails = smith_membership(matrix, right)
    assert len(pivots) == DEGREE + 1 == 121
    assert max(pivots) == 26 < PRECISION
    positive_count = DEGREE + 1 + 3
    assert all(failure is None for failure in failures[:positive_count])
    assert not any(tails[:positive_count])
    assert failures[positive_count] is not None
    assert failures[positive_count + 1] is None and tails[positive_count + 1] == 0

    # The ordinary leading coefficient independently forbids the negative
    # control: packet columns have denominator valuation at most 8*11=88,
    # whereas the normalized degree-120 Newton polynomial has value -112.
    derivative_valuation = sum((abs(root) & -abs(root)).bit_length() - 1
                               for root in range(-LENGTH, LENGTH) if root)
    newton_exponent = sum(DEGREE // (1 << level) for level in range(1, 5))
    assert LENGTH * derivative_valuation == 88 < newton_exponent == 112

    root = Path(__file__).resolve().parent.parent
    sources = ['tools/check_cauchy_packet_lattice.py', 'tools/cauchy_packet_polynomials.py',
               'research/eight-wire-cauchy-lattice.md']
    return {'status': 'pass', 'L': LENGTH, 'fixed_row_state': 0,
            'polynomial_degree_limit': DEGREE, 'matrix_shape': [len(matrix), len(matrix[0])],
            'jet_order_per_center': LENGTH, 'modular_precision': PRECISION,
            'proved_rational_rank': len(pivots), 'maximum_pivot_valuation': max(pivots),
            'pivot_valuation_multiplicities': dict(sorted(Counter(pivots).items())),
            'integral_monomials_realized': DEGREE + 1,
            'wire_scale_tapered_witnesses_realized': [
                {'M': LENGTH, 'u': levels, 'degree': tapered_jet(0, levels)[1]}
                for levels in (1, 2, 3)],
            'proper_lattice_negative_control': failures[positive_count],
            'negative_control_leading_denominator_exponent': newton_exponent,
            'maximum_packet_leading_denominator_exponent': LENGTH * derivative_valuation,
            'scaled_negative_control_realized': 26,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in sources},
            'scope': 'An exact finite eight-wire realization theorem over rationals with '
                     'odd denominators. It provides no all-size construction, no canonical '
                     'low-precision realization, and no factoring exponent improvement.'}


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
