#!/usr/bin/env python3
"""Exact polynomial blocks and right-inverse controls for channel grouping."""

from itertools import product

import check_cancellation_module_automata as arithmetic
import check_minimal_boundary_modules as linear
import check_exact_module_orbits as lattices


ZERO, ONE = arithmetic.ZERO, arithmetic.ONE


def polynomial_matrix_product(left, right):
    rows, inner, columns = len(left), len(right), len(right[0])
    assert all(len(row) == inner for row in left)
    answer = [[0] * columns for _ in range(rows)]
    for i, k, j in product(range(rows), range(inner), range(columns)):
        answer[i][j] ^= arithmetic.pmul(left[i][k], right[k][j])
    return answer


def block_local_crossing(local, alphabet, degree, length):
    """Group length successive sites; base-alphabet digits are low first."""
    if alphabet < 2 or degree < 0 or length < 1:
        raise ValueError("invalid polynomial block parameters")
    dimension = len(next(iter(local.values())))
    assert set(local) == set(product(range(alphabet), repeat=4))
    block_alphabet = alphabet**length
    digits = [tuple((value // alphabet**position) % alphabet for position in range(length))
              for value in range(block_alphabet)]
    answer = {}
    for labels in product(range(block_alphabet), repeat=4):
        matrix = [[int(i == j) for j in range(dimension)] for i in range(dimension)]
        for position in range(length):
            letter = tuple(digits[label][position] for label in labels)
            matrix = polynomial_matrix_product(matrix, local[letter])
        assert all(value.bit_length() <= length*degree+1 for row in matrix for value in row)
        answer[labels] = matrix
    return answer


def transition_product(left, right):
    answer = []
    for row in left:
        result = {}
        for middle, first in row.items():
            for target, second in right[middle].items():
                result[target] = result.get(target, 0) + first*second
        answer.append(result)
    return answer


def regular_right_inverse(matrix):
    """A regular right inverse from a unit residue-column minor."""
    rank, width = len(matrix), len(matrix[0])
    assert all(len(row) == width for row in matrix)
    pivots, selected = {}, []
    for column in range(width):
        value = sum(arithmetic.jet(matrix[row][column], 1) << row for row in range(rank))
        for pivot in sorted(pivots, reverse=True):
            if value >> pivot & 1:
                value ^= pivots[pivot]
        if value:
            pivots[value.bit_length()-1] = value
            selected.append(column)
            if len(selected) == rank:
                break
    assert len(selected) == rank, "right factor is not surjective over the regular ring"
    minor = tuple(tuple(row[column] for column in selected) for row in matrix)
    inverse = linear.inverse(minor)
    assert all(lattices.regular(value) for row in inverse for value in row)
    answer = [[ZERO] * rank for _ in range(width)]
    for index, column in enumerate(selected):
        answer[column] = list(inverse[index])
    product_matrix = linear.matrix_product(matrix, answer)
    assert product_matrix == tuple(tuple(ONE if i == j else ZERO for j in range(rank))
                                   for i in range(rank))
    return tuple(tuple(row) for row in answer)


def conditional_right_product(first, seconds):
    """Rows (first channel, next channel), columns (first word, next word).

    Next-state channel counts may differ. The output is rectangular and
    retains a regular right inverse whenever every input does.
    """
    first_rank, width = len(first), len(first[0])
    assert len(seconds) == first_rank
    assert all(len(row) == width for matrix in seconds for row in matrix)
    answer = []
    for channel, second in enumerate(seconds):
        for row in second:
            answer.append(tuple(arithmetic.rmul(first[channel][a], row[b])
                                for a in range(width) for b in range(width)))
    return tuple(answer)
