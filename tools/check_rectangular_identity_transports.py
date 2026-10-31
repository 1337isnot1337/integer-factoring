#!/usr/bin/env python3
"""General constant core transports with exact identity ramified endpoint.

The validator checks C F(alpha,beta) W=I. It does not use a triangular
flag. Individual dressed maps may be rational; paired-type polynomiality
follows from the identity endpoint and the original height dressing.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random

from check_rectangular_f4_repair import (
    ZERO, add, scale, times_binary, binary_power, matrix_multiply,
    residual, full_inverse_repair, partial_repair, build, fmul, finv,
)
from check_rectangular_refined_repairs import height, _verified_repair


IDENTITY = [[int(i == j) for j in range(12)] for i in range(12)]


def checked_matrix(matrix):
    if len(matrix) != 12 or any(len(row) != 12 for row in matrix):
        raise ValueError('A constant core map must be twelve by twelve.')
    if any(value not in (0, 1, 2, 3) for row in matrix for value in row):
        raise ValueError('Constant core coefficients must lie in F4.')
    return [list(row) for row in matrix]


def inverse_matrix(matrix):
    matrix = checked_matrix(matrix)
    augmented = [row+IDENTITY[i][:] for i, row in enumerate(matrix)]
    for column in range(12):
        pivot = next((row for row in range(column, 12) if augmented[row][column]), None)
        if pivot is None:
            raise ValueError('The proposed constant core map is singular.')
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scalar = finv(augmented[column][column])
        augmented[column] = [fmul(scalar, value) for value in augmented[column]]
        for row in range(12):
            if row != column and augmented[row][column]:
                scalar = augmented[row][column]
                augmented[row] = [value ^ fmul(scalar, source)
                                  for value, source in zip(augmented[row], augmented[column])]
    inverse = [row[12:] for row in augmented]
    assert matrix_multiply(matrix, inverse) == IDENTITY
    assert matrix_multiply(inverse, matrix) == IDENTITY
    return inverse


def validate_identity_transport(alpha, beta, left, right):
    """Validate C F W=I and return independent copies of the two matrices."""
    return validate_typed_transport(alpha, beta, left, right, IDENTITY)


def validate_typed_transport(alpha, beta, left, right, target):
    """Validate a height-unitriangular target and the exact identity C F W=target."""
    if alpha not in (0, 1, 2, 3) or beta not in (0, 1, 2, 3):
        raise ValueError('The code parameters must lie in F4.')
    left, right, target = checked_matrix(left), checked_matrix(right), checked_matrix(target)
    for i, row in enumerate(target):
        if row[i] != 1 or any(value and i != j and height(i) >= height(j) for j, value in enumerate(row)):
            raise ValueError('The target must have diagonal one and only strict actual-height-lowering extras.')
    endpoint = matrix_multiply(matrix_multiply(left, residual(alpha, beta)), right)
    if endpoint != target:
        raise ValueError('The constant core maps do not satisfy C F(alpha,beta) W=target.')
    return left, right


def conjugation_spec(alpha, beta, conjugator):
    """Return C=V F^-1 and W=V^-1 for any invertible constant V."""
    conjugator = checked_matrix(conjugator)
    left = matrix_multiply(conjugator, full_inverse_repair(alpha, beta))
    right = inverse_matrix(conjugator)
    return validate_identity_transport(alpha, beta, left, right)


def typed_conjugation_spec(alpha, beta, conjugator):
    """Return C,W,U preserving the original partial-repair endpoint U."""
    conjugator = checked_matrix(conjugator)
    base = partial_repair(alpha, beta)
    target = matrix_multiply(base, residual(alpha, beta))
    left = matrix_multiply(conjugator, base)
    right = matrix_multiply(matrix_multiply(inverse_matrix(target), inverse_matrix(conjugator)), target)
    validate_typed_transport(alpha, beta, left, right, target)
    return left, right, target


def conjugate_seed(seed, conjugator, check_types=True):
    """Dress an independently built full-inverse seed by V and V^-1."""
    if seed['repair_kind'] != 'mixed_height_inverse' or seed['ramified_endpoint_matrix'] != IDENTITY:
        raise ValueError('Start from the independently constructed full-inverse seed.')
    conjugator = checked_matrix(conjugator)
    inverse = inverse_matrix(conjugator)
    return _apply_transport_seed(seed, conjugator, inverse, IDENTITY, check_types)


def typed_conjugate_seed(seed, conjugator, check_types=True):
    """Preserve a partial seed's original endpoint while conjugating its core."""
    if seed['repair_kind'] != 'mixed_height_partial':
        raise ValueError('Start from an independently constructed partial-repair seed.')
    left, right, target = typed_conjugation_spec(seed['alpha'], seed['beta'], conjugator)
    return _apply_transport_seed(seed, conjugator, right, target, check_types)


def _apply_transport_seed(seed, left_change, right_change, target, check_types):
    conjugator, inverse = checked_matrix(left_change), checked_matrix(right_change)
    left = matrix_multiply(conjugator, seed['repair_matrix'])
    right = inverse
    validate_typed_transport(seed['alpha'], seed['beta'], left, right, target)
    A, B, bulk = seed['q_x'], seed['q_u'], seed['bulk_size']
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    left_gate = [[(source % 4, source // 4, coefficient,
                   binary_power(3, 2+height(source)-height(target)))
                  for source, coefficient in enumerate(row) if coefficient]
                 for target, row in enumerate(conjugator)]
    middle = [ZERO]*len(seed['values'])
    for x, y, u, v in product_letters(A, B):
        answer = ZERO
        for vs, ys, coefficient, polynomial in left_gate[v//bulk+4*(y//bulk)]:
            source = seed['values'][position(x, ys*bulk+y % bulk, u, vs*bulk+v % bulk)]
            answer = add(answer, scale(times_binary(source, polynomial), coefficient))
        middle[position(x, y, u, v)] = answer
    right_gate = [[(source % 4, source // 4, inverse[source][target],
                    binary_power(3, 2+height(target)-height(source)))
                   for source in range(12) if inverse[source][target]] for target in range(12)]
    values = [ZERO]*len(middle)
    for x, y, u, v in product_letters(A, B):
        answer = ZERO
        for xs, us, coefficient, polynomial in right_gate[x//bulk+4*(u//bulk)]:
            source = middle[position(xs*bulk+x % bulk, y, us*bulk+u % bulk, v)]
            answer = add(answer, scale(times_binary(source, polynomial), coefficient))
        values[position(x, y, u, v)] = answer
    uncancelled_degree = seed['denominator_degree']+4
    provisional = dict(seed, denominator_degree=uncancelled_degree)
    result = _verified_repair(provisional, values, target, [(height(i), 0) for i in range(12)], check_types)
    assert result['denominator_degree'] <= 6
    result.update({'repair_kind': 'mixed_typed_transport', 'repair_matrix': left,
                   'input_correction': right, 'constant_left_change': conjugator,
                   'uncancelled_denominator_degree': uncancelled_degree,
                   'cancelled_c_order': uncancelled_degree-result['denominator_degree']})
    return result


def product_letters(A, B):
    for x in range(A):
        for y in range(B):
            for u in range(B):
                for v in range(A):
                    yield x, y, u, v


def build_identity_transport(left, right, indices=(0, 1, 2), ambient=4,
                             alpha=1, beta=2, check_types=True, quartic_terms=None):
    """Independent polynomial-numerator constructor for a validated C,F,W triple."""
    left, right = validate_identity_transport(alpha, beta, left, right)
    conjugator = matrix_multiply(left, residual(alpha, beta))
    assert inverse_matrix(conjugator) == right
    base = build(indices, ambient, alpha, beta, check_types=False,
                 independent_exterior=len(indices) <= 3, quartic_terms=quartic_terms, repair='inverse')
    return conjugate_seed(base, conjugator, check_types)


def build_typed_transport(left, right, target, indices=(0, 1, 2), ambient=4,
                          alpha=1, beta=2, check_types=True, quartic_terms=None):
    """Independent numerator constructor for any validated height-unitriangular target."""
    left, right = validate_typed_transport(alpha, beta, left, right, target)
    base = build(indices, ambient, alpha, beta, check_types=False,
                 independent_exterior=len(indices) <= 3, quartic_terms=quartic_terms, repair='inverse')
    left_change = matrix_multiply(left, residual(alpha, beta))
    return _apply_transport_seed(base, left_change, right, target, check_types)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = build((0, 1), 2, 1, 1, repair='inverse')
    cases = []
    for row in range(12):
        for column in range(12):
            if row == column:
                continue
            for coefficient in (1, 2, 3):
                V = [line[:] for line in IDENTITY]
                V[row][column] = coefficient
                seed = conjugate_seed(base, V)
                cases.append({'row': row, 'column': column, 'coefficient': coefficient,
                              'actual_height_change': height(row)-height(column),
                              'denominator_degree': seed['denominator_degree'],
                              'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                              'two_copy_type_checks': seed['two_copy_type_checks'],
                              'coefficient_sha256': seed['coefficient_sha256']})
    dense = []
    rng = random.Random(1930107)
    for alpha, beta in ((1, 1), (1, 2), (0, 0)):
        V = [line[:] for line in IDENTITY]
        for _ in range(30):
            row, column = rng.sample(range(12), 2)
            coefficient = rng.randrange(1, 4)
            V[row] = [a ^ fmul(coefficient, b) for a, b in zip(V[row], V[column])]
        V[0] = [fmul(2, value) for value in V[0]]
        left, right = conjugation_spec(alpha, beta, V)
        seed = build_identity_transport(left, right, alpha=alpha, beta=beta)
        dense.append({key: value for key, value in seed.items() if key not in ('values', 'residual_matrix', 'c_order_histogram')})
    partial_base = build((0, 1), 2, 1, 2)
    typed = []
    for row in range(12):
        for column in range(12):
            if row == column:
                continue
            V = [line[:] for line in IDENTITY]
            V[row][column] = 2
            seed = typed_conjugate_seed(partial_base, V)
            typed.append({'row': row, 'column': column, 'coefficient': 2,
                          'actual_height_change': height(row)-height(column),
                          'denominator_degree': seed['denominator_degree'],
                          'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                          'two_copy_type_checks': seed['two_copy_type_checks'],
                          'coefficient_sha256': seed['coefficient_sha256']})
    try:
        validate_identity_transport(1, 1, IDENTITY, IDENTITY)
    except ValueError:
        rejected = True
    else:
        raise AssertionError('An incorrect identity transport was accepted.')
    result = {'status': 'exact arbitrary constant core transports with certified typed endpoints; no pressure claim',
              'coupled_transvection_cases': cases, 'dense_transport_cases': dense,
              'partial_target_transvection_cases': typed,
              'incorrect_endpoint_rejected': rejected}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
