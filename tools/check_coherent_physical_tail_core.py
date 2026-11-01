#!/usr/bin/env python3
"""Independent full-kernel controls for the eight-k physical tail core."""

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from coherent_physical_tail_core import concentration_numerator, polynomial_images, finite_field_images
from coherent_physical_tail_low_band import MUL, F4, mul, power
from reflection_dual_rail_gaussian import gauged_kernel


def evaluate(polynomial, t):
    value = 0
    for bit in range(polynomial.bit_length() - 1, -1, -1):
        value = mul(value, t) ^ ((polynomial >> bit) & 1)
    return value


def matrix_product(A, B):
    A, B = np.asarray(A, dtype=np.uint8), np.asarray(B, dtype=np.uint8)
    return np.bitwise_xor.reduce(MUL[A[:, :, None], B[None, :, :]], axis=1)


def inverse(A):
    n = len(A)
    augmented = np.concatenate((np.asarray(A, dtype=np.uint8), np.eye(n, dtype=np.uint8)), axis=1)
    for column in range(n):
        choices = np.flatnonzero(augmented[column:, column])
        assert len(choices)
        pivot = column + int(choices[0])
        augmented[[column, pivot]] = augmented[[pivot, column]]
        augmented[column] = MUL[power(int(augmented[column, column]), -1), augmented[column]]
        rows = np.flatnonzero(augmented[:, column])
        rows = rows[rows != column]
        augmented[rows] ^= MUL[augmented[rows, column, None], augmented[column, None, :]]
    assert np.array_equal(augmented[:, :n], np.eye(n, dtype=np.uint8))
    return augmented[:, n:]


def double(A):
    n = len(A)
    answer = np.zeros((2 * n, 2 * n), dtype=np.uint8)
    answer[:n, :n] = answer[n:, n:] = A
    return answer


def half(polynomial):
    assert all(not (polynomial >> j & 1) for j in range(1, polynomial.bit_length(), 2))
    return sum(((polynomial >> j) & 1) << (j // 2) for j in range(0, polynomial.bit_length(), 2))


def literal_kernel(M):
    source = gauged_kernel(2 * M, 'zero')
    # Ascending physical heights; the last two groups are swapped.
    order = [leg * M + M - 1 - a for leg in (0, 1, 3, 2) for a in range(M)]
    return [[half(source[i][j]) for j in order] for i in order]


def linear_data(M, k, t):
    n = 2 * M
    K = np.array([[evaluate(v, t) for v in row] for row in literal_kernel(M)], dtype=np.uint8)
    A, Fkernel = K[:n, :n], K[n:, n:]
    F = np.array([[evaluate(3 << (2 * (a & ~b)), t) for b in range(M)]
                  for a in range(M)], dtype=np.uint8)
    L, R = np.eye(n, dtype=np.uint8), np.eye(n, dtype=np.uint8)
    L[M:, :M] = F
    R[:M, M:] = F.T
    numerator, denominator = concentration_numerator(M, k)
    C = np.array([[mul(evaluate(v, t), power(evaluate(denominator, t), -1)) for v in row]
                  for row in numerator], dtype=np.uint8)
    U = inverse(C)
    J = np.array([power(t, 2 * a + 1) for a in range(M)], dtype=np.uint8)
    J_inverse = np.array([power(int(a), -1) for a in J], dtype=np.uint8)
    V = MUL[MUL[J_inverse[:, None], C], J[None, :]]
    Vinverse = MUL[MUL[J_inverse[:, None], U], J[None, :]]
    diagonal = np.diag(np.array([power(t, 2 * a + 1) for a in range(M)]
                     + [power(t, 2 * M + 2 * a + 1) for a in range(M)], dtype=np.uint8))
    assert np.array_equal(matrix_product(matrix_product(L, diagonal), R), K[:n, n:])
    assert np.array_equal(matrix_product(C, U), np.eye(M, dtype=np.uint8))
    assert np.array_equal(matrix_product(U, C), np.eye(M, dtype=np.uint8))
    return A, Fkernel, L, R, double(C), double(U), double(V), double(Vinverse)


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit in pivots:
                row ^= pivots[bit]
            else:
                pivots[bit] = row
                break
    return len(pivots)


def concentration_controls():
    counts = dict(weighted_polynomial_entries=0, constant_unit_ranks=0,
                  full_operator_entries=0, bivariate_specialization_entries=0,
                  literal_internal_kernel_entries=0)
    for k, M in ((1, 4), (1, 8), (2, 8), (2, 16), (2, 32), (4, 16), (4, 32)):
        C, denominator = concentration_numerator(M, k)
        assert denominator & 1
        weighted = []
        for i in range(M):
            row = []
            for j in range(M):
                value, shift = C[i][j], 2 * (j - i)
                if shift < 0:
                    assert not value & ((1 << -shift) - 1)
                    value >>= -shift
                else:
                    value <<= shift
                row.append(value)
                counts['weighted_polynomial_entries'] += 1
            weighted.append(row)
        for matrix in (C, weighted):
            assert binary_rank([sum((v & 1) << j for j, v in enumerate(row)) for row in matrix]) == M
            counts['constant_unit_ranks'] += 1
        K = literal_kernel(M)
        T = 1 << (2 * M)
        for a, b in product(range(2 * M), repeat=2):
            expected = ((1 ^ T) if a % M + b % M == M - 1 else 0) if a // M == b // M \
                       else 3 << (2 * ((a % M) & (b % M)))
            assert K[a][b] == K[2 * M + a][2 * M + b] == expected
            counts['literal_internal_kernel_entries'] += 2
        complete, orders, spectators = polynomial_images(M, k)
        assert sorted(orders + spectators) == list(range(1, 4 * M, 2))
        retained = list(range(2 * k)) + list(range(M - 2 * k, M))
        active = retained + [M + j for j in retained]
        for t in (2, 3):
            expected = finite_field_images(k, t, power(t, 2 * M), mul, power)
            A, Fkernel, L, R, C, U, V, Vinverse = linear_data(M, k, t)
            for group in range(4):
                for a in range(k):
                    physical = a + (M if group % 2 else 0)
                    unit = np.eye(2 * M, dtype=np.uint8)[:, physical]
                    zero = np.zeros(2 * M, dtype=np.uint8)
                    pre = ((unit, zero), (A[:, physical], unit)) if group < 2 \
                          else ((unit, Fkernel[physical]), (zero, unit))
                    for which, (creator, annihilator) in enumerate(pre):
                        if group < 2:
                            creator = matrix_product(U, matrix_product(L, creator[:, None]))[:, 0]
                            annihilator = matrix_product(matrix_product(annihilator[None, :], L), C)[0]
                        else:
                            creator = matrix_product(Vinverse, matrix_product(R, creator[:, None]))[:, 0]
                            annihilator = matrix_product(matrix_product(annihilator[None, :], R), V)[0]
                        for part, actual in enumerate((creator, annihilator)):
                            target = np.zeros(2 * M, dtype=np.uint8)
                            target[active] = expected[group][a][which][part]
                            assert np.array_equal(actual, target), (M, k, t, group, a, which, part)
                            polynomial = complete[group][a][which][part]
                            assert [evaluate(v, t) for v in polynomial] == expected[group][a][which][part]
                            counts['full_operator_entries'] += 2 * M
                            counts['bivariate_specialization_entries'] += 8 * k
    return counts


def action(matrix, operator, right=False):
    if right:
        return action(matrix.T, (operator[1], operator[0])).T
    result = np.zeros_like(matrix)
    size = len(matrix)
    for mode, (creator, annihilator) in enumerate(zip(*operator)):
        low = np.array([i for i in range(size) if not i >> mode & 1])
        high = low ^ (1 << mode)
        if creator:
            result[high] ^= MUL[creator, matrix[low]]
        if annihilator:
            result[low] ^= MUL[annihilator, matrix[high]]
    return result


def physical_operator(matrix, coefficients, images, right=False):
    creator, annihilator = images
    number = action(action(matrix, creator, True), annihilator, True) if right \
             else action(action(matrix, annihilator), creator)
    answer = MUL[int(coefficients[0][0]), matrix]
    answer ^= MUL[int(coefficients[0][0]) ^ int(coefficients[1][1]), number]
    answer ^= MUL[int(coefficients[0][1]), action(matrix, annihilator, right)]
    answer ^= MUL[int(coefficients[1][0]), action(matrix, creator, right)]
    return answer


def compound(T):
    n = len(T)
    size = 1 << n
    output = np.zeros((size, size), dtype=np.uint8)
    output[0, 0] = 1
    for state in range(1, size):
        bit = state & -state
        source = output[:, state ^ bit]
        column = bit.bit_length() - 1
        for mode in range(n):
            coefficient = int(T[mode, column])
            if coefficient:
                low = np.array([i for i in range(size) if not i >> mode & 1])
                output[low ^ (1 << mode), state] ^= MUL[coefficient, source[low]]
    return output


def gaussian_multiplier(A):
    n = len(A)
    output = np.eye(1 << n, dtype=np.uint8)
    zero = [0] * n
    for i in range(n):
        first = zero.copy()
        first[i] = 1
        for j in range(i + 1, n):
            if A[i, j]:
                second = zero.copy()
                second[j] = 1
                output ^= MUL[int(A[i, j]), action(action(output, (second, zero)), (first, zero))]
    return output


def literal_control():
    """All matrix entries, with a coherent nonhomogeneous one-bit code.

    At infinity, the ambient height corrections occur on row u and
    column v.  They must not be assigned to both row groups.
    """
    from eight_letter_reflection_retraction import ambient_normalized_tensor
    M, k = 4, 1
    tensor = ambient_normalized_tensor()
    reverse = lambda word: sum(((word >> i) & 1) << (3 - i) for i in range(4))
    P = np.array([1, F4[2]], dtype=np.uint8)
    Lmap = np.array([F4[2], F4[2]], dtype=np.uint8)
    assert int(np.bitwise_xor.reduce(MUL[P, Lmap])) == 1
    Pi = MUL[P[:, None], Lmap[None, :]]
    checked = 0
    for t in (2, 3):
        cache = {value: evaluate(value, t) for value in set(tensor)}
        ambient = np.array([[cache[tensor[reverse(row & 15) | reverse(row >> 4) << 4
                             | reverse(column >> 4) << 8 | reverse(column & 15) << 12]]
                             for column in range(256)] for row in range(256)], dtype=np.uint8)
        A, Fkernel, L, R, C, U, V, Vinverse = linear_data(M, k, t)
        GA, GF = gaussian_multiplier(A), gaussian_multiplier(Fkernel)
        left = matrix_product(matrix_product(compound(U), compound(L)), GA)
        right = matrix_product(matrix_product(GF.T, compound(R)), compound(V))
        groups = finite_field_images(k, t, power(t, 2 * M), mul, power)
        _, orders, unused = polynomial_images(M, k)
        diagonal = np.diag(np.array([power(t, sum(d for i, d in enumerate(orders) if state >> i & 1))
                           for state in range(256)], dtype=np.uint8))
        direct_images = []
        for mode in (0, M):
            creator, annihilator = [0] * 8, [0] * 8
            creator[mode] = annihilator[mode] = 1
            direct_images.append(((creator, [0] * 8), ([0] * 8, annihilator)))
        D = np.diag(np.array([1, power(t, -1)], dtype=np.uint8))
        for place in ('zero', 'infinity'):
            coefficients = [matrix_product(Pi, D), matrix_product(Pi, D), Pi, Pi] if place == 'zero' \
                else [Pi, matrix_product(Pi.T, D), matrix_product(D, Pi), Pi.T]
            actual, expected = ambient.copy(), diagonal.copy()
            for group, coefficient in enumerate(coefficients):
                actual = physical_operator(actual, coefficient, direct_images[group % 2], group >= 2)
                expected = physical_operator(expected, coefficient, groups[group][0], group >= 2)
            actual = matrix_product(matrix_product(left, actual), right)
            assert np.array_equal(actual, expected), (t, place, np.count_nonzero(actual ^ expected))
            checked += actual.size
    return dict(literal_transformed_matrix_entries=checked, physical_code_rank=1,
                physical_code_is_nonhomogeneous=True, opposite_cut_hot_weights_checked=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='passed', factoring_exponent_improved=False,
                  concentration=concentration_controls(), literal_control=literal_control())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
