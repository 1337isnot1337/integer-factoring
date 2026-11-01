#!/usr/bin/env python3
"""Exhaust arbitrary binary polynomial nilpotent frames of the q=4 seed."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_triangular_auxiliary_row_fiber import HEIGHT
from joint_height_conjugations import binary_product, restored_endpoints, source_matrix, xor_sum


ROOT = Path(__file__).resolve().parents[1]


def small_product(left, right):
    return tuple(xor_sum(right[k] for k in range(len(right)) if row >> k & 1) for row in left)


def nilpotent_height_matrices():
    small_two = [rows for rows in product(range(4), repeat=2)
                 if not any(small_product(rows, rows))]
    small_four = [rows for rows in product(range(16), repeat=4)
                  if not any(small_product(rows, rows))]
    assert len(small_two) == 4 and len(small_four) == 316
    middle = (3, 12, 6, 9)
    pairs = ((2, 8), (1, 4), (11, 14), (7, 13))
    for block in small_four:
        for smaller in product(small_two, repeat=4):
            rows = [0]*16
            for i, ambient in enumerate(middle):
                rows[ambient] = xor_sum(1 << middle[j] for j in range(4) if block[i] >> j & 1)
            for indices, local in zip(pairs, smaller):
                for i, ambient in enumerate(indices):
                    rows[ambient] = xor_sum(1 << indices[j] for j in range(2) if local[i] >> j & 1)
            yield tuple(rows)


def binary_rank(columns):
    basis = {}
    for value in columns:
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def packed_matrix(rows):
    return sum(row << (16*i) for i, row in enumerate(rows))


def frame_columns(commutators, products, degree):
    # Negative powers -1,...,-2D each receive one disjoint 256-bit block.
    # The final bit is f(1)=0. Squaring is F_2-linear on the coefficients.
    columns = []
    for j in range(degree+1):
        value = 1 << (512*degree)
        for ell in range(min(degree+1, j)):
            value ^= commutators[ell] << (256*(j-ell-1))
        for ell in range(min(degree+1, 2*j)):
            value ^= products[ell] << (256*(2*j-ell-1))
        columns.append(value)
    return columns


def direct_negative_vector(polynomial, commutators, products, degree):
    support = [j for j in range(degree+1) if polynomial >> j & 1]
    squared = set()
    for i, j in product(support, repeat=2):
        exponent = i+j
        if exponent in squared:
            squared.remove(exponent)
        else:
            squared.add(exponent)
    correction = {}
    for j in support:
        for ell, value in enumerate(commutators):
            correction[ell-j] = correction.get(ell-j, 0) ^ value
    for j in squared:
        for ell, value in enumerate(products):
            correction[ell-j] = correction.get(ell-j, 0) ^ value
    vector = (len(support) % 2) << (512*degree)
    for exponent, value in correction.items():
        if exponent < 0:
            vector ^= value << (256*(-exponent-1))
    return vector


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = source_matrix()
    zero, infinity = restored_endpoints(source)
    degree = max(value.bit_length()-1 for row in infinity for value in row)
    assert degree == 8
    coefficients = []
    for ell in range(degree+1):
        coefficients.append(tuple(sum(
            ((source[i][j] >> (HEIGHT[j//4]+HEIGHT[j%4]+4-ell)) & 1) << j
            for j in range(16) if HEIGHT[j//4]+HEIGHT[j%4]+4 >= ell)
            for i in range(16)))
    restoration_checks = 0
    for x, y, u, v in product(range(4), repeat=4):
        restored = sum(((coefficients[ell][4*x+y] >> (4*u+v)) & 1) << ell
                       for ell in range(degree+1))
        assert restored == infinity[4*x+u][4*y+v]
        restoration_checks += 1
    digest = sha256()
    total = eligible = direct_controls = 0
    unchanged = []
    ranks = Counter()
    for matrix in nilpotent_height_matrices():
        total += 1
        assert not any(binary_product(matrix, matrix))
        for i, row in enumerate(matrix):
            assert all(HEIGHT[i//4]+HEIGHT[i%4] == HEIGHT[j//4]+HEIGHT[j%4]
                       for j in range(16) if row >> j & 1)
        digest.update(b''.join(row.to_bytes(2, 'little') for row in matrix))
        if any(binary_product(binary_product(matrix, coefficients[0]), matrix)):
            continue
        eligible += 1
        commutators, products = [], []
        for coefficient in coefficients:
            left = binary_product(matrix, coefficient)
            right = binary_product(coefficient, matrix)
            commutators.append(packed_matrix(tuple(a ^ b for a,b in zip(left,right))))
            products.append(packed_matrix(binary_product(left, matrix)))
        if not any(commutators) and not any(products):
            unchanged.append(matrix)
            continue
        columns = frame_columns(commutators, products, degree)
        rank = binary_rank(columns)
        ranks[rank] += 1
        assert rank == degree+1
        if direct_controls < 256:
            for polynomial in (0, 1, 3, 5, 17, 255, 256, 257, 511):
                encoded = xor_sum(columns[j] for j in range(degree+1) if polynomial >> j & 1)
                assert encoded == direct_negative_vector(polynomial, commutators, products, degree)
                direct_controls += 1
    assert total == 80896 and eligible == 6400 and len(unchanged) == 32
    assert dict(ranks) == {9: 6368}
    result = {'status': 'pass',
              'scope': 'All F_2 polynomial frames I+f(z)N with f(1)=0, constant height-preserving N^2=0, and the original infinity charges',
              'height_preserving_square_zero_matrices': total,
              'necessary_leading_pole_passes': eligible,
              'changing_matrices_with_full_rank_frame_constraints': sum(ranks.values()),
              'constraint_rank_histogram': dict(ranks),
              'unchanged_matrices': unchanged,
              'maximum_nontrivial_polynomial_degree_needed': degree,
              'infinity_restoration_entries': restoration_checks,
              'direct_laurent_convolution_controls': direct_controls,
              'enumeration_sha256': digest.hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'dependencies_sha256': {path: sha256((ROOT/path).read_bytes()).hexdigest() for path in (
                  'tools/joint_height_conjugations.py', 'tools/check_triangular_auxiliary_row_fiber.py')},
              'nontrivial_regular_frames': 0, 'factoring_exponent_improved': False}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'unchanged_matrices'}, indent=2))


if __name__ == '__main__':
    main()
