#!/usr/bin/env python3
"""Replay the exact common-family cube construction and its rational bound."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, log, log2, prod
from pathlib import Path
import time

from cube_gaussian_crossings import endpoint, polynomial_product, rank_factorization
from cube_homogeneous_gaussian import (
    augmentation_crossblock, axis_matrix, cayley_orders, determinantal_orders,
    exact_monomial_smith, gaussian_cut_factors, gaussian_tensor, kernel,
    matrix_product, monomial_profile, multiply, serializable,
)


def binary_minor(matrix, rows, columns):
    if not rows:
        return 1
    result = 0
    for index, column in enumerate(columns):
        result ^= polynomial_product(matrix[rows[0]][column],
                                     binary_minor(matrix, rows[1:], columns[:index] + columns[index + 1:]))
    return result


def literal_local_orders(degree, first, second, triple):
    a, b, c, d = [1 << value for value in (degree, first, second, triple)]
    matrix = [[a, 0, 0, 0], [b, a, 0, 0], [c, 0, a, 0], [d, c, b, a]]
    previous, orders, count = 0, [], 0
    for size in range(1, 5):
        valuations = []
        for rows in combinations(range(4), size):
            for columns in combinations(range(4), size):
                value = binary_minor(matrix, rows, columns)
                if value:
                    valuations.append((value & -value).bit_length() - 1)
                count += 1
        current = min(valuations)
        orders.append(current - previous)
        previous = current
    return orders, count


def finite_typed_certificate(tensor):
    # Five sites and the two weight-one letters 00001 and 00010. The three
    # vacuum sites contribute one, leaving the two complementary factors.
    typed = tuple(multiply(tensor[mask], tensor[mask ^ 255]) for mask in range(256))
    assert all(len(value) % 2 == endpoint(mask) for mask, value in enumerate(typed))
    assert all(sum(exponent) == 24 for value in typed for exponent in value)
    matrices = [axis_matrix(typed, axis) for axis in range(3)]
    exponents = sorted(set().union(*typed))
    pieces = []
    for exponent in exponents:
        ranks, factors = [], []
        for matrix in matrices:
            rows = [sum(int(exponent in value) << j for j, value in enumerate(row)) for row in matrix]
            left, right = rank_factorization(rows, 16)
            ranks.append(len(right))
            factors.append({'left_rows': left, 'right_rows': right})
        assert all(ranks)
        pieces.append({'exponent': list(exponent), 'ranks': ranks, 'factors': factors})
    return {'sites': 5, 'weight': 1, 'selected_alphabet': [1, 2],
            'homogeneous_degree': 24, 'coefficient_count': len(pieces),
            'endpoint_entries': 256, 'factorization_entries': len(pieces) * 3 * 256,
            'maximum_rank_product': max(prod(piece['ranks']) for piece in pieces),
            'pieces': pieces}


def pressure_floor_certificate():
    """All-scale rational certificate for the scoped common-profile floor 3/5."""
    def product(first, second):
        result = [0] * (len(first) + len(second) - 1)
        for i, a in enumerate(first):
            for j, b in enumerate(second):
                result[i + j] += a * b
        return result

    partition = [1]
    for degree in (3, 6, 12, 15):
        factor = [0] * (degree + 1)
        factor[0] = factor[degree] = 1
        partition = product(partition, factor)
    difference = [1]
    for _ in range(5):
        difference = product(difference, partition)
    for index in range(13):
        difference[5 * index] -= comb(12, index)
    vanishing_order = 0
    while difference[0] == 0:
        difference = difference[1:]
        vanishing_order += 1
    assert vanishing_order == 3
    degree = len(difference) - 1
    transformed = [sum(coefficient * comb(degree - i, k - i)
                       for i, coefficient in enumerate(difference[:k + 1]))
                   for k in range(degree + 1)]
    assert all(coefficient > 0 for coefficient in transformed)
    return {
        'scoped_rho_floor': '3/5',
        'difference': '[(1+s^3)(1+s^6)(1+s^12)(1+s^15)]^5 - (1+s^5)^12',
        'vanishing_order': vanishing_order,
        'quotient_degree': degree,
        'positive_transform': '(1+u)^177 * (difference/s^3)(u/(1+u))',
        'positive_transform_coefficients': transformed,
        'minimum_transform_coefficient': min(transformed),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='evidence/cube-homogeneous-gaussian.json')
    parser.add_argument('--local-degree', type=int, default=8)
    parser.add_argument('--global-degree', type=int, default=12)
    args = parser.parse_args()
    started = time.monotonic()
    source = kernel()
    tensor = gaussian_tensor(source)
    assert all(len(value) % 2 == endpoint(mask) for mask, value in enumerate(tensor))
    assert all(sum(exponent) == 3 * mask.bit_count()
               for mask, value in enumerate(tensor) for exponent in value)
    axis_records, complete_gaussian_entries, minor_count = [], 0, 0
    for axis in range(3):
        block = augmentation_crossblock(source, axis)
        frame = exact_monomial_smith(block, axis)
        literal, checked = determinantal_orders(block, axis)
        assert literal == frame['orders'] == [2, 4, 8, 10]
        minor_count += checked
        left, cross_compound, right, _ = gaussian_cut_factors(source, axis)
        assert matrix_product(matrix_product(left, cross_compound), right) == axis_matrix(tensor, axis)
        complete_gaussian_entries += 256
        fock_orders = sorted(sum(frame['orders'][i] for i in range(4) if mask >> i & 1)
                             for mask in range(16))
        axis_records.append({'axis': axis, 'crossblock': block, 'smith_frame': frame,
                             'fock_orders': fock_orders})

    local_cases, local_minors = 0, 0
    for degree in range(1, args.local_degree + 1):
        for first in range(degree + 1):
            for second in range(degree + 1):
                for triple in range(degree + 1):
                    actual, checked = literal_local_orders(degree, first, second, triple)
                    assert actual == cayley_orders(degree, first, second, triple)
                    local_cases += 1
                    local_minors += checked

    global_cases, equality_profiles = 0, 0
    for degree in range(1, args.global_degree + 1):
        for xy_x in range(degree + 1):
            for xz_x in range(degree + 1):
                for yz_y in range(degree + 1):
                    for triple_x in range(degree + 1):
                        for triple_y in range(degree - triple_x + 1):
                            profiles = monomial_profile(degree, xy_x, xz_x, yz_y, triple_x, triple_y)
                            assert all(profile == sorted(profile) and sum(profile) == 4 * degree
                                       for profile in profiles)
                            assert sum(profile[0] for profile in profiles) <= degree
                            assert sum(profile[0] + profile[1] for profile in profiles) <= 3 * degree
                            if all(profile == [degree // 3, 2 * degree // 3, 4 * degree // 3, 5 * degree // 3]
                                   for profile in profiles) and degree % 3 == 0:
                                equality_profiles += 1
                            global_cases += 1

    scale, frequency = Fraction(3, 7), Fraction(1, 5)
    partition = prod(1 + scale ** exponent for exponent in (1, 2, 4, 5))
    left = partition ** 10 * 2 ** 40
    right = 5 ** 25 * scale ** 8
    assert left < right
    ratio = left / right
    entropy = log2(5) - 8 / 5
    numerical_rho = (log2(float(partition)) - 4 / 5 * log2(float(scale))) / (4 * entropy)
    assert numerical_rho < 5 / 8
    typed = finite_typed_certificate(tensor)
    pressure_floor = pressure_floor_certificate()
    source_paths = [Path(__file__), Path(__file__).with_name('cube_homogeneous_gaussian.py'),
                    Path(__file__).with_name('cube_gaussian_crossings.py')]
    evidence = {
        'schema': 'cube-homogeneous-gaussian-v1',
        'coefficient_field': 'F2',
        'kernel_homogeneous_degree': 6,
        'one_particle_endpoint_entries': 64,
        'gaussian_endpoint_entries': 256,
        'gaussian_monomial_terms': sum(map(len, tensor)),
        'complete_gaussian_factorization_entries': complete_gaussian_entries,
        'crossblock_minors': minor_count,
        'axes': axis_records,
        'local_degree_through': args.local_degree,
        'literal_local_profiles': local_cases,
        'literal_local_minors': local_minors,
        'global_degree_through': args.global_degree,
        'common_monomial_lists_checked': global_cases,
        'balanced_equality_profiles': equality_profiles,
        'typed_binary_example': typed,
        'common_scale_cayley_pressure_floor': pressure_floor,
        'pressure': {
            'physical_frequency': str(frequency), 'scale_after_deflation': str(scale),
            'partition_function': '(1+r)(1+r^2)(1+r^4)(1+r^5)',
            'exact_inequality': 'Z(r)^10 * 2^40 < 5^25 * r^8',
            'left_over_right_numerator': ratio.numerator,
            'left_over_right_denominator': ratio.denominator,
            'certified_rho_upper_bound': '5/8',
            'numerical_rho_upper_bound': numerical_rho,
            'cube_transfer_alpha_from_rational_bound': log(8 / 5) / log(3),
        },
        'new_uniform_family': True,
        'new_factoring_score': False,
        'scope': 'Actual common eight-leg family; its exponent is below the existing one-half result.',
        'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in source_paths},
        'elapsed_seconds': time.monotonic() - started,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(serializable(evidence), indent=2) + '\n')
    print(json.dumps({'axes': [record['smith_frame']['orders'] for record in axis_records],
                      'literal_local_profiles': local_cases, 'common_monomial_lists': global_cases,
                      'typed_coefficients': typed['coefficient_count'],
                      'rho_bound': numerical_rho, 'seconds': evidence['elapsed_seconds']}))


if __name__ == '__main__':
    main()
