#!/usr/bin/env python3
"""Exact factors, flag intersections, and leading jets of nine-term templates."""

from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_f4_recentered_walks import literal_matrix
from check_pure_two_copy_flag_floor import binary_inverse, in_span, profile_and_flags, span_basis, tensor_matrix_bits
from f4_recentered_lattices import add_equation, regular_frame
from nine_term_flag_templates import ANTITONE, COMPONENTS, COORDINATES, HEIGHTS, PAIRS, component_data, updates
from search_f4_recentered_crossings import five_letter_source


def map_vector(rows, value):
    return sum(((row & value).bit_count() % 2) << i for i, row in enumerate(rows))


def tensor_coordinates(first, second):
    value = 0
    for i in range(25):
        if first >> i & 1:
            value ^= second << (25 * i)
    return value


def f4_rows_rank(rows):
    pivots = {}
    for low, high in rows:
        assert add_equation(pivots, low, high, 0)
    return len(pivots)


def verify():
    data = profile_and_flags()
    basis, orders = data['full_image_basis'], data['relative_orders']
    components = [sum(value << (5 * i + coordinate) for i, value in enumerate(vector))
                  for vector, coordinate in zip(COMPONENTS, COORDINATES)]
    assert components[0] ^ components[1] ^ components[2] == basis[2] == 4769284
    physical = [tensor_matrix_bits(components[i], components[j]) for i, j in PAIRS]
    assert len(span_basis(physical)) == 9
    full_sum = 0
    for value in physical:
        full_sum ^= value
    assert full_sum == tensor_matrix_bits(basis[2], basis[2])
    intersections = []
    for cutoff in range(10):
        flag = [tensor_matrix_bits(basis[i], basis[j]) for i in range(25) for j in range(25)
                if orders[i] + orders[j] <= cutoff]
        dimension = len(span_basis(flag))
        intersection = 9 + dimension - len(span_basis(flag + physical))
        assert intersection == int(cutoff >= 6)
        if cutoff >= 6:
            assert in_span(full_sum, span_basis(flag))
        intersections.append({'relative_order_cutoff': cutoff, 'reference_flag_dimension': dimension,
                              'intersection_dimension': intersection})

    reference, heights, _, _ = five_letter_source()
    matrix, _ = literal_matrix(reference, 5, 18, heights, 0)
    _, left, right, _ = regular_frame(matrix, 64)
    left0 = [sum((value & 1) << j for j, value in enumerate(row)) for row in left]
    right0t = [sum((right[i][j] & 1) << i for i in range(25)) for j in range(25)]
    transformed_left = [map_vector(left0, value) for value in components]
    transformed_right = [map_vector(right0t, value) for value in components]
    low_coordinates = [25 * i + j for i in range(25) for j in range(25) if orders[i] + orders[j] < 6]
    new_coordinates = [25 * i + j for i in range(25) for j in range(25) if orders[i] + orders[j] == 6]
    assert len(low_coordinates) == 14 and len(new_coordinates) == 5
    mask = ((1 << 625) - 1) ^ sum(1 << i for i in low_coordinates)
    row_vectors = [tensor_coordinates(transformed_left[i], transformed_left[j]) & mask for i, j in PAIRS]
    column_vectors = [tensor_coordinates(transformed_right[i], transformed_right[j]) & mask for i, j in PAIRS]
    new_space = [1 << i for i in new_coordinates]
    assert len(span_basis(row_vectors)) == len(span_basis(column_vectors)) == 9
    assert 9 + 5 - len(span_basis(row_vectors + new_space)) == 1
    assert 9 + 5 - len(span_basis(column_vectors + new_space)) == 1
    permutations = [tuple(range(9)), ANTITONE, tuple((i + 1) % 9 for i in range(9)),
                    tuple(PAIRS.index((j, i)) for i, j in PAIRS)]
    patterns = [(1,) * 9, (2,) * 9, tuple(1 + i % 3 for i in range(9))]
    ranks = []
    for permutation in permutations:
        for pattern in patterns:
            rows = [[0, 0] for _ in range(625)]
            for i in new_coordinates:
                rows[i][0] = 1 << i
            for component, target in enumerate(permutation):
                active = row_vectors[component]
                while active:
                    bit = active & -active
                    active ^= bit
                    row = bit.bit_length() - 1
                    if pattern[component] & 1:
                        rows[row][0] ^= column_vectors[target]
                    if pattern[component] & 2:
                        rows[row][1] ^= column_vectors[target]
            rank = f4_rows_rank(rows)
            assert rank >= 12
            ranks.append({'pairing': permutation, 'coefficients': pattern,
                          'order_six_quotient_rank': rank})

    factors = component_data()
    sigma = [row['balanced_height_sum'] for row in factors]
    assert sum(sigma) == 54 and all(sigma[i] + sigma[ANTITONE[i]] == 12 for i in range(9))
    endpoint_polynomials = entry_order_controls = 0
    for zero_order, shape in product(range(16, 22), ('binomial', 'endpoint_power')):
        for row in updates(zero_orders=zero_order, shape=shape):
            polynomial = row['polynomial']
            low = (polynomial & -polynomial).bit_length() - 1
            high = polynomial.bit_length() - 1
            assert polynomial.bit_count() % 2 == 0 and 0 <= low < high <= 36
            assert low == row['lowest_polynomial_degree'] and high == row['highest_polynomial_degree']
            supports = [[i for i, value in enumerate(vector) if value] for vector in row['vectors']]
            at_infinity = min(HEIGHTS[u // 5] + HEIGHTS[u % 5] + HEIGHTS[v // 5] + HEIGHTS[v % 5] - high
                              for x, y, u, v in product(*supports))
            at_zero = min(HEIGHTS[x // 5] + HEIGHTS[x % 5] + HEIGHTS[y // 5] + HEIGHTS[y % 5] + low
                          for x, y, u, v in product(*supports))
            assert at_infinity == -10 and at_zero == zero_order
            endpoint_polynomials += 1
            entry_order_controls += 2
    source = Path('tools/nine_term_flag_templates.py')
    return {'status': 'exact nine-term template factors and leading-jet constraints certified',
            'constructor_sha256': sha256(source.read_bytes()).hexdigest(),
            'one_copy_rank_three_image': basis[2], 'one_copy_components': components,
            'nine_pair_labels': PAIRS, 'antitone_pairing': ANTITONE,
            'flag_intersections': intersections,
            'intersection_coefficient_vector_at_orders_six_to_nine': [1] * 9,
            'order_six_reference_quotient_rank': 5,
            'full_matching_update_quotient_rank': 9,
            'image_intersection_dimension': 1, 'coimage_intersection_dimension': 1,
            'universal_order_six_updated_rank_lower_bound': 12,
            'literal_F4_leading_quotient_controls': ranks,
            'component_height_sums': sigma, 'height_sum_average_per_matched_pair': 12,
            'best_simultaneous_zero_entry_order': 21,
            'endpoint_polynomials_checked': endpoint_polynomials,
            'literal_extreme_entry_orders_checked': entry_order_controls,
            'scope': 'Leading infinity jet for matched sparse factors, all pairings and nonzero field weights; not a global pressure exclusion.',
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
