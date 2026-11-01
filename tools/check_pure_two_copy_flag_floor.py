#!/usr/bin/env python3
"""Exact F4 image-flag obstruction and the uniform pure two-copy5/7 floor."""

from collections import Counter
from fractions import Fraction
from itertools import product
import json

from check_f4_recentered_walks import literal_matrix
from f4_recentered_lattices import add_equation, assignment, regular_frame, unpack
from search_f4_recentered_crossings import digest, five_letter_source
from source_slab_module_saturation import insert_binary


def binary_inverse(rows, size):
    work = [row | (1 << (size + i)) for i, row in enumerate(rows)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i] >> column & 1)
        work[column], work[pivot] = work[pivot], work[column]
        for i in range(size):
            if i != column and work[i] >> column & 1:
                work[i] ^= work[column]
    assert all((row & ((1 << size) - 1)) == 1 << i for i, row in enumerate(work))
    return [row >> size for row in work]


def span_basis(rows):
    pivots = {}
    for row in rows:
        insert_binary(pivots, row)
    return pivots


def in_span(value, pivots):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            return False
        value ^= pivots[pivot]
    return True


def annihilator(rows, width=25):
    pivots = span_basis(rows)
    result = []
    for variable in range(width):
        if variable in pivots:
            continue
        value = 1 << variable
        for pivot, row in sorted(pivots.items()):
            value |= ((row & value).bit_count() % 2) << pivot
        assert all((value & row).bit_count() % 2 == 0 for row in rows)
        result.append(value)
    return result


def right_factor_system(left, checks):
    pivots, selected = {}, []
    for index, check in enumerate(checks):
        coefficients = []
        for column in range(5):
            value = 0
            for row in range(5):
                if check >> (5 * row + column) & 1:
                    value ^= left[row]
            coefficients.append(value)
        low = sum((value & 1) << i for i, value in enumerate(coefficients))
        high = sum((value >> 1) << i for i, value in enumerate(coefficients))
        before = len(pivots)
        assert add_equation(pivots, low, high, 0)
        if len(pivots) != before:
            selected.append(index)
    return pivots, selected


def f4_multiply(a, b):
    low = ((a & 1) & (b & 1)) ^ ((a >> 1) & (b >> 1))
    high = ((a & 1) & (b >> 1)) ^ ((a >> 1) & (b & 1)) ^ ((a >> 1) & (b >> 1))
    return low | (high << 1)


def pure_outer(left, right):
    values = [f4_multiply(a, b) for a in left for b in right]
    return sum((value & 1) << i for i, value in enumerate(values)), sum((value >> 1) << i for i, value in enumerate(values))


def tensor_matrix_bits(first, second):
    result = 0
    for i in range(25):
        if not first >> i & 1:
            continue
        x1, y1 = divmod(i, 5)
        for j in range(25):
            if second >> j & 1:
                x2, y2 = divmod(j, 5)
                result ^= 1 << (25 * (5 * x1 + x2) + 5 * y1 + y2)
    return result


def profile_and_flags():
    values, heights, energy_x, energy_y = five_letter_source()
    profiles, flags, all_images, extremals = [], [], [], []
    for side in range(2):
        matrix, scalar = literal_matrix(values, 5, 18, heights, side)
        orders, left, right, _ = regular_frame(matrix, 64)
        profile = tuple(order - scalar for order in orders)
        profiles.append(profile)
        constants = []
        for frame in (left, right):
            for row in frame:
                assert all(unpack(value)[1] == 0 for value in row)
            constants.append([sum((value & 1) << j for j, value in enumerate(row)) for row in frame])
        inverse_left, inverse_right = [binary_inverse(rows, 25) for rows in constants]
        images = [sum(((row >> i) & 1) << j for j, row in enumerate(inverse_left)) for i in range(25)]
        coimages = inverse_right
        indices = [i for i, value in enumerate(profile) if value <= min(profile) + 9]
        assert indices == list(range(11))
        flags.extend((images[:11], coimages[:11]))
        all_images.append(images)
        head = [sum(((value >> (2 * min(orders))) & 3) << j for j, value in enumerate(row)) for row in matrix]
        nonzero = [row for row in head if row]
        assert all(row == nonzero[0] for row in nonzero)
        support = sum(int(bool(row)) << i for i, row in enumerate(head))
        assert support == nonzero[0] == 17246593
        assert len(span_basis([(support >> (5 * i)) & 31 for i in range(5)])) == 5
        extremals.append({'minimum_order': min(profile), 'rank_one_cut_vector': support,
                          'physical_reshape_rank': 5})
    common = span_basis(flags[0])
    assert len(common) == 11
    for flag in flags:
        assert len(span_basis(flag)) == 11 and all(in_span(value, common) for value in flag)
    relative = [[value - min(profile) for value in profile] for profile in profiles]
    assert relative[0] == relative[1]
    coefficients = [Counter(relative[0])[i] for i in range(21)]
    assert energy_x == energy_y == (0, 10, 6, 2, 8)
    return {'reference_digest': digest(values), 'profiles': profiles, 'relative_orders': relative[0],
            'relative_profile_polynomial': coefficients, 'common_flag_basis': flags[0],
            'full_image_basis': all_images[0], 'extremals': extremals}


def integer_pressure_certificate(coefficients):
    convolution = [sum(coefficients[i] * coefficients[k - i]
                       for i in range(21) if 0 <= k - i < 21) for k in range(41)]
    convolution[10] -= 1
    assert min(convolution) >= 0 and convolution[0] == 1 and convolution[1] == 2
    assert 80 * 4097 ** 19 < 81 * 4096 ** 19
    assert 14 * 16 > 81

    def numerator_a(m, depth):
        denominator = 1 << depth
        x, y = m ** 10, denominator ** 10
        numerator = sum(value * x ** i * y ** (20 - i) for i, value in enumerate(coefficients) if value)
        return numerator * numerator - m ** 100 * denominator ** 300

    def numerator_q(m, depth):
        denominator = 1 << depth
        return (denominator ** 70 + m ** 14 * denominator ** 56 + m ** 42 * denominator ** 28
                + m ** 56 * denominator ** 14 + m ** 70)

    pending, leaves = [(1, 1)], []
    while pending:
        m, depth = pending.pop()
        denominator = 1 << depth
        if numerator_a(m, depth) ** 7 > numerator_q(m + 1, depth) ** 20 * denominator ** 1400:
            leaves.append((m, depth))
        else:
            assert depth < 20
            pending.extend(((2 * m, depth + 1), (2 * m + 1, depth + 1)))
    leaves.sort(key=lambda item: Fraction(item[0], 1 << item[1]))
    endpoint = Fraction(1, 2)
    for m, depth in leaves:
        assert Fraction(m, 1 << depth) == endpoint
        endpoint = Fraction(m + 1, 1 << depth)
    assert endpoint == 1 and len(leaves) == 60 and max(depth for _, depth in leaves) == 10

    # An overly ambitious3/4 lower rate is not supported by this bound.
    t = Fraction(4, 5)
    p = sum(value * t ** (3 * i) for i, value in enumerate(coefficients))
    a = p * p - t ** 30
    q = 1 + t ** 4 + t ** 12 + t ** 16 + t ** 20
    assert a * a < q ** 6
    assert 2 * 5 ** 2 > 7 ** 2
    return {'rate_floor': [5, 7], 'strictly_above_one_over_sqrt_two': True,
            'nonnegative_partition_bound_coefficients': convolution,
            'small_interval_upper_endpoint': [1, 2],
            'small_interval_integer_comparison': [80, 4097, 81, 4096, 19],
            'dyadic_intervals': leaves, 'interval_count': len(leaves), 'maximum_depth': 10,
            'overclaimed_three_quarters_rejected': True}


def verify():
    data = profile_and_flags()
    basis = data['common_flag_basis']
    checks = annihilator(basis)
    assert len(checks) == 14
    projective = [(0,) * pivot + (1,) + tail for pivot in range(5)
                  for tail in product(range(4), repeat=4 - pivot)]
    selected_rows = []
    for left in projective:
        pivots, selected = right_factor_system(left, checks)
        assert len(pivots) == len(selected) == 5
        selected_rows.append(selected)
    assert len(projective) == 341

    left = (0, 1, 1, 1, 0)
    larger = data['full_image_basis'][:12]
    pivots, _ = right_factor_system(left, annihilator(larger))
    free = next(i for i in range(5) if i not in pivots)
    low, high = assignment(pivots, 1 << free, 0)
    right = tuple(((low >> i) & 1) | (((high >> i) & 1) << 1) for i in range(5))
    outer = pure_outer(left, right)
    assert any(right) and all(in_span(value, span_basis(larger)) for value in outer)
    assert not all(in_span(value, span_basis(basis)) for value in outer)

    tensor_basis = span_basis([tensor_matrix_bits(a, b) for a in basis for b in basis])
    assert len(tensor_basis) == 121
    relative = data['relative_orders']
    low_pairs = [(i, j) for i in range(25) for j in range(25) if relative[i] + relative[j] <= 9]
    for i, j in low_pairs:
        assert i < 11 and j < 11
        assert in_span(tensor_matrix_bits(data['full_image_basis'][i], data['full_image_basis'][j]), tensor_basis)
    pressure = integer_pressure_certificate(data['relative_profile_polynomial'])
    return {'status': 'exact fixed-natural-gauge F4 pure two-copy5/7 pressure floor certified',
            'reference_digest': data['reference_digest'],
            'complete_reference_profiles': data['profiles'],
            'relative_profile_polynomial': data['relative_profile_polynomial'],
            'four_equal_constant_image_and_coimage_flags': True,
            'physical_image_flag_basis': basis,
            'flag_annihilator_basis': checks,
            'projective_left_factors_exhausted': len(projective),
            'full_rank_constraint_minor_rows': selected_rows,
            'first_rank_one_vector': {'flag_dimension': 12, 'relative_order': 10,
                                      'left_factor': left, 'right_factor': right},
            'extremal_coefficient_controls': data['extremals'],
            'tensor_square_flag_dimension': 121,
            'two_copy_orders_at_relative_cutoff9': len(low_pairs),
            'pressure_certificate': pressure,
            'scope': 'One physically pure F4(z) correction to the tensor square of the retained seed; fixed natural physical heights; arbitrary Laurent leading orders are allowed.',
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
