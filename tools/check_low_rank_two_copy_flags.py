#!/usr/bin/env python3
"""Exact low-rank image flags and pressure floors for two-copy corrections.

The complete degree-five copy flag has 89,478,485 projective F4 vectors.
Its native enumeration normally takes a few minutes. Every field operation
is exact; pressure inequalities use integer endpoint comparisons only.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

from check_f4_recentered_walks import literal_matrix
from check_pure_two_copy_flag_floor import (
    binary_inverse, in_span, profile_and_flags, span_basis, tensor_matrix_bits,
)
from f4_recentered_lattices import add_equation, assignment, f4_rank, regular_frame
from search_f4_recentered_crossings import five_letter_source


def grassmann(n, dimension):
    for pivots in combinations(range(n), dimension):
        free = [(i, j) for i, pivot in enumerate(pivots)
                for j in range(pivot + 1, n) if j not in pivots]
        for coefficients in product(range(4), repeat=len(free)):
            rows = [[int(j == pivot) for j in range(n)] for pivot in pivots]
            for (i, j), coefficient in zip(free, coefficients):
                rows[i][j] = coefficient
            yield rows


def contract(forms, matrix):
    low = high = 0
    for a, row in enumerate(forms):
        low_row = high_row = 0
        for i, coefficient in enumerate(row):
            source = (matrix >> (5 * i)) & 31
            if coefficient & 1:
                low_row ^= source
            if coefficient & 2:
                high_row ^= source
        low |= low_row << (5 * a)
        high |= high_row << (5 * a)
    return low, high


def first_dependency(forms, basis, limit):
    pivots = {}
    for j, matrix in enumerate(basis[:limit]):
        before = len(pivots)
        low, high = contract(forms, matrix)
        assert add_equation(pivots, low, high, 0)
        if len(pivots) == before:
            return j + 1
    return None


def literal_rank(low, high, size):
    return f4_rank([[((low >> (size * i + j)) & 1) | (((high >> (size * i + j)) & 1) << 1)
                    for j in range(size)] for i in range(size)])


def combine(basis, coefficients):
    low = high = 0
    for matrix, coefficient in zip(basis, coefficients):
        if coefficient & 1:
            low ^= matrix
        if coefficient & 2:
            high ^= matrix
    return low, high


def one_copy_frontiers(data):
    basis, orders = data['full_image_basis'], data['relative_orders']
    result = []
    for rank, expected_count in ((1, 341), (2, 5797), (3, 5797), (4, 341)):
        earliest, best_forms, count = 26, None, 0
        for forms in grassmann(5, 5 - rank):
            count += 1
            at = first_dependency(forms, basis, earliest - 1)
            if at is not None:
                earliest, best_forms = at, forms
        assert count == expected_count
        # A separate pass checks the final excluded prefix directly.
        assert all(first_dependency(forms, basis, earliest - 1) is None
                   for forms in grassmann(5, 5 - rank))
        images = [contract(best_forms, matrix) for matrix in basis[:earliest]]
        constraints = {}
        for coordinate in range(5 * (5 - rank)):
            low = sum(((a >> coordinate) & 1) << j for j, (a, b) in enumerate(images))
            high = sum(((b >> coordinate) & 1) << j for j, (a, b) in enumerate(images))
            assert add_equation(constraints, low, high, 0)
        free = next(j for j in range(earliest) if j not in constraints)
        low, high = assignment(constraints, 1 << free, 0)
        coefficients = [((low >> j) & 1) | (((high >> j) & 1) << 1) for j in range(earliest)]
        witness = combine(basis[:earliest], coefficients)
        assert literal_rank(*witness, 5) == rank
        result.append({'rank_threshold': rank, 'grassmann_subspaces_exhausted': count,
                       'first_flag_dimension': earliest, 'first_relative_order': orders[earliest - 1],
                       'annihilating_forms_of_witness_column_space': best_forms,
                       'witness_coefficients': coefficients, 'witness_bit_planes': witness})
    assert [(row['first_flag_dimension'], row['first_relative_order']) for row in result] == [(12, 10), (6, 6), (3, 3), (2, 1)]
    return result


def all_constant_flags_equal(data):
    values, heights, _, _ = five_letter_source()
    reference = data['full_image_basis']
    checked = 0
    for side in range(2):
        matrix, scalar = literal_matrix(values, 5, 18, heights, side)
        orders, left, right, _ = regular_frame(matrix, 64)
        assert [order - min(orders) for order in orders] == data['relative_orders']
        for kind, frame in enumerate((left, right)):
            inverse = binary_inverse([sum((value & 1) << j for j, value in enumerate(row)) for row in frame], 25)
            actual = inverse if kind else [sum(((row >> i) & 1) << j for j, row in enumerate(inverse)) for i in range(25)]
            for size in range(1, 26):
                pivots = span_basis(reference[:size])
                assert len(pivots) == size
                assert all(in_span(value, pivots) for value in actual[:size])
                checked += 1
    return checked


def native_executable():
    source = Path(__file__).with_name('low_rank_copy_flag_enumeration.cpp')
    compiler = shutil.which('c++')
    if compiler is None:
        raise RuntimeError('The complete finite rank enumeration requires a C++17 compiler.')
    fingerprint = sha256(source.read_bytes()).hexdigest()
    executable = Path(tempfile.gettempdir()) / ('integer-factoring-copy-flag-ranks-' + fingerprint[:24])
    if not executable.exists():
        with tempfile.TemporaryDirectory(prefix='integer-factoring-copy-flag-build-') as directory:
            temporary = Path(directory) / 'enumerator'
            subprocess.run([compiler, '-O3', '-std=c++17', str(source), '-o', str(temporary)], check=True)
            temporary.replace(executable)
    return executable, fingerprint


def native_enumerate(executable, basis):
    assert 1 <= len(basis) <= 14 and all(0 <= value < 1 << 625 for value in basis)
    rows = [str(len(basis))]
    rows.extend(' '.join(str((value >> (25 * i)) & ((1 << 25) - 1)) for i in range(25)) for value in basis)
    payload = ('\n'.join(rows) + '\n').encode()
    result = json.loads(subprocess.check_output([str(executable)], input=payload))
    assert result['dimension'] == len(basis)
    assert result['projective_vectors_checked'] == (4 ** len(basis) - 1) // 3
    coefficients = result['witness_coefficients']
    assert len(coefficients) == len(basis) and any(coefficients)
    assert next(value for value in coefficients if value) == 1
    witness = combine(basis, coefficients)
    assert literal_rank(*witness, 25) == result['minimum_rank']
    result['input_sha256'] = sha256(payload).hexdigest()
    result['witness_bit_planes'] = witness
    return result


def native_controls(executable):
    identity = sum(1 << (26 * i) for i in range(25))
    # This binary pencil is nonsingular on every nonzero F2 combination,
    # but has rank-one members over F4 because x^2+x+1 splits there.
    small_identity = 1 | (1 << 26)
    quadratic = (1 << 1) | (1 << 25) | (1 << 26)
    rng = random.Random(90793)
    controls = [[identity], [1], [identity, identity], [small_identity, quadratic]]
    controls += [[rng.getrandbits(625) for _ in range(3)] for _ in range(3)]
    vectors = 0
    for basis in controls:
        expected = 26
        for first in range(len(basis)):
            for tail in product(range(4), repeat=len(basis) - first - 1):
                coefficients = [0] * first + [1] + list(tail)
                expected = min(expected, literal_rank(*combine(basis, coefficients), 25))
                vectors += 1
        actual = native_enumerate(executable, basis)
        assert actual['minimum_rank'] == expected
    assert all(literal_rank(*combine([small_identity, quadratic], coefficients), 25) == 2
               for coefficients in ((1, 0), (0, 1), (1, 1)))
    assert native_enumerate(executable, [small_identity, quadratic])['minimum_rank'] == 1
    return {'independent_literal_spaces': len(controls), 'independent_projective_rank_comparisons': vectors,
            'nonbinary_rank_drop_positive_control': True, 'dependent_basis_zero_control': True}


def pressure_certificate(p, rank, cutoff, numerator, denominator):
    m, n = numerator, denominator
    coefficients = [sum(p[i] * p[k - i] for i in range(21) if 0 <= k - i < 21) for k in range(41)]
    coefficients[cutoff] -= rank
    negative = [i for i, value in enumerate(coefficients) if value < 0]
    # The only signed case retained below is -3s^6, absorbed by 4s^5.
    if negative:
        assert negative == [6] and coefficients[5] + coefficients[6] >= 0
        assert 5 * coefficients[5] + 6 * coefficients[6] >= 0
    assert coefficients[0] == 1 and coefficients[1] == 2
    assert n > m
    bound = 16 * m + 1
    unit_denominator = 1 << (2 * n)
    assert 16 * m * (unit_denominator + 4) ** (4 * m - 1) < bound * unit_denominator ** (4 * m - 1)
    assert 2 * n * (1 << (2 * (n - m))) > bound

    @lru_cache(None)
    def power_a(a, depth):
        d = 1 << depth
        x, y = a ** (2 * m), d ** (2 * m)
        value = sum(coefficient * x ** i * y ** (20 - i) for i, coefficient in enumerate(p))
        value = value * value - rank * a ** (2 * m * cutoff) * d ** (2 * m * (40 - cutoff))
        assert value > 0
        return value ** n

    @lru_cache(None)
    def power_q(a, depth):
        d = 1 << depth
        return (d ** (10 * n) + a ** (2 * n) * d ** (8 * n) + a ** (6 * n) * d ** (4 * n)
                + a ** (8 * n) * d ** (2 * n) + a ** (10 * n)) ** (4 * m)

    pending, leaves = [(1, 1)], []
    while pending:
        a, depth = pending.pop()
        d = 1 << depth
        if power_a(a, depth) > power_q(a + 1, depth) * d ** (40 * m * n):
            leaves.append((a, depth))
        else:
            assert depth < 20
            pending.extend(((2 * a, depth + 1), (2 * a + 1, depth + 1)))
    leaves.sort(key=lambda pair: Fraction(pair[0], 1 << pair[1]))
    current = Fraction(1, 2)
    for a, depth in leaves:
        assert Fraction(a, 1 << depth) == current
        current = Fraction(a + 1, 1 << depth)
    assert current == 1
    assert 2 * m * m > n * n
    return {'pure_term_count': rank, 'relative_tail_cutoff': cutoff, 'rate_floor': [m, n],
            'lower_partition_coefficients': coefficients, 'monotonicity_negative_degrees': negative,
            'small_interval_derivative_bound': bound, 'dyadic_intervals': leaves,
            'interval_count': len(leaves), 'maximum_depth': max(depth for _, depth in leaves),
            'strictly_above_one_over_sqrt_two': True}


def verify():
    data = profile_and_flags()
    prefix_checks = all_constant_flags_equal(data)
    frontiers = one_copy_frontiers(data)
    executable, native_hash = native_executable()
    controls = native_controls(executable)
    basis, orders = data['full_image_basis'], data['relative_orders']
    copy_flags = []
    for cutoff, dimension, expected_rank in ((4, 10, 12), (5, 14, 12)):
        pairs = [(i, j) for i in range(25) for j in range(25) if orders[i] + orders[j] <= cutoff]
        matrices = [tensor_matrix_bits(basis[i], basis[j]) for i, j in pairs]
        assert len(matrices) == len(span_basis(matrices)) == dimension
        result = native_enumerate(executable, matrices)
        assert result['minimum_rank'] == expected_rank
        copy_flags.append({'relative_order_cutoff': cutoff, 'basis_pairs': pairs, **result})
    rank_nine = tensor_matrix_bits(basis[2], basis[2])
    assert orders[2] + orders[2] == 6 and literal_rank(rank_nine, 0, 25) == 9
    pressure = [pressure_certificate(data['relative_profile_polynomial'], *parameters)
                for parameters in ((2, 6, 5, 7), (3, 5, 5, 7), (5, 6, 5, 7), (8, 6, 17, 24))]
    return {'status': 'exact low-rank two-copy flags and fixed-gauge pressure floors certified',
            'reference_digest': data['reference_digest'],
            'native_enumerator_sha256': native_hash,
            'constant_image_and_coimage_prefix_equalities': prefix_checks,
            'one_copy_rank_frontiers': frontiers, 'native_independent_controls': controls,
            'two_copy_flags': copy_flags,
            'rank_nine_order_six_positive_control': {'basis_pair': [2, 2], 'rank': 9, 'bit_matrix': rank_nine},
            'pressure_certificates': pressure,
            'scope': 'Sums of at most eight physically pure F4(z) corrections to the retained tensor square at fixed natural heights; arbitrary Laurent leading orders.',
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
