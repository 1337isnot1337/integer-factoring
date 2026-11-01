#!/usr/bin/env python3
"""Independent coefficient, quotient-rank, and retained-code block controls."""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

from check_schur_gaussian_reduction import polynomial_multiply
from hot_code_charge_profiles import compound_block, natural_profiles, one_particle
from native_matrix_smith_f4 import MatrixOracle


def exact_minor(matrix, rows, columns):
    """Independent characteristic-two Laplace recurrence on complete polynomials."""
    rows, columns = tuple(rows), tuple(columns)

    @lru_cache(None)
    def visit(index, remaining):
        if index == len(rows):
            return 1
        result = 0
        for position, column in enumerate(remaining):
            result ^= polynomial_multiply(matrix[rows[index]][column],
                                          visit(index + 1, remaining[:position] + remaining[position + 1:]))
        return result

    return visit(0, columns)


def binary_rank(rows):
    basis = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def quotient_rank(matrix, length):
    n = len(matrix)
    rows = []
    for coefficient in range(length):
        for i in range(n):
            value = 0
            for j in range(n):
                for k in range(coefficient + 1):
                    if matrix[i][j] >> k & 1:
                        value ^= 1 << ((coefficient - k) * n + j)
            rows.append(value)
    return binary_rank(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    counts = {'independent_complete_determinants': 0, 'block_endpoint_entries': 0,
              'independent_quotient_ranks': 0, 'retained_profile_channels': 0,
              'one_particle_channels': 0, 'invalid_inputs_rejected': 0}
    small_controls = []
    for first, second in (((0, 1), (0, 1, 2)), ((0, 1, 2), (0, 1)),
                          ((0, 1, 2, 3), (0, 1, 4)), ((0, 4), (0, 1, 2, 3))):
        result = natural_profiles(first, second, keep_blocks=True)
        modes = result['modes']
        source = one_particle(modes)
        for block in result['blocks']:
            charge = block['charge']
            pairs = [(x, y) for x in first for y in second if x.bit_count() + y.bit_count() == charge]
            row_masks = [x | (y << modes) for x, y in pairs]
            column_masks = [y | (x << modes) for x, y in pairs]
            matrix = block['matrix']
            for i, row_mask in enumerate(row_masks):
                rr = [a for a in range(2 * modes) if row_mask >> a & 1]
                for j, column_mask in enumerate(column_masks):
                    cc = [a for a in range(2 * modes) if column_mask >> a & 1]
                    assert matrix[i][j] == exact_minor(source, rr, cc)
                    assert matrix[i][j].bit_count() % 2 == int(i == j)
                    counts['independent_complete_determinants'] += 1
                    counts['block_endpoint_entries'] += 1
            orders = block['orders_t']
            lengths = sorted({1, max(orders) + 1} | {value for value in orders if value > 0})
            for length in lengths:
                assert quotient_rank(matrix, length) == sum(max(length - value, 0) for value in orders)
                counts['independent_quotient_ranks'] += 1
        small_controls.append({'x_letters': first, 'y_letters': second,
                               'orders_w': result['orders_w'],
                               'block_sizes': [block['size'] for block in result['blocks']]})
    one_particle_controls = []
    for modes in range(1, 17):
        matrix = [[(value, 0) for value in row] for row in one_particle(modes)]
        orders = MatrixOracle(matrix, 2 * modes).profile()['orders']
        assert orders == list(range(1, 4 * modes, 2))
        counts['one_particle_channels'] += len(orders)
        one_particle_controls.append({'modes': modes, 'orders_t': list(map(int, orders))})
    archive = Path(__file__).resolve().parents[1] / 'evidence/implicit-reflection-code-search-replay.json'
    raw = archive.read_bytes()
    retained_controls = []
    for previous in json.loads(raw)['natural_profiles']:
        letters = [int(format(value, '08b')[::-1], 2) for value in previous['letters']]
        result = natural_profiles(letters)
        expected = previous['integer_profiles_infinity_then_zero']
        assert all(result['orders_w'] == [value - min(orders) for value in orders] for orders in expected)
        assert min(expected[0]) + min(expected[1]) == 0
        counts['retained_profile_channels'] += sum(map(len, expected))
        retained_controls.append({'size': len(letters), 'hot_letters': letters,
                                  'tensor_sha256': previous['tensor_sha256'],
                                  'orders_w': result['orders_w'],
                                  'block_sizes': [block['size'] for block in result['blocks']],
                                  'seconds': result['seconds']})
    bad_calls = [lambda: natural_profiles((1,)),
                 lambda: natural_profiles((0, 3)),
                 lambda: natural_profiles((0, 1, 1)),
                 lambda: compound_block(one_particle(2), (1,), (3,)),
                 lambda: compound_block(one_particle(2), (1,), (1,))]
    for function in bad_calls:
        try:
            function()
        except (ValueError, ArithmeticError):
            counts['invalid_inputs_rejected'] += 1
        else:
            raise AssertionError('An invalid or endpoint-incompatible input was accepted.')
    names = (Path(__file__).name, 'hot_code_charge_profiles.py', 'hot_code_charge_profiles.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'check_schur_gaussian_reduction.py')
    report = {'scope': 'Exact selected compound polynomials and complete natural profiles; finite one-particle pattern only, no all-size Smith formula or factoring improvement.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'retained_replay_sha256': sha256(raw).hexdigest(), 'counts': counts,
              'small_controls': small_controls, 'one_particle_controls': one_particle_controls,
              'retained_controls': retained_controls, 'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': counts}))


if __name__ == '__main__':
    main()
