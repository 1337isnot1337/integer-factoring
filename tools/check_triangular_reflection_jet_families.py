#!/usr/bin/env python3
"""Exact controls for a uniform large linear family of reflection frames."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, coupled_admissible_poles
from reflection_low_jet_frames import reflection_jets, matrix_coefficient, pole_conditions
from square_zero_coupled_reflection_frames import matrix_terms, rank, product


def triangular_entry(interval, row, column):
    """Entry of (1+w) tensor_r [[1+w^(2^(r+1)),0],[w^(2^(r+1)),1]]."""
    value = 1
    for bit in range(interval.bit_length() - 1):
        a, b = (row >> bit) & 1, (column >> bit) & 1
        power = 1 << (bit + 1)
        if a == 0 and b == 1:
            return 0
        if a == b == 0:
            value ^= value << power
        elif a == 1 and b == 0:
            value <<= power
    return value ^ (value << 1)


def upper_block(q, dense=False):
    half = q // 2
    generator = [[0] * q for _ in range(q)]
    for row in range(half):
        for column in range(half):
            if row.bit_count() % 2 == column.bit_count() % 2:
                continue
            if dense:
                generator[row][half + column] = 1 + (7 * row + 3 * column) % 3
            elif column == row ^ 1:
                generator[row][half + column] = 1
    return generator


def negative_matrix_coefficients(jets, generators):
    size, maximum = jets['alphabet'], len(generators)
    coefficients = [matrix_coefficient(jets['B'], r) for r in range(2 * maximum)]
    result = {}
    for exponent in range(-2 * maximum, 0):
        output = [[0] * size for _ in range(size)]

        def add(value):
            for i in range(size):
                for j in range(size):
                    output[i][j] ^= value[i][j]

        for j, generator in enumerate(generators, 1):
            if exponent + j >= 0:
                coefficient = coefficients[exponent + j]
                add(product(generator, coefficient))
                add(product(coefficient, generator))
        for i, first in enumerate(generators, 1):
            for j, second in enumerate(generators, 1):
                if exponent + i + j >= 0:
                    coefficient = coefficients[exponent + i + j]
                    add(product(product(first, coefficient), second))
        result[exponent] = output
    return result


def independent_laurent_controls():
    from check_dense_coupled_reflection_frames import coefficients, apply_matrix

    bell, _ = bell_normalize(tensor(4))
    pairs = [
        ([[1, 2], [3, 1]], [[2, 1], [1, 3]]),
        ([[0, 3], [2, 1]], [[1, 1], [2, 0]]),
    ]
    entries = 0
    for blocks in pairs:
        generators = []
        for block in blocks:
            generator = [[0] * 4 for _ in range(4)]
            for i in range(2):
                for j in range(2):
                    generator[i][j + 2] = block[i][j]
            generators.append(generator)
        original = [coefficients(value) for value in bell['values']]
        values = original
        for repetition in range(2):
            for leg in range(4):
                for power, generator in enumerate(generators, 1):
                    matrix = generator if leg in (0, 2) else list(map(list, zip(*generator)))
                    values = apply_matrix(values, matrix, leg, 4, -power, identity=True)
            if repetition == 0:
                assert all(all(power >= 0 for power in value) for value in values)
                assert max(max(value, default=-1) for value in values) == 24
                for index, value in enumerate(values):
                    endpoint = 0
                    for coefficient in value.values():
                        endpoint ^= coefficient
                    x, y, u, v = (index // 4 ** leg % 4 for leg in range(4))
                    assert endpoint == int(x == v and y == u)
            else:
                assert values == original
        entries += len(values)
    jets = reflection_jets(8, cutoff=8)
    q, half = 16, 8
    first = [[0] * q for _ in range(q)]
    second = [[0] * q for _ in range(q)]
    for i in range(half):
        for j in range(half):
            first[i][half + j] = (i + 2 * j) % 4
            second[i][half + j] = (3 * i + j + 1) % 4
    generators = [first, second, upper_block(q, True), upper_block(q, False)]
    negative = negative_matrix_coefficients(jets, generators)
    assert all(not any(value for row in matrix for value in row) for matrix in negative.values())
    invalid = [[row[:] for row in generator] for generator in generators]
    invalid[2][0][8] ^= 1  # This changes total parity.
    negative_invalid = negative_matrix_coefficients(jets, invalid)
    assert any(value for row in negative_invalid[-1] for value in row)
    return {'four_letter_independent_coefficient_cases': len(pairs),
            'four_letter_endpoint_degree_involution_entries': entries,
            'sixteen_letter_mixed_negative_coefficient_entries': 2 * 8 * q * q,
            'nonbinary_upper_blocks': pairs,
            'odd_parity_highest_coefficient_rejected': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    one_particle_entries = 0
    for interval in (4, 8, 16, 32, 64):
        for row in range(interval):
            for column in range(interval):
                value = 0
                for middle in range(interval):
                    if row & middle == 0:
                        value ^= 1 << (2 * (interval - 1 - (middle | column)))
                value ^= value << 1
                assert value == triangular_entry(interval, row, column)
                one_particle_entries += 1
        middle = interval // 2 - 1
        for row in range(interval):
            value = triangular_entry(interval, row, middle)
            if row == middle:
                assert value == (1 ^ 2 ^ (1 << interval) ^ (1 << (interval + 1)))
            elif row == interval - 1:
                assert value == ((1 << interval) ^ (1 << (interval + 1)))
            else:
                assert value == 0
    cases = []
    for interval in (4, 8, 16):
        jets = reflection_jets(interval, cutoff=interval)
        q, half = jets['alphabet'], jets['alphabet'] // 2
        parity = [[int(i == j and i.bit_count() % 2) for j in range(q)] for i in range(q)]
        assert matrix_coefficient(jets['B'], 2) == parity
        assert jets['Q0'][0][q - 1] == 1
        for i in range(q):
            for j in range(q):
                assert not jets['Q0'][i][j] or i.bit_count() % 2 == j.bit_count() % 2
        assert all(jets['B'][i][j] == 0 for i in range(half, q) for j in range(half))
        dimension = sum(i.bit_count() % 2 == j.bit_count() % 2
                        for i in range(half) for j in range(half, q))
        assert dimension == q * q // 8
        record = {'interval': interval, 'alphabet': q, 'parity_even_parameter_dimension': dimension,
                  'B2_entries_checked': q * q,
                  'invariant_subspace_entries_checked_through_order_M_minus_1': q * q // 4,
                  'B_truncated_sha256': sha256(json.dumps(jets['B'], separators=(',', ':')).encode()).hexdigest()}
        if interval <= 8:
            bell, _ = bell_normalize(tensor(interval))
            controls = []
            for dense in (False, True):
                generator = upper_block(q, dense)
                maximum = min(4, interval // 2)
                predicted = pole_conditions(jets, generator, maximum)
                assert all(case['polynomial'] for case in predicted)
                terms = matrix_terms(bell, generator)
                literal = coupled_admissible_poles(terms, maximum)
                assert literal == list(range(1, maximum + 1))
                controls.append({'dense': dense, 'rank': rank(generator),
                                 'allowed_poles_checked': literal, 'generator': generator})
            record['full_tensor_controls'] = controls
        cases.append(record)
    independent = independent_laurent_controls()
    names = (Path(__file__).name, 'reflection_low_jet_frames.py',
             'reflection_dual_rail_gaussian.py', 'reflection_laurent_physical_frames.py',
             'square_zero_coupled_reflection_frames.py', 'check_dense_coupled_reflection_frames.py')
    result = {'status': 'passed', 'scope': 'uniform triangular proof with exact finite physical-jet and arbitrary-rank full-tensor controls; no pressure improvement',
              'one_particle_product_entries': one_particle_entries, 'cases': cases,
              'B2_entries_checked': sum(case['B2_entries_checked'] for case in cases),
              'invariant_subspace_entries_checked': sum(case['invariant_subspace_entries_checked_through_order_M_minus_1'] for case in cases),
              'independent_laurent_coefficients': independent,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
