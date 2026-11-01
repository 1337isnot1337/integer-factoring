#!/usr/bin/env python3
"""Independent scalar and full-tensor controls for opposite Bell CRT jets."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from bell_jet_crt_crossings import construct, scalar_data
from check_dense_coupled_reflection_frames import coefficients, matrix_product, multiply
from f16_identity_operator_profiles import inverse_matrix
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from screen_affine_finite_places import bitplane_rank
from search_entrywise_reflection_jets import basis_catalogue, bell_values


def field_power(value, exponent):
    result = 1
    for _ in range(exponent):
        result = multiply(result, value)
    return result


def hasse(polynomial, point, order):
    result = 0
    for exponent, coefficient in polynomial.items():
        if exponent >= order and comb(exponent, order) % 2:
            result ^= multiply(coefficient, field_power(point, exponent - order))
    return result


def add_polynomial(target, source):
    for exponent, coefficient in source.items():
        value = target.get(exponent, 0) ^ coefficient
        if value:
            target[exponent] = value
        elif exponent in target:
            del target[exponent]


def independent_check(seed, row_matrix, column_matrix):
    entries = list(map(coefficients, seed['values']))
    inverse_row, inverse_column = inverse_matrix(row_matrix), inverse_matrix(column_matrix)
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    assert matrix_product(row_matrix, inverse_row) == identity
    assert matrix_product(column_matrix, inverse_column) == identity
    traces = [[{} for _ in range(16)] for _ in range(2)]
    endpoints = derivatives = jets = 0
    for index, polynomial in enumerate(entries):
        x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
        assert hasse(polynomial, 1, 0) == int(x == v and y == u)
        assert hasse(polynomial, 1, 1) == 0
        endpoints += 1
        derivatives += 1
        expected_row = multiply(3, multiply(row_matrix[x][y], inverse_row[u][v]))
        expected_column = multiply(2, multiply(column_matrix[x][u], inverse_column[y][v]))
        for point, expected in ((2, expected_row), (3, expected_column)):
            for order in range(7):
                assert hasse(polynomial, point, order) == (expected if order == 0 else 0)
                jets += 1
        if x == v:
            add_polynomial(traces[0][4 * y + u], polynomial)
        if y == u:
            add_polynomial(traces[1][4 * x + v], polynomial)
    for trace in traces:
        for index, polynomial in enumerate(trace):
            expected = {0: 1, 16: 1} if index // 4 == index % 4 else {}
            assert polynomial == expected
    assert max(max(value, default=-1) for value in entries) == seed['degree'] == 16
    return endpoints, derivatives, jets, 32


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    scalar = scalar_data()
    # The constructor's bitwise Hasse rows are replayed with integer
    # binomial coefficients and an independent F4 multiplication.
    scalar_controls = 0
    expected_values = {'endpoint': {1: 1, 2: 0, 3: 0},
                       'row': {1: 0, 2: 3, 3: 0},
                       'column': {1: 0, 2: 0, 3: 2},
                       'kernel': {1: 0, 2: 0, 3: 0}}
    for name, value in scalar.items():
        polynomial = coefficients(value)
        for point, length in ((1, 2), (2, 7), (3, 7)):
            for order in range(length):
                expected = expected_values[name][point] if order == 0 else 0
                assert hasse(polynomial, point, order) == expected
                scalar_controls += 1
    assert max(coefficients(scalar['kernel'])) == 16
    matrices = basis_catalogue(tensor(4))
    matrices['alternating'] = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
    matrices['nonsymmetric'] = [[1, 2, 0, 0], [0, 1, 3, 0], [0, 0, 1, 1], [0, 0, 0, 1]]
    pairs = [('identity', 'identity'), ('source_Bell', 'identity'),
             ('identity', 'source_Bell'), ('source_Bell', 'source_Bell_inverse'),
             ('alternating', 'identity'), ('nonsymmetric', 'identity'),
             ('nonsymmetric', 'alternating'), ('evaluation', 'Pascal')]
    records, counts = [], [0, 0, 0, 0]
    for row_name, column_name in pairs:
        row, column = matrices[row_name], matrices[column_name]
        perturbed = bell_values(row, 0)
        # One entry with x!=v and y!=u has both paired traces zero.
        perturbed[0 + 4 * 1 + 16 * 2 + 64 * 1] ^= 2
        for name, carrier in (('row', 'row'), ('column', 'column'),
                              ('minimal', 'minimal'), ('perturbed', perturbed)):
            seed = construct(row, column, carrier)
            actual = independent_check(seed, row, column)
            counts = [a + b for a, b in zip(counts, actual)]
            records.append({'row_matrix': row_name, 'column_matrix': column_name,
                            'carrier': name, 'tensor_sha256': digest(seed['values']),
                            'degree': seed['degree']})
    trace_rows = [[0, 0] for _ in range(32)]
    for entry in range(256):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        if x == v:
            trace_rows[4 * y + u][0] ^= 1 << entry
        if y == u:
            trace_rows[16 + 4 * x + v][0] ^= 1 << entry
    trace_rank = bitplane_rank(trace_rows)
    assert trace_rank == 31 and 256 - trace_rank == 225
    names = (Path(__file__).name, 'bell_jet_crt_crossings.py',
             'check_dense_coupled_reflection_frames.py', 'entrywise_reflection_constraints.py',
             'search_entrywise_reflection_jets.py', 'screen_affine_finite_places.py')
    report = {'status': 'passed',
              'scope': 'Exact scalar CRT and full tensor endpoint, derivative, traces, and seventh-order Bell jets; no pressure or optimality claim.',
              'scalar_Hasse_controls': scalar_controls, 'full_tensor_cases': len(records),
              'independent_endpoint_entries': counts[0], 'independent_derivative_entries': counts[1],
              'independent_finite_Hasse_entries': counts[2], 'independent_polynomial_traces': counts[3],
              'constant_trace_rank': trace_rank, 'carrier_affine_dimension': 256 - trace_rank,
              'matrices': matrices, 'cases': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('matrices', 'cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
