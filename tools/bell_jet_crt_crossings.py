#!/usr/bin/env python3
"""Explicit degree16 crossings with opposite seventh-order Bell jets."""
from __future__ import annotations

from functools import lru_cache

from entrywise_reflection_constraints import hasse_row, local_affine, verify_target
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, degree, scale
from search_entrywise_reflection_jets import bell_jets, bell_values


def binary_product(first, second):
    result = 0
    while second:
        bit = second & -second
        second ^= bit
        result ^= first << (bit.bit_length() - 1)
    return result


@lru_cache(maxsize=1)
def scalar_data():
    coefficients = []
    for selected in (1, 2, 3):
        equations = []
        for point, multiplicity in ((1, 2), (2, 7), (3, 7)):
            for jet in range(multiplicity):
                equations.append((hasse_row(point, jet, 15), int(point == selected and jet == 0)))
        solved = local_affine(equations, 16)
        assert solved is not None and solved[1] == []
        coefficients.append(solved[0])
    h = 1
    for _ in range(7):
        h = binary_product(h, 7)
    k = binary_product(h, 5)
    e1, e2, e3 = coefficients
    assert e1 == (h << 1, 0)
    assert add(add(e1, e2), e3) == (1, 0)
    assert (e2[0] & 1, e2[1] & 1) == (0, 1)
    assert (e3[0] & 1, e3[1] & 1) == (1, 1)
    row, column = scale(e2, 3), scale(e3, 2)
    assert add(row, column) == ((1 | (1 << 16)) ^ k, 0)
    return {'endpoint': e1, 'row': row, 'column': column, 'kernel': (k, 0)}


def construct(row_matrix, column_matrix, carrier='column'):
    scalar = scalar_data()
    row, column = bell_values(row_matrix, 0), bell_values(column_matrix, 1)
    if carrier == 'row':
        correction = row
    elif carrier == 'column':
        correction = column
    elif carrier == 'minimal':
        correction = []
        for entry in range(256):
            x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
            correction.append(int(x == v == 0 and y == u) ^ int(x == v and y == u == 0))
    else:
        correction = list(carrier)
        if len(correction) != 256 or any(value not in range(4) for value in correction):
            raise ValueError('A custom trace carrier has 256 F4 entries.')
    values = []
    for entry in range(256):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        value = scale(scalar['endpoint'], int(x == v and y == u))
        value = add(value, scale(scalar['row'], row[entry]))
        value = add(value, scale(scalar['column'], column[entry]))
        value = add(value, scale(scalar['kernel'], correction[entry]))
        values.append(value)
    seed = {'q_x': 4, 'q_u': 4, 'dimensions': [4] * 4, 'degree': max(map(degree, values)),
            'values': values, 'strict_endpoint': True, 'denominator_degree': 0,
            'construction': 'opposite seventh-order Bell CRT jets with a constant trace carrier'}
    assert seed['degree'] == 16
    prescribed = [bell_jets(row_matrix, 0, 2, 7), bell_jets(column_matrix, 1, 3, 7)]
    verify_target(tensor(4), seed, 16, [], 0, 0, 16, 'exact', prescribed)
    if carrier == 'row':
        assert [pair[0] & 1 | ((pair[1] & 1) << 1) for pair in values] == column
    if carrier == 'column':
        assert [pair[0] & 1 | ((pair[1] & 1) << 1) for pair in values] == row
    return seed
