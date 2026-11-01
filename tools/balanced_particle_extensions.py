#!/usr/bin/env python3
"""Exact balanced four-particle tensors on separate physical alphabets.

This generalizes the literal q5 construction to paired alphabet sizes
A,B,B,A without altering the frozen q5 constructor or search sources.
"""
from __future__ import annotations

from balanced_four_particle_codes import FOUR_PARTICLE_STATES, REFERENCE_CODE, compress_square, grade, natural_energies
from native_matrix_smith_f4 import MatrixOracle
from saturated_dual_kernel_directions import divide, polynomial_gcd
from selected_gaussian_coefficients import selected_tensor


def canonical(code):
    code = tuple(code)
    if not 2 <= len(code) <= 64 or len(set(code)) != len(code) or any(mask not in FOUR_PARTICLE_STATES for mask in code):
        raise ValueError('Choose between two and64 distinct four-particle masks on eight modes.')
    return tuple(sorted(code, key=lambda mask: (grade(mask), mask)))


def coordinates(index, a, b):
    x = index % a
    index //= a
    y = index % b
    index //= b
    u = index % b
    return x, y, u, index//b


def construct(code_x, code_y=None):
    code_x = canonical(code_x)
    code_y = code_x if code_y is None else canonical(code_y)
    a, b = len(code_x), len(code_y)
    original = selected_tensor(tuple(range(8)), 8, code_x, code_y)
    compressed, common = [], 0
    for low, high in original['values']:
        assert not high
        value = compress_square(low)
        compressed.append(value)
        common = polynomial_gcd(common, value)
    assert common and common.bit_count() % 2 == 1
    values = []
    for index, value in enumerate(compressed):
        quotient, remainder = divide(value, common)
        assert not remainder
        x, y, u, v = coordinates(index, a, b)
        assert quotient.bit_count() % 2 == int(x == v and y == u)
        values.append((quotient, 0))
    return {'values': tuple(values), 'code_x': code_x, 'code_y': code_y,
            'q_x': a, 'q_u': b, 'dimensions': [a, b, b, a],
            'energies_x': natural_energies(code_x), 'energies_y': natural_energies(code_y),
            'entry_gcd_in_z': common, 'entry_gcd_degree': common.bit_length()-1,
            'degree': max(value.bit_length()-1 for value, _ in values),
            'raw_degree_in_z': max(value.bit_length()-1 for value in compressed),
            'strict_endpoint': True, 'denominator_degree': 0,
            'native_complete_degree_bound_in_w': original['coefficient_degree_bound'],
            'native_exact_precision': original['exact_precision']}


def matrices(seed):
    a, b = seed['dimensions'][:2]
    n = a*b
    result = [[[(0, 0)]*n for _ in range(n)] for _ in range(2)]
    for index, pair in enumerate(seed['values']):
        x, y, u, v = coordinates(index, a, b)
        result[0][b*x+u][a*y+v] = tuple(pair)
        result[1][b*x+y][a*u+v] = tuple(pair)
    return result


def complete_profiles(seed):
    a, b = seed['dimensions'][:2]
    hx = tuple(value//2 for value in seed['energies_x'])
    hy = tuple(value//2 for value in seed['energies_y'])
    first, second = matrices(seed)
    infinity = MatrixOracle(first, a*b).profile(
        row_weights=[hy[u] for x in range(a) for u in range(b)],
        column_weights=[hx[v] for y in range(b) for v in range(a)], infinity=True)['orders']
    zero = MatrixOracle(second, a*b).profile(
        row_weights=[hx[x]+hy[y] for x in range(a) for y in range(b)])['orders']
    assert all(abs(value-round(value)) < 1e-8 for value in infinity+zero)
    return tuple(map(round, infinity)), tuple(map(round, zero))


def extension_cases():
    reference = canonical(REFERENCE_CODE)
    yield 'reference', None, reference, reference
    omitted = sorted(set(FOUR_PARTICLE_STATES)-set(reference), key=lambda mask: (grade(mask), mask))
    for mask in omitted:
        extended = canonical(reference+(mask,))
        yield 'symmetric_6x6', mask, extended, extended
        yield 'rectangular_5x6', mask, reference, extended
        yield 'rectangular_6x5', mask, extended, reference


def start_params(seed, scale, energy_x=None, energy_y=None):
    a, b = seed['dimensions'][:2]
    alphabets = (a, b, b, a)
    offsets, stride = [], 0
    for size in alphabets:
        offsets.append(stride)
        stride += size-1
    ex = seed['energies_x'] if energy_x is None else energy_x
    ey = seed['energies_y'] if energy_y is None else energy_y
    result = [0.]*(2*stride)+[float(scale)]
    for (side, leg), energy in {(0, 2): ey, (0, 3): ex, (1, 0): ex, (1, 1): ey}.items():
        for symbol in range(1, alphabets[leg]):
            result[side*stride+offsets[leg]+symbol-1] = (energy[symbol]-energy[0])*scale/2
    return result
