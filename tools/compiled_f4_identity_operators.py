#!/usr/bin/env python3
"""Compile affine F4 operator families into complete native bit-plane matrices."""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

from coupled_identity_operator_filters import states
from coupled_identity_weyl_profiles import reduced_components, weyl_data
from f4_identity_operator_filters import normalized_kernel
from f4_repaired_identity_filters import MUL, add, multiply, evaluate_one
from native_matrix_smith_f4 import MatrixOracle


def _scaled_planes(planes, color):
    first, second = planes
    return ((0, 0), (first, second), (second, first ^ second),
            (first ^ second, first))[color]


def combine_kernel(constant, perturbations, colors):
    result = {index: dict(value) for index, value in constant.items()}
    for perturbation, color in zip(perturbations, colors):
        if color:
            for index, value in perturbation.items():
                result[index] = add(result.get(index, {}), value, color)
    return {index: value for index, value in result.items() if value}


def _compile_matrices(sources):
    size = len(sources[0])
    shifts, entries, bits = [], [], 1
    for i in range(size):
        row_shifts, row_entries = [], []
        for j in range(size):
            minimum = min((degree for matrix in sources for degree in matrix[i][j]), default=0)
            planes = []
            for matrix in sources:
                value = matrix[i][j]
                packed = tuple(sum(1 << (degree - minimum) for degree, coefficient in value.items()
                                   if coefficient & (1 << plane)) for plane in range(2))
                bits = max(bits, *(value.bit_length() for value in packed))
                planes.append(packed)
            row_shifts.append(minimum)
            row_entries.append(tuple(planes))
        shifts.append(row_shifts)
        entries.append(row_entries)
    precision = max(512, 1 << (bits - 1).bit_length())
    if precision > 8192:
        raise ArithmeticError('A complete affine source exceeds the supported native window.')
    return shifts, entries, precision


def _orders(compiled, colors):
    shifts, entries, initial_precision = compiled
    packed = []
    for row in entries:
        values = []
        for sources in row:
            first, second = sources[0]
            for source, color in zip(sources[1:], colors):
                a, b = _scaled_planes(source, color)
                first ^= a
                second ^= b
            values.append((first, second))
        packed.append(values)
    if len(packed) == 1:
        support = packed[0][0][0] | packed[0][0][1]
        if not support:
            raise ArithmeticError('The promised full-rank scalar block vanished.')
        return (shifts[0][0] + (support & -support).bit_length() - 1,)
    precision = initial_precision
    while precision <= 8192:
        try:
            orders = MatrixOracle(packed, len(packed), shifts, precision).profile()['orders']
            if any(order != round(order) for order in orders):
                raise ArithmeticError('A Laurent polynomial matrix returned a noninteger order.')
            return tuple(int(order) for order in orders)
        except ArithmeticError:
            precision *= 2
    raise ArithmeticError('The complete affine matrix did not resolve every promised pivot.')


def compile_family(n, ramification, modes, constant, perturbations, geometry='mixed'):
    """Return evaluate(colors); every call validates the repaired endpoint."""
    alphabet = states(2, modes)
    data = weyl_data(2, n, modes, geometry)
    selected = data['selected']
    components = reduced_components(2, n, modes, geometry)
    zero_sources, infinity_sources = [], []
    endpoint_sources = []
    for kernel in (constant,) + tuple(perturbations):
        zero = [[{} for column in alphabet] for row in alphabet]
        infinity = [[{} for column in alphabet] for row in alphabet]
        for index, scalar in kernel.items():
            if index not in components:
                raise ValueError('Affine perturbations must preserve total binary occupation.')
            if any(not isinstance(degree, int) or not isinstance(value, int) or not 0 <= value < 4
                   for degree, value in scalar.items()):
                raise ValueError('Affine sources must be complete F4 Laurent polynomials.')
            for (row, column), polynomial in components[index].items():
                weight = ramification * sum(j * digit for j, digit in zip(selected, alphabet[row]))
                value = {ramification * degree + weight: coefficient for degree, coefficient in polynomial.items()}
                zero[row][column] = add(zero[row][column], multiply(value, scalar))
        for i, row in enumerate(alphabet):
            for j, column in enumerate(alphabet):
                weight = sum(ramification * (2 * mode + 1 - 2 * n) * digit
                             for mode, digit in zip(selected, column))
                infinity[i][j] = {weight - 2 * degree: coefficient
                                  for degree, coefficient in kernel.get((row, column), {}).items()}
        zero_sources.append(zero)
        infinity_sources.append(infinity)
        endpoint_sources.append(tuple(evaluate_one(kernel.get((row, column), {}))
                                      for row in alphabet for column in alphabet))
    zero_compiled, infinity_compiled = [], []
    for occupation in range(modes + 1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        for sources, result in ((zero_sources, zero_compiled), (infinity_sources, infinity_compiled)):
            result.append((occupation, _compile_matrices(
                [[[matrix[i][j] for j in indices] for i in indices] for matrix in sources])))
    auxiliary = [(sum(state), -2 * ramification * sum(j * digit for j, digit in zip(selected, state)))
                 for state in alphabet]
    endpoint_changes = tuple((i, tuple((j, value) for j, value in enumerate(source) if value))
                             for i, source in enumerate(endpoint_sources[1:]) if any(source))

    @lru_cache(maxsize=4096)
    def validate_endpoint(values):
        endpoint = {(row, column): {0: values[i*len(alphabet)+j]}
                    for i, row in enumerate(alphabet) for j, column in enumerate(alphabet)
                    if values[i*len(alphabet)+j]}
        normalized_kernel(modes, endpoint, selected)
        return True

    def evaluate(colors):
        if len(colors) != len(perturbations) or any(color not in (0, 1, 2, 3) for color in colors):
            raise ValueError('Use one F4 color per affine perturbation.')
        endpoint = list(endpoint_sources[0])
        for i, source in endpoint_changes:
            color = colors[i]
            if color:
                for j, value in source:
                    endpoint[j] ^= MUL[color][value]
        validate_endpoint(tuple(endpoint))
        zero = defaultdict(list)
        for occupation, compiled in zero_compiled:
            orders = tuple(2 * order for order in _orders(compiled, colors))
            for auxiliary_occupation, shift in auxiliary:
                zero[occupation - auxiliary_occupation].extend(order + shift for order in orders)
        infinity = tuple(sorted(order for occupation, compiled in infinity_compiled
                                for order in _orders(compiled, colors)))
        blocks = tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                       for charge, orders in sorted(zero.items()))
        return blocks, infinity

    return evaluate
