#!/usr/bin/env python3
"""Complete F4 occupation-transfer kernels with a filtered endpoint repair."""
from __future__ import annotations

from collections import defaultdict

from coupled_identity_mode_filters import zero_charge_blocks
from coupled_identity_operator_filters import (
    _geometry, active_components, original_components, states,
    profile_data as binary_profile_data,
)
from coupled_identity_weyl_profiles import reduced_components, weyl_data
from f4_repaired_identity_filters import MUL, INV, add, multiply, evaluate_one, laurent_orders
from prime_state_filter_profiles import validate


def inverse_matrix(matrix):
    size = len(matrix)
    rows = [list(row) + [int(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column]), None)
        if pivot is None:
            raise ValueError('The F4 endpoint matrix must be invertible.')
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = INV[rows[column][column]]
        rows[column] = [MUL[inverse][value] for value in rows[column]]
        for i in range(size):
            if i != column and rows[i][column]:
                coefficient = rows[i][column]
                rows[i] = [a ^ MUL[coefficient][b] for a, b in zip(rows[i], rows[column])]
    return tuple(tuple(row[size:]) for row in rows)


def normalized_kernel(modes, kernel, selected):
    alphabet = states(2, modes)
    energies = {state: sum(j * digit for j, digit in zip(selected, state)) for state in alphabet}
    result = {}
    for (row, column), polynomial in kernel.items():
        if row not in energies or column not in energies or sum(row) != sum(column):
            raise ValueError('Use valid binary states and preserve total occupation.')
        if any(not isinstance(degree, int) or not isinstance(value, int) or not 0 <= value < 4
               for degree, value in polynomial.items()):
            raise ValueError('Use integer Laurent powers and F4 coefficients encoded by 0,1,2,3.')
        value = {degree: coefficient for degree, coefficient in polynomial.items() if coefficient}
        if value:
            if evaluate_one(value) and energies[row] < energies[column]:
                raise ValueError('The endpoint must preserve the physical Y energy filtration.')
            result[row, column] = value
    endpoint = [[evaluate_one(result.get((row, column), {})) for column in alphabet] for row in alphabet]
    inverse = inverse_matrix(endpoint)
    return result, inverse


def active_tensor(n, ramification, modes, kernel, geometry='initial'):
    selected = _geometry(2, n, modes, geometry)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    components = active_components(2, n, ramification, modes, geometry)
    result = {}
    for index, scalar in kernel.items():
        scalar = {2 * degree: value for degree, value in scalar.items()}
        for powers, value in components[index].items():
            result[powers] = add(result.get(powers, {}), multiply(value, scalar))
    return {powers: value for powers, value in result.items() if value}


def full_zero_profiles(n, ramification, modes, kernel, geometry='initial'):
    tensor = active_tensor(n, ramification, modes, kernel, geometry)
    result = []
    for charge, matrix in zero_charge_blocks(2, modes, tensor):
        orders, precision = laurent_orders(matrix)
        result.append({'charge': charge, 'orders': orders})
    return tuple(result)


def zero_profiles(n, ramification, modes, kernel, geometry='initial'):
    validate(2, n, ramification, (0, 0))
    data = weyl_data(2, n, modes, geometry)
    kernel, inverse = normalized_kernel(modes, kernel, data['selected'])
    alphabet = states(2, modes)
    matrix = [[{} for column in alphabet] for row in alphabet]
    components = reduced_components(2, n, modes, geometry)
    for index, scalar in kernel.items():
        for (row, column), value in components[index].items():
            ramified = {ramification * degree: coefficient for degree, coefficient in value.items()}
            matrix[row][column] = add(matrix[row][column], multiply(ramified, scalar))
    for row, state in enumerate(alphabet):
        weight = ramification * sum(j * digit for j, digit in zip(data['selected'], state))
        matrix[row] = [{degree + weight: value for degree, value in entry.items()} for entry in matrix[row]]
    reduced = []
    for occupation in range(modes + 1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        orders, precision = laurent_orders([[matrix[i][j] for j in indices] for i in indices])
        reduced.append((occupation, tuple(2 * order for order in orders)))
    result = defaultdict(list)
    for auxiliary in alphabet:
        shift = -2 * ramification * sum(j * digit for j, digit in zip(data['selected'], auxiliary))
        for occupation, orders in reduced:
            result[occupation - sum(auxiliary)].extend(order + shift for order in orders)
    return tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                 for charge, orders in sorted(result.items()))


def infinity_profiles(n, ramification, modes, kernel, geometry='initial'):
    validate(2, n, ramification, (0, 0))
    selected = _geometry(2, n, modes, geometry)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    alphabet = states(2, modes)
    matrix = []
    for row in alphabet:
        values = []
        for column in alphabet:
            weight = sum(ramification * (2 * j + 1 - 2 * n) * digit
                         for j, digit in zip(selected, column))
            values.append({weight - 2 * degree: coefficient
                           for degree, coefficient in kernel.get((row, column), {}).items()})
        matrix.append(values)
    return laurent_orders(matrix)[0]


def profile_data(n, ramification, charge, modes, kernel, geometry='initial', blocks=None):
    if blocks is None:
        blocks = zero_profiles(n, ramification, modes, kernel, geometry)
    identity = {(state, state): {0: 1} for state in states(2, modes)}
    data = binary_profile_data(2, n, ramification, charge, modes, identity, geometry, blocks)
    data['infinity_active'] = infinity_profiles(n, ramification, modes, kernel, geometry)
    data['coefficient_field'] = 'F4'
    return data


def polynomial_seed(n, ramification, modes, kernel, geometry='initial', repair_endpoint=True):
    validate(2, n, ramification, (0, 0))
    selected = _geometry(2, n, modes, geometry)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    components, gaussian = original_components(2, n, modes, geometry)
    clearing = max(0, -min(degree for value in kernel.values() for degree in value))
    tensor = {}
    for index, scalar in kernel.items():
        scalar = {degree + clearing: coefficient for degree, coefficient in scalar.items()}
        for powers, value in components[index].items():
            ramified = {ramification * degree: coefficient for degree, coefficient in value.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(ramified, scalar))
    tensor = {powers: value for powers, value in tensor.items() if value}
    if repair_endpoint:
        alphabet = states(2, modes)
        lookup = {state: i for i, state in enumerate(alphabet)}
        repaired = {}
        for powers, value in tensor.items():
            column = lookup[tuple(powers[n + j] for j in selected)]
            for row, state in enumerate(alphabet):
                if inverse[row][column]:
                    target = list(powers)
                    for j, digit in zip(selected, state):
                        target[n + j] = digit
                    target = tuple(target)
                    repaired[target] = add(repaired.get(target, {}), value, inverse[row][column])
        tensor = {powers: value for powers, value in repaired.items() if value}
        assert all(evaluate_one(tensor.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
                   for powers in set(tensor) | set(gaussian))
    degree = max(degree for value in tensor.values() for degree in value)
    bound = ramification * (max(degree for value in gaussian.values() for degree in value)
                           + 2 * n * modes) + max(degree + clearing for value in kernel.values() for degree in value)
    assert degree <= bound
    return tensor, clearing, gaussian
