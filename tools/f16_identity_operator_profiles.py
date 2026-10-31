#!/usr/bin/env python3
"""Complete F16 operator profiles by restriction of scalars to native F4.

For F16=F4[u]/(u^2+u+omega), multiplication by a+u*b has matrix
[[a,omega*b],[b,a+b]]. The expanded matrix has every F16 Smith order
twice. Complete Laurent polynomials are retained before the native call.
"""
from __future__ import annotations

from collections import defaultdict

from coupled_identity_operator_filters import (
    _geometry, original_components, states, profile_data as binary_profile_data,
)
from coupled_identity_weyl_profiles import reduced_components, weyl_data
from f4_lattice_linear import MUL as MUL4, f16_multiply
from f4_repaired_identity_filters import laurent_orders as f4_orders
from prime_state_filter_profiles import validate


MUL = tuple(tuple(f16_multiply(a, b) for b in range(16)) for a in range(16))
INV = (0,) + tuple(next(b for b in range(1, 16) if MUL[a][b] == 1) for a in range(1, 16))


def add(first, second, scalar=1):
    result = dict(first)
    for degree, value in second.items():
        result[degree] = result.get(degree, 0) ^ MUL[scalar][value]
    return {degree: value for degree, value in result.items() if value}


def multiply(first, second):
    result = {}
    for a, x in first.items():
        for b, y in second.items():
            result[a+b] = result.get(a+b, 0) ^ MUL[x][y]
    return {degree: value for degree, value in result.items() if value}


def evaluate_one(polynomial):
    result = 0
    for coefficient in polynomial.values():
        result ^= coefficient
    return result


def inverse_matrix(matrix):
    size = len(matrix)
    rows = [list(row)+[int(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column]), None)
        if pivot is None:
            raise ValueError('The F16 endpoint matrix must be invertible.')
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scalar = INV[rows[column][column]]
        rows[column] = [MUL[scalar][value] for value in rows[column]]
        for i in range(size):
            if i != column and rows[i][column]:
                scalar = rows[i][column]
                rows[i] = [a ^ MUL[scalar][b] for a, b in zip(rows[i], rows[column])]
    return tuple(tuple(row[size:]) for row in rows)


def normalized_kernel(modes, kernel, selected):
    alphabet = states(2, modes)
    energies = {state: sum(j*digit for j, digit in zip(selected, state)) for state in alphabet}
    result = {}
    for (row, column), polynomial in kernel.items():
        if row not in energies or column not in energies or sum(row) != sum(column):
            raise ValueError('The F16 operator must preserve valid binary occupation.')
        if any(not isinstance(degree, int) or not isinstance(value, int) or not 0 <= value < 16
               for degree, value in polynomial.items()):
            raise ValueError('Use complete F16 Laurent coefficients encoded by integers 0..15.')
        value = {degree: coefficient for degree, coefficient in polynomial.items() if coefficient}
        if value:
            if evaluate_one(value) and energies[row] < energies[column]:
                raise ValueError('The endpoint must preserve the physical Y energy filtration.')
            result[row, column] = value
    endpoint = [[evaluate_one(result.get((row, column), {})) for column in alphabet] for row in alphabet]
    return result, inverse_matrix(endpoint)


def laurent_orders(matrix):
    size = len(matrix)
    expanded = [[{} for _ in range(2*size)] for _ in range(2*size)]
    for i, row in enumerate(matrix):
        for j, polynomial in enumerate(row):
            for degree, value in polynomial.items():
                a, b = value & 3, value >> 2
                coefficients = ((a, MUL4[2][b]), (b, a ^ b))
                for r in range(2):
                    for c in range(2):
                        if coefficients[r][c]:
                            expanded[2*i+r][2*j+c][degree] = coefficients[r][c]
    orders, precision = f4_orders(expanded)
    if any(orders[2*i] != orders[2*i+1] for i in range(size)):
        raise ArithmeticError('Restriction of scalars did not duplicate every Smith order.')
    return tuple(orders[::2]), precision


def zero_profiles(n, ramification, modes, kernel, geometry='initial'):
    validate(2, n, ramification, (0, 0))
    data = weyl_data(2, n, modes, geometry)
    kernel, inverse = normalized_kernel(modes, kernel, data['selected'])
    alphabet = states(2, modes)
    matrix = [[{} for _ in alphabet] for _ in alphabet]
    for index, scalar in kernel.items():
        for (row, column), polynomial in reduced_components(2, n, modes, geometry)[index].items():
            weight = ramification*sum(j*digit for j, digit in zip(data['selected'], alphabet[row]))
            value = {ramification*degree+weight: coefficient for degree, coefficient in polynomial.items()}
            matrix[row][column] = add(matrix[row][column], multiply(value, scalar))
    reduced = []
    for occupation in range(modes+1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        orders, precision = laurent_orders([[matrix[i][j] for j in indices] for i in indices])
        reduced.append((occupation, tuple(2*value for value in orders)))
    result = defaultdict(list)
    for auxiliary in alphabet:
        shift = -2*ramification*sum(j*digit for j, digit in zip(data['selected'], auxiliary))
        for occupation, orders in reduced:
            result[occupation-sum(auxiliary)].extend(value+shift for value in orders)
    return tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                 for charge, orders in sorted(result.items()))


def infinity_profiles(n, ramification, modes, kernel, geometry='initial'):
    validate(2, n, ramification, (0, 0))
    selected = _geometry(2, n, modes, geometry)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    alphabet = states(2, modes)
    result = []
    for occupation in range(modes+1):
        sector = [state for state in alphabet if sum(state) == occupation]
        matrix = []
        for row in sector:
            values = []
            for column in sector:
                weight = ramification*sum((2*j+1-2*n)*digit for j, digit in zip(selected, column))
                values.append({weight-2*degree: coefficient
                               for degree, coefficient in kernel.get((row, column), {}).items()})
            matrix.append(values)
        result.extend(laurent_orders(matrix)[0])
    return tuple(sorted(result))


def profile_data(n, ramification, charge, modes, kernel, geometry='initial', blocks=None):
    if blocks is None:
        blocks = zero_profiles(n, ramification, modes, kernel, geometry)
    identity = {(state, state): {0: 1} for state in states(2, modes)}
    data = binary_profile_data(2, n, ramification, charge, modes, identity, geometry, blocks)
    data['infinity_active'] = infinity_profiles(n, ramification, modes, kernel, geometry)
    data['coefficient_field'] = 'F16'
    return data


def polynomial_seed(n, ramification, modes, kernel, geometry='initial', repair_endpoint=True):
    validate(2, n, ramification, (0, 0))
    selected = _geometry(2, n, modes, geometry)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    components, gaussian = original_components(2, n, modes, geometry)
    clearing = max(0, -min(degree for value in kernel.values() for degree in value))
    tensor = {}
    for index, scalar in kernel.items():
        scalar = {degree+clearing: value for degree, value in scalar.items()}
        for powers, polynomial in components[index].items():
            value = {ramification*degree: coefficient for degree, coefficient in polynomial.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(value, scalar))
    tensor = {powers: value for powers, value in tensor.items() if value}
    if repair_endpoint:
        alphabet = states(2, modes)
        lookup = {state: i for i, state in enumerate(alphabet)}
        repaired = {}
        for powers, value in tensor.items():
            column = lookup[tuple(powers[n+j] for j in selected)]
            for row, state in enumerate(alphabet):
                if inverse[row][column]:
                    target = list(powers)
                    for j, digit in zip(selected, state):
                        target[n+j] = digit
                    target = tuple(target)
                    repaired[target] = add(repaired.get(target, {}), value, inverse[row][column])
        tensor = {powers: value for powers, value in repaired.items() if value}
        assert all(evaluate_one(tensor.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
                   for powers in set(tensor)|set(gaussian))
    degree = max(degree for value in tensor.values() for degree in value)
    bound = ramification*(max(degree for value in gaussian.values() for degree in value)+2*n*modes)
    bound += max(degree+clearing for value in kernel.values() for degree in value)
    assert degree <= bound
    return tensor, clearing, gaussian
