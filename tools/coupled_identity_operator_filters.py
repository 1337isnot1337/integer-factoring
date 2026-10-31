#!/usr/bin/env python3
"""Endpoint-identity operator kernels on selected polynomial identity modes.

Unlike a diagonal occupation filter, the selected kernel may transfer
occupation between modes. This implementation keeps equal total input
and output occupation, so the existing zero charge blocks remain valid.
All scalar coefficients are Laurent polynomials in q=-v.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
from math import factorial

from coupled_identity_mode_filters import (
    add, multiply, algebra_multiply, concentration_data,
    zero_charge_blocks, active_power_kernels, profile_data as diagonal_profile_data,
)
from prime_state_filter_profiles import validate


def states(prime, modes):
    return tuple(product(range(prime), repeat=modes))


def identity_kernel(prime, modes):
    return {(state, state): {0: 1} for state in states(prime, modes)}


def _inverse_matrix(matrix, prime):
    size = len(matrix)
    rows = [[value % prime for value in row] + [int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column]), None)
        if pivot is None:
            raise ValueError('The endpoint kernel must be invertible.')
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = pow(rows[column][column], -1, prime)
        rows[column] = [inverse * value % prime for value in rows[column]]
        for i in range(size):
            if i != column and rows[i][column]:
                factor = rows[i][column]
                rows[i] = [(a - factor * b) % prime for a, b in zip(rows[i], rows[column])]
    return tuple(tuple(row[size:]) for row in rows)


def endpoint_inverse(prime, modes, kernel):
    alphabet = states(prime, modes)
    matrix = [[sum(kernel.get((row, column), {}).values()) % prime for column in alphabet]
              for row in alphabet]
    return _inverse_matrix(matrix, prime)


def normalized_kernel(prime, modes, kernel, endpoint='identity', selected=None):
    alphabet = set(states(prime, modes))
    result = {}
    for (row, column), polynomial in kernel.items():
        if row not in alphabet or column not in alphabet:
            raise ValueError('Every kernel index must be a valid occupation state.')
        value = {power: coefficient % prime for power, coefficient in polynomial.items()
                 if coefficient % prime}
        if any(not isinstance(power, int) for power in value):
            raise ValueError('Use integer Laurent powers in q=-v.')
        if value and sum(row) != sum(column):
            raise ValueError('This operator family preserves total occupation.')
        if value:
            result[row, column] = value
    if endpoint == 'identity':
        for row in alphabet:
            for column in alphabet:
                if sum(result.get((row, column), {}).values()) % prime != int(row == column):
                    raise ValueError('The selected operator kernel must equal the identity at q=1.')
    elif endpoint == 'filtered':
        if selected is None:
            selected = tuple(range(modes))
        if len(selected) != modes:
            raise ValueError('One physical mode index is required per selected direction.')
        energy = {state: sum(j * digit for j, digit in zip(selected, state)) for state in alphabet}
        for (row, column), value in result.items():
            if sum(value.values()) % prime and energy[row] < energy[column]:
                raise ValueError('The endpoint kernel must preserve the physical Y energy filtration.')
        endpoint_inverse(prime, modes, result)
    else:
        raise ValueError('Endpoint policy must be identity or filtered.')
    return result


def _geometry(prime, n, modes, geometry):
    if geometry == 'initial':
        data = concentration_data(prime, n, modes)
        return {'selected': tuple(range(modes)), 'retained': data['top'],
                'active_x': data['active_x'], 'active_v': data['active_v']}
    from coupled_identity_extreme_filters import geometry_data
    return geometry_data(prime, n, modes, geometry)


def _scaled_algebra(value, scalar, prime):
    return {powers: {degree: coefficient * scalar % prime
                     for degree, coefficient in polynomial.items()
                     if coefficient * scalar % prime}
            for powers, polynomial in value.items() if scalar % prime}


def _algebra_add(first, second, prime, coefficient=1):
    result = dict(first)
    for powers, polynomial in second.items():
        value = add(result.get(powers, {}), polynomial, prime, coefficient)
        if value:
            result[powers] = value
        else:
            result.pop(powers, None)
    return result


def _removed_gaussian(gaussian, left, right, prime):
    result = gaussian
    zero = (0,) * len(next(iter(gaussian)))
    for first, second in zip(left, right):
        q = algebra_multiply(first, second, prime)
        exponential = {zero: {0: 1}}
        power = exponential
        for digit in range(1, prime):
            power = algebra_multiply(power, q, prime)
            exponential = _algebra_add(exponential, power, prime,
                                       (-1)**digit * pow(factorial(digit), -1, prime))
        result = algebra_multiply(result, exponential, prime)
    return result


def _components(gaussian, left, right, prime):
    modes = len(left)
    forms = tuple(left) + tuple(right)
    cache = {(0,) * (2 * modes): _removed_gaussian(gaussian, left, right, prime)}
    def component(powers):
        if powers not in cache:
            index = next(i for i, value in enumerate(powers) if value)
            previous = list(powers)
            previous[index] -= 1
            cache[powers] = algebra_multiply(forms[index], component(tuple(previous)), prime)
        return cache[powers]
    result = {}
    for row in states(prime, modes):
        for column in states(prime, modes):
            if sum(row) != sum(column):
                continue
            normalization = 1
            for digit in column:
                normalization = normalization * pow(factorial(digit), -1, prime) % prime
            result[row, column] = _scaled_algebra(component(row + column), normalization, prime)
    return result


@lru_cache(maxsize=16)
def active_components(prime, n, ramification=1, modes=2, geometry='initial'):
    validate(prime, n, ramification, (0,) * prime)
    data = _geometry(prime, n, modes, geometry)
    if geometry == 'initial':
        gaussian = active_power_kernels(prime, n, ramification, modes)[(0,) * modes]
    else:
        from coupled_identity_extreme_filters import active_power_kernels as extreme_kernels
        gaussian = extreme_kernels(prime, n, ramification, modes, geometry)[(0,) * modes]
    left, right = [], []
    zero = (0,) * (4 * modes)
    for j, selected in enumerate(data['selected']):
        order = 2 * ramification * selected
        powers = list(zero)
        powers[modes + j] = 1
        first = {tuple(powers): {order: 1}}
        powers = list(zero)
        powers[2 * modes + j] = 1
        second = {tuple(powers): {0: 1}}
        for h in range(modes):
            powers = list(zero)
            powers[h] = 1
            coefficient = multiply({0: 1, 1: -1 % prime}, data['active_v'][h][j], prime)
            first[tuple(powers)] = {2 * ramification * degree: value
                                    for degree, value in coefficient.items()}
            powers = list(zero)
            powers[3 * modes + h] = 1
            second[tuple(powers)] = {2 * ramification * degree - order: -value % prime
                                     for degree, value in data['active_x'][h][j].items()}
        left.append(first)
        right.append(second)
    return _components(gaussian, left, right, prime)


def active_tensor(prime, n, ramification, modes, kernel, geometry='initial', endpoint='identity'):
    kernel = normalized_kernel(prime, modes, kernel, endpoint,
                               _geometry(prime, n, modes, geometry)['selected'])
    components = active_components(prime, n, ramification, modes, geometry)
    result = {}
    for index, scalar in kernel.items():
        scalar = {2 * power: value for power, value in scalar.items()}
        for variables, polynomial in components[index].items():
            value = multiply(polynomial, scalar, prime)
            result[variables] = add(result.get(variables, {}), value, prime)
    return {variables: polynomial for variables, polynomial in result.items() if polynomial}


def _laurent_orders(matrix, prime):
    if prime == 2:
        from native_matrix_smith_f4 import MatrixOracle
        shifts = [[min(value) if value else 0 for value in row] for row in matrix]
        packed = [[(sum(1 << (degree - shifts[i][j]) for degree, value in entry.items() if value), 0)
                   for j, entry in enumerate(row)] for i, row in enumerate(matrix)]
        bits = max(1, max(value[0].bit_length() for row in packed for value in row))
        precision = max(512, 1 << (bits - 1).bit_length())
        while precision <= 8192:
            try:
                orders = MatrixOracle(packed, len(matrix), shifts, precision).profile()['orders']
                assert all(order == round(order) for order in orders)
                return tuple(int(order) for order in orders)
            except ArithmeticError:
                precision *= 2
        raise ArithmeticError('The complete Laurent matrix exceeds the supported Smith window.')
    from check_all_field_pascal_crossings import smith_orders
    minimum = min(degree for row in matrix for value in row for degree in value)
    maximum = max(degree for row in matrix for value in row for degree in value)
    dense = [[tuple(value.get(degree, 0) for degree in range(minimum, max(value) + 1))
              if value else () for value in row] for row in matrix]
    orders, precision = smith_orders(dense, prime,
                                     maximum_precision=max(512, 2 * len(matrix) * (maximum - minimum + 1)))
    return tuple(order + minimum for order in orders)


def zero_profiles(prime, n, ramification, modes, kernel, geometry='initial', endpoint='identity'):
    tensor = active_tensor(prime, n, ramification, modes, kernel, geometry, endpoint)
    return tuple({'charge': charge, 'orders': _laurent_orders(matrix, prime)}
                 for charge, matrix in zero_charge_blocks(prime, modes, tensor))


def infinity_profiles(prime, n, ramification, modes, kernel, geometry='initial', endpoint='identity'):
    selected = _geometry(prime, n, modes, geometry)['selected']
    kernel = normalized_kernel(prime, modes, kernel, endpoint, selected)
    alphabet = states(prime, modes)
    matrix = []
    for row in alphabet:
        values = []
        for column in alphabet:
            weight = -n * ramification * sum(row)
            weight += sum(ramification * (2 * j - n + 1) * digit
                          for j, digit in zip(selected, column))
            values.append({weight - 2 * degree: coefficient
                           for degree, coefficient in kernel.get((row, column), {}).items()})
        matrix.append(values)
    return _laurent_orders(matrix, prime)


def profile_data(prime, n, ramification, charge, modes, kernel, geometry='initial', blocks=None,
                 endpoint='identity'):
    if blocks is None:
        blocks = zero_profiles(prime, n, ramification, modes, kernel, geometry, endpoint)
    diagonal = {state: {0: 1} for state in states(prime, modes)}
    if geometry == 'initial':
        data = diagonal_profile_data(prime, n, ramification, charge, modes, diagonal, blocks)
    else:
        from coupled_identity_extreme_filters import profile_data as extreme_data
        data = extreme_data(prime, n, ramification, charge, modes, diagonal, geometry, blocks)
    data['infinity_active'] = infinity_profiles(prime, n, ramification, modes, kernel, geometry, endpoint)
    data['operator_kernel'] = kernel
    data['endpoint_policy'] = endpoint
    return data


@lru_cache(maxsize=8)
def original_components(prime, n, modes, geometry='initial'):
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel
    data = _geometry(prime, n, modes, geometry)
    matrix = pascal_kernel(range(n), 2 * n, prime)
    gaussian = {variables: {degree: value for degree, value in enumerate(polynomial) if value}
                for variables, polynomial in direct_kernel(matrix, prime).items()}
    left, right = [], []
    for j in data['selected']:
        first, second = {}, {}
        for index, polynomial in [(n + j, (1,))] + [(i, matrix[i][j]) for i in range(n)]:
            if polynomial:
                variables = [0] * (4 * n)
                variables[index] = 1
                first[tuple(variables)] = {degree: value for degree, value in enumerate(polynomial) if value}
        for index, polynomial in [(2 * n + j, (1,))] + [(3 * n + i, matrix[n + j][n + i]) for i in range(n)]:
            if polynomial:
                variables = [0] * (4 * n)
                variables[index] = 1
                second[tuple(variables)] = {degree: value for degree, value in enumerate(polynomial) if value}
        left.append(first)
        right.append(second)
    return _components(gaussian, left, right, prime), gaussian


def polynomial_seed(prime, n, ramification, modes, kernel, geometry='initial', endpoint='identity'):
    selected = _geometry(prime, n, modes, geometry)['selected']
    kernel = normalized_kernel(prime, modes, kernel, endpoint, selected)
    components, gaussian = original_components(prime, n, modes, geometry)
    clearing = max(0, -min(degree for value in kernel.values() for degree in value))
    result = {}
    for index, scalar in kernel.items():
        scalar = {degree + clearing: coefficient * (-1)**(degree + clearing) % prime
                  for degree, coefficient in scalar.items()}
        for variables, value in components[index].items():
            ramified = {ramification * degree: coefficient for degree, coefficient in value.items()}
            result[variables] = add(result.get(variables, {}), multiply(ramified, scalar, prime), prime)
    result = {variables: value for variables, value in result.items() if value}
    if endpoint == 'filtered':
        alphabet = states(prime, modes)
        lookup = {state: i for i, state in enumerate(alphabet)}
        inverse = endpoint_inverse(prime, modes, kernel)
        repaired = {}
        for variables, value in result.items():
            column = lookup[tuple(variables[n + j] for j in selected)]
            for i, row in enumerate(alphabet):
                if not inverse[i][column]:
                    continue
                target = list(variables)
                for j, digit in zip(selected, row):
                    target[n + j] = digit
                target = tuple(target)
                repaired[target] = add(repaired.get(target, {}), value, prime, inverse[i][column])
        result = {variables: value for variables, value in repaired.items() if value}
    gaussian_degree = max(degree for value in gaussian.values() for degree in value)
    coefficient_degree = max(degree + clearing for value in kernel.values() for degree in value)
    bound = ramification * (gaussian_degree + 2 * n * modes * (prime - 1)) + coefficient_degree
    assert max(degree for value in result.values() for degree in value) <= bound
    endpoint = lambda value: sum(coefficient * (-1)**degree for degree, coefficient in value.items()) % prime
    assert all(endpoint(result.get(variables, {})) == endpoint(gaussian.get(variables, {}))
               for variables in set(result) | set(gaussian))
    return result, clearing
