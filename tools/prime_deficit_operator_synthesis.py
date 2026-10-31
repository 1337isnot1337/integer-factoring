#!/usr/bin/env python3
"""Exact prime-field affine kernel synthesis in finite polynomial row frames.

Only coefficients below an explicit bound enter the constraints. The row
frames themselves are polynomial matrices with invertible constant terms;
no unknown source tail is used to assert a valuation.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from coupled_identity_mode_filters import add
from coupled_identity_operator_filters import states, normalized_kernel, _inverse_matrix, _laurent_orders
from prime_deficit_identity_operators import operator_geometry, reduced_components
from prime_state_filter_profiles import validate


def flag_kernel(prime, selected, basis, pole, gains):
    """B^-1 diag(q^-g) + q^-H(I-B^-1), with endpoint exactly I."""
    alphabet = states(prime, len(selected))
    if len(basis) != len(alphabet) or len(gains) != len(alphabet) or pole < 0 or min(gains) < 0:
        raise ValueError('Use a square occupation frame and nonnegative pole/gain exponents.')
    inverse = _inverse_matrix(basis, prime)
    kernel = {}
    for i, row in enumerate(alphabet):
        for j, column in enumerate(alphabet):
            value = add({-gains[j]: inverse[i][j]}, {-pole: (int(i == j)-inverse[i][j]) % prime}, prime)
            if value:
                kernel[row, column] = value
    return normalized_kernel(prime, len(selected), kernel, 'identity', selected)


def product_below(first, second, prime, precision):
    result = {}
    right = sorted(second.items())
    for a, x in first.items():
        for b, y in right:
            if a+b >= precision:
                break
            result[a+b] = (result.get(a+b, 0)+x*y) % prime
    return {degree: value for degree, value in result.items() if value}


def matrix_product_below(first, second, prime, precision):
    result = [[{} for _ in second[0]] for _ in first]
    for i, row in enumerate(first):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(second[k]):
                    if other:
                        result[i][j] = add(result[i][j], product_below(value, other, prime, precision), prime)
    return result


def regular_row_frame(source, prime, precision):
    """A polynomial unit U with Gaussian row floors equal to exact Smith orders."""
    size = len(source)
    matrix = [[{degree: value for degree, value in entry.items() if degree < precision}
               for entry in row] for row in source]
    left = [[{0: 1} if i == j else {} for j in range(size)] for i in range(size)]
    permutation, orders = list(range(size)), []
    for k in range(size):
        candidates = [(min(matrix[i][j]), i, j) for i in range(k, size)
                      for j in range(k, size) if matrix[i][j]]
        if not candidates:
            raise ArithmeticError('The row frame needs more known coefficients or the source is singular.')
        order, row, column = min(candidates)
        matrix[k], matrix[row] = matrix[row], matrix[k]
        left[k], left[row] = left[row], left[k]
        for line in matrix:
            line[k], line[column] = line[column], line[k]
        permutation[k], permutation[column] = permutation[column], permutation[k]
        unit = {degree-order: value for degree, value in matrix[k][k].items()}
        assert unit.get(0, 0)
        for i in range(k+1, size):
            if not matrix[i][k]:
                continue
            quotient = {degree-order: value for degree, value in matrix[i][k].items()}
            assert min(quotient) >= 0
            for values in (matrix, left):
                values[i] = [add(product_below(unit, a, prime, precision),
                                 product_below(quotient, b, prime, precision), prime, -1)
                             for a, b in zip(values[i], values[k])]
            assert not matrix[i][k]
        orders.append(order)
    _inverse_matrix([[entry.get(0, 0) for entry in row] for row in left], prime)
    transformed = matrix_product_below(left, source, prime, precision)
    assert [[row[j] for j in permutation] for row in transformed] == matrix
    assert tuple(orders) == _laurent_orders(source, prime)
    assert tuple(min(degree for value in row for degree in value) for row in transformed) == tuple(orders)
    return tuple(orders), left


@lru_cache(maxsize=16)
def row_source(prime, parent, selected, P=1, precision_margin=16):
    """Complete input sources, transformed modulo a declared positive precision."""
    geometry = operator_geometry(prime, parent, 1, 0, len(selected), selected)
    validate(prime, geometry['n'], P, (0,)*prime)
    augmented = states(prime, len(selected)+1)
    actual = states(prime, len(selected))
    components = reduced_components(prime, parent, 1, 0, len(selected), selected)
    sources = {}
    for index, entries in components.items():
        matrix = [[{} for _ in augmented] for _ in augmented]
        for (row, column), value in entries.items():
            weight = P*sum(j*digit for j, digit in zip(geometry['augmented'], augmented[row]))
            matrix[row][column] = {P*degree+weight: coefficient for degree, coefficient in value.items()}
        sources[index] = matrix
    gaussian = [[{} for _ in augmented] for _ in augmented]
    for state in actual:
        gaussian = [[add(a, b, prime) for a, b in zip(first, second)]
                    for first, second in zip(gaussian, sources[state, state])]
    groups = []
    maximum = 0
    for occupation in range((len(selected)+1)*(prime-1)+1):
        indices = tuple(i for i, state in enumerate(augmented) if sum(state) == occupation)
        block = [[gaussian[i][j] for j in indices] for i in indices]
        expected = _laurent_orders(block, prime)
        maximum = max(maximum, max(expected))
        groups.append((occupation, indices, block))
    if precision_margin < 2:
        raise ValueError('Use a positive reserve beyond the Gaussian orders.')
    precision = maximum+precision_margin
    transformed, orders, occupations, frames, offset = {index: {} for index in sources}, [], [], [], 0
    for occupation, indices, block in groups:
        block_orders, left = regular_row_frame(block, prime, precision)
        frames.append({'indices': indices, 'row_offset': offset, 'left': left})
        for index, source in sources.items():
            value = matrix_product_below(left, [[source[i][j] for j in indices] for i in indices], prime, precision)
            transformed[index].update({(offset+i, j): entry for i, row in enumerate(value)
                                       for j, entry in enumerate(row) if entry})
        orders.extend(block_orders)
        occupations.extend([occupation]*len(indices))
        offset += len(indices)
    infinity = tuple(P*sum((2*j+1-2*geometry['n'])*digit for j, digit in zip(selected, state))
                     for state in actual)
    return {'prime': prime, 'parent': parent, 'selected': selected, 'ramification': P,
            'states': actual, 'augmented_states': augmented, 'geometry': geometry,
            'sources': transformed, 'zero_orders': tuple(orders), 'row_occupations': tuple(occupations),
            'infinity_orders': infinity, 'known_below': precision, 'frames': tuple(frames)}


@lru_cache(maxsize=16)
def coefficient_system(prime, parent, selected, P=1, lower=-1, upper=1, exponents=None):
    if lower > upper:
        raise ValueError('The Laurent interval is empty.')
    degrees = tuple(range(lower, upper+1)) if exponents is None else tuple(exponents)
    if not degrees or degrees != tuple(sorted(set(degrees))) or min(degrees) < lower or max(degrees) > upper:
        raise ValueError('Use distinct increasing allowed exponents inside the stated Laurent interval.')
    data = row_source(prime, parent, selected, P, max(16, 16-lower))
    variables = tuple((index, degree) for index in data['sources'] for degree in degrees)
    count = len(variables)
    equations = {}
    for variable, (index, degree) in enumerate(variables):
        for (row, column), polynomial in data['sources'][index].items():
            for power, coefficient in polynomial.items():
                key = row, column, degree+power
                if key not in equations:
                    equations[key] = np.zeros(count+1, dtype=np.int64)
                equations[key][variable] = coefficient
    endpoint = []
    for index in data['sources']:
        equation = np.zeros(count+1, dtype=np.int64)
        for variable, (other, degree) in enumerate(variables):
            if other == index:
                equation[variable] = 1
        equation[-1] = int(index[0] == index[1])
        endpoint.append(equation)
    return {'data': data, 'variables': variables, 'zero_constraints': equations,
            'endpoint_constraints': endpoint, 'lower': lower, 'upper': upper,
            'row_prefix_cache': {}}


def homogeneous_row_basis(equations, count, prime):
    pivots = {}
    for source in equations:
        row = np.array(source, dtype=np.int64) % prime
        assert row[-1] == 0
        while True:
            indices = np.flatnonzero(row[:count])
            if not len(indices):
                break
            pivot = int(indices[-1])
            if pivot in pivots:
                row = (row-row[pivot]*pivots[pivot]) % prime
            else:
                pivots[pivot] = row*pow(int(row[pivot]), -1, prime) % prime
                break
    return tuple(pivots.values())


def solve_affine(equations, count, prime):
    pivots = {}
    for source in equations:
        row = np.array(source, dtype=np.int64) % prime
        while True:
            indices = np.flatnonzero(row[:count])
            if not len(indices):
                if row[-1]:
                    return None
                break
            pivot = int(indices[-1])
            if pivot in pivots:
                row = (row-row[pivot]*pivots[pivot]) % prime
            else:
                row = row*pow(int(row[pivot]), -1, prime) % prime
                pivots[pivot] = row
                break
    value = np.zeros(count, dtype=np.int64)
    for pivot in sorted(pivots):
        row = pivots[pivot]
        value[pivot] = (row[-1]-np.dot(row[:count], value)) % prime
    return value, len(pivots)


def constraints(system, zero_targets, infinity_targets, infinity_basis=None):
    data = system['data']
    if len(zero_targets) != len(data['zero_orders']) or len(infinity_targets) != len(data['states']):
        raise ValueError('Use one target for every reduced row at its cut.')
    if any(not isinstance(value, (int, np.integer)) for value in tuple(zero_targets)+tuple(infinity_targets)):
        raise ValueError('Row targets must be integral.')
    if max(zero_targets) > data['known_below']+system['lower']:
        raise ArithmeticError('The requested target exceeds the proved coefficient precision.')
    equations = list(system['endpoint_constraints'])
    count = len(system['variables'])
    for row, target in enumerate(zero_targets):
        key = row, int(target)
        if key not in system['row_prefix_cache']:
            source = (value for (i, column, power), value in system['zero_constraints'].items()
                      if i == row and power < target)
            system['row_prefix_cache'][key] = homogeneous_row_basis(source, count, data['prime'])
        equations.extend(system['row_prefix_cache'][key])
    lookup = {state: i for i, state in enumerate(data['states'])}
    size = len(lookup)
    if infinity_basis is None:
        infinity_basis = tuple(tuple(int(i == j) for j in range(size)) for i in range(size))
    if len(infinity_basis) != size or any(len(row) != size for row in infinity_basis):
        raise ValueError('The infinity row frame has the wrong size.')
    if any(coefficient % data['prime'] and sum(data['states'][i]) != sum(data['states'][j])
           for i, row in enumerate(infinity_basis) for j, coefficient in enumerate(row)):
        raise ValueError('The infinity row frame must preserve occupation.')
    _inverse_matrix(infinity_basis, data['prime'])
    forbidden = {}
    for variable, ((row, column), degree) in enumerate(system['variables']):
        original, j = lookup[row], lookup[column]
        for i in range(size):
            coefficient = infinity_basis[i][original] % data['prime']
            if coefficient and data['infinity_orders'][j]-2*degree < infinity_targets[i]:
                key = i, j, degree
                if key not in forbidden:
                    forbidden[key] = np.zeros(count+1, dtype=np.int64)
                forbidden[key][variable] = coefficient
    equations.extend(forbidden.values())
    return equations


def synthesize(prime, parent, selected, P, lower, upper, zero_targets, infinity_targets, exponents=None, infinity_basis=None):
    system = coefficient_system(prime, parent, selected, P, lower, upper, exponents)
    equations = constraints(system, zero_targets, infinity_targets, infinity_basis)
    answer = solve_affine(equations, len(system['variables']), prime)
    if answer is None:
        return None
    values, rank = answer
    assert all((np.dot(row[:-1], values)-row[-1]) % prime == 0 for row in equations)
    kernel = {}
    for coefficient, (index, degree) in zip(values, system['variables']):
        if coefficient:
            kernel.setdefault(index, {})[degree] = int(coefficient)
    normalized_kernel(prime, len(selected), kernel, 'identity', selected)
    return {'kernel': kernel, 'variables': len(values), 'affine_rank': rank,
            'constraints': len(equations), 'zero_targets': tuple(zero_targets),
            'infinity_targets': tuple(infinity_targets), 'infinity_basis': infinity_basis,
            'known_below': system['data']['known_below']}
