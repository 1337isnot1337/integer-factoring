#!/usr/bin/env python3
"""Exact prime-field occupation operators on fixed-deficit identity cores.

Complete sparse polynomials feed the existing exact Smith routine. The
auxiliary dummy Gaussian factors are removed before physical costs.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import factorial

import numpy as np
from scipy.special import logsumexp

from coupled_identity_mode_filters import add, multiply, shift, algebra_multiply
from coupled_identity_operator_filters import (
    states, normalized_kernel, endpoint_inverse, _components, _laurent_orders,
)
from deficit_identity_mode_geometry import geometry_data, weyl_data, remove_dummy_factor
from prime_state_filter_pressure import gaussian_pressure
from prime_state_filter_profiles import validate


@lru_cache(maxsize=32)
def parameter(prime, parent, deficit, initial_modes, terminal_modes):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    actual = weyl_data(prime, parent, deficit, initial_modes, terminal_modes)
    if deficit != 1 or initial_modes not in (0, 1):
        return actual
    c, d = {0: 1, 1: prime-1}, {0: 1, parent: prime-1}
    e = add(d, multiply(c, {j: 1 for j in data['augmented']}, prime), prime, -1)
    h = min(e)
    e0 = shift(e, -h, prime)
    delta = multiply(d, e0, prime)
    numerator = tuple(tuple(add(e0 if i == j else {},
        shift(c, parent-h+selected, prime), prime)
        for j in range(data['augmented_modes']))
        for i, selected in enumerate(data['augmented']))
    for i, row in enumerate(numerator):
        for j, value in enumerate(row):
            assert multiply(value, actual['denominator_unit'], prime) == multiply(
                actual['inverse_parameter_numerator'][i][j], delta, prime)
    return {'denominator_unit': delta, 'inverse_parameter_numerator': numerator}


@lru_cache(maxsize=16)
def reduced_components(prime, parent, deficit, initial_modes, terminal_modes):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    local = parameter(prime, parent, deficit, initial_modes, terminal_modes)
    modes, size = data['modes'], data['augmented_modes']
    alphabet, actual_states = states(prime, size), states(prime, modes)
    lookup = {state: i for i, state in enumerate(alphabet)}
    delta, numerator = local['denominator_unit'], local['inverse_parameter_numerator']
    maximum = size*(prime-1)
    delta_powers = [{0: 1}]
    for _ in range(maximum):
        delta_powers.append(multiply(delta_powers[-1], delta, prime))
    forms = []
    for j in range(size):
        form = {}
        for i in range(size):
            powers = [0]*size
            powers[i] = 1
            if numerator[j][i]:
                form[tuple(powers)] = shift(numerator[j][i], 0, prime, -1)
        forms.append(form)
    cache = {(0,)*size: {(0,)*size: {0: 1}}}

    def monomial(powers):
        if powers not in cache:
            j = next(i for i, digit in enumerate(powers) if digit)
            previous = list(powers)
            previous[j] -= 1
            cache[powers] = algebra_multiply(forms[j], monomial(tuple(previous)), prime)
        return cache[powers]

    result = {}
    for a, b in product(actual_states, repeat=2):
        if sum(a) != sum(b):
            continue
        matrix = {}
        for dummy in states(prime, deficit):
            aa, bb = a+dummy, b+dummy
            exponents = tuple(prime-1-digit for digit in bb)
            polynomial = {powers: multiply(value, delta_powers[sum(bb)], prime)
                          for powers, value in monomial(exponents).items()}
            derivatives = tuple(prime-1-digit for digit in aa)
            scalar = (-1)**size
            for digit in aa:
                scalar *= factorial(digit)
            for column, source in enumerate(alphabet):
                for powers, value in polynomial.items():
                    total = tuple(x+y for x, y in zip(source, powers))
                    if any(digit >= prime for digit in total):
                        continue
                    if any(digit < derivative for digit, derivative in zip(total, derivatives)):
                        continue
                    target = tuple(digit-derivative for digit, derivative in zip(total, derivatives))
                    coefficient = scalar
                    for digit, derivative in zip(total, derivatives):
                        coefficient *= factorial(digit)//factorial(digit-derivative)
                    coefficient %= prime
                    if coefficient:
                        entry = (lookup[target], column)
                        matrix[entry] = add(matrix.get(entry, {}), value, prime, coefficient)
        result[a, b] = {entry: value for entry, value in matrix.items() if value}
    return result


def reduced_matrix(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint='identity'):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    validate(prime, data['n'], P, (0,)*prime)
    kernel = normalized_kernel(prime, data['modes'], kernel, endpoint, data['selected'])
    alphabet = states(prime, data['augmented_modes'])
    matrix = [[{} for _ in alphabet] for _ in alphabet]
    components = reduced_components(prime, parent, deficit, initial_modes, terminal_modes)
    for index, scalar in kernel.items():
        for (row, column), polynomial in components[index].items():
            weight = P*sum(j*digit for j, digit in zip(data['augmented'], alphabet[row]))
            value = {P*degree+weight: coefficient for degree, coefficient in polynomial.items()}
            matrix[row][column] = add(matrix[row][column], multiply(value, scalar, prime), prime)
    return matrix


def zero_profiles(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint='identity'):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    size, modes = data['augmented_modes'], data['modes']
    alphabet = states(prime, size)
    matrix = reduced_matrix(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint)
    result = defaultdict(list)
    for occupation in range(size*(prime-1)+1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        orders = tuple(2*value for value in _laurent_orders([[matrix[i][j] for j in indices] for i in indices], prime))
        for auxiliary in alphabet:
            offset = -2*P*sum(j*digit for j, digit in zip(data['augmented'], auxiliary))
            result[occupation-sum(auxiliary)].extend(value+offset for value in orders)
    result = dict(result)
    for j in reversed(data['dummy']):
        result = remove_dummy_factor(result, prime, 2*P*j, min(result), max(result))
    assert sum(map(len, result.values())) == prime**(modes+size)
    return tuple({'charge': charge, 'orders': tuple(sorted(values))}
                 for charge, values in sorted(result.items()))


def infinity_profiles(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint='identity'):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    validate(prime, data['n'], P, (0,)*prime)
    kernel = normalized_kernel(prime, data['modes'], kernel, endpoint, data['selected'])
    alphabet = states(prime, data['modes'])
    orders = []
    for occupation in range(data['modes']*(prime-1)+1):
        labels = [state for state in alphabet if sum(state) == occupation]
        matrix = [[{P*sum((2*j+1-2*data['n'])*digit for j, digit in zip(data['selected'], column))-2*degree:
                    value for degree, value in kernel.get((row, column), {}).items()}
                   for column in labels] for row in labels]
        orders.extend(_laurent_orders(matrix, prime))
    return tuple(sorted(orders))


@lru_cache(maxsize=8)
def active_components(prime, parent, deficit, initial_modes, terminal_modes, P):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    modes, size = data['modes'], data['augmented_modes']
    zero = (0,)*(2*(modes+size))
    gaussian = {zero: {0: 1}}
    for index, selected in enumerate(data['selected']):
        factor = {}
        for digit in range(prime):
            powers = list(zero)
            powers[size+index] = powers[size+modes+index] = digit
            factor[tuple(powers)] = {2*P*selected*digit: pow(factorial(digit), -1, prime)}
        gaussian = algebra_multiply(gaussian, factor, prime)
    for index, retained in enumerate(data['retained']):
        factor, power = {}, {0: 1}
        for digit in range(prime):
            powers = list(zero)
            powers[index] = powers[size+2*modes+index] = digit
            factor[tuple(powers)] = {2*P*(degree+(retained+1)*digit): coefficient*pow(factorial(digit), -1, prime) % prime
                                    for degree, coefficient in power.items()}
            power = multiply(power, data['unit'], prime)
        gaussian = algebra_multiply(gaussian, factor, prime)
    left, right = [], []
    for j, selected in enumerate(data['selected']):
        powers = list(zero)
        powers[size+j] = 1
        first = {tuple(powers): {2*P*selected: 1}}
        powers = list(zero)
        powers[size+modes+j] = 1
        second = {tuple(powers): {0: 1}}
        for i in range(size):
            powers = list(zero)
            powers[i] = 1
            value = multiply({0: 1, 1: prime-1}, data['augmented_v'][i][j], prime)
            first[tuple(powers)] = {2*P*degree: coefficient for degree, coefficient in value.items()}
            powers = list(zero)
            powers[size+2*modes+i] = 1
            second[tuple(powers)] = {2*P*(degree-selected): -coefficient % prime
                                    for degree, coefficient in data['augmented_x'][i][j].items()}
        left.append(first)
        right.append(second)
    return _components(gaussian, left, right, prime)


def full_zero_profiles(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint='identity'):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    validate(prime, data['n'], P, (0,)*prime)
    modes, size = data['modes'], data['augmented_modes']
    kernel = normalized_kernel(prime, modes, kernel, endpoint, data['selected'])
    components = active_components(prime, parent, deficit, initial_modes, terminal_modes, P)
    tensor = {}
    for index, scalar in kernel.items():
        scalar = {2*degree: value for degree, value in scalar.items()}
        for powers, value in components[index].items():
            tensor[powers] = add(tensor.get(powers, {}), multiply(value, scalar, prime), prime)
    result = []
    for charge in range(-size*(prime-1), modes*(prime-1)+1):
        labels = [(y, x) for y, x in product(states(prime, modes), states(prime, size)) if sum(y)-sum(x) == charge]
        matrix = [[tensor.get(v+y+u+x, {}) for u, v in labels] for y, x in labels]
        result.append({'charge': charge, 'orders': _laurent_orders(matrix, prime)})
    return tuple(result)


@lru_cache(maxsize=4)
def original_components(prime, n, selected):
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel
    matrix = pascal_kernel(range(n), 2*n, prime)
    gaussian = {powers: {degree: value for degree, value in enumerate(polynomial) if value}
                for powers, polynomial in direct_kernel(matrix, prime).items()}
    left, right = [], []
    for j in selected:
        first, second = {}, {}
        for index, polynomial in [(n+j, (1,))]+[(i, matrix[i][j]) for i in range(n)]:
            if polynomial:
                powers = [0]*(4*n)
                powers[index] = 1
                first[tuple(powers)] = {degree: value for degree, value in enumerate(polynomial) if value}
        for index, polynomial in [(2*n+j, (1,))]+[(3*n+i, matrix[n+j][n+i]) for i in range(n)]:
            if polynomial:
                powers = [0]*(4*n)
                powers[index] = 1
                second[tuple(powers)] = {degree: value for degree, value in enumerate(polynomial) if value}
        left.append(first)
        right.append(second)
    return _components(gaussian, left, right, prime), gaussian


def polynomial_seed(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint='identity'):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    validate(prime, data['n'], P, (0,)*prime)
    modes, n, selected = data['modes'], data['n'], data['selected']
    kernel = normalized_kernel(prime, modes, kernel, endpoint, selected)
    components, gaussian = original_components(prime, n, selected)
    clearing = max(0, -min(degree for value in kernel.values() for degree in value))
    tensor = {}
    for index, scalar in kernel.items():
        scalar = {degree+clearing: coefficient*(-1)**(degree+clearing) % prime
                  for degree, coefficient in scalar.items()}
        for powers, value in components[index].items():
            ramified = {P*degree: coefficient for degree, coefficient in value.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(ramified, scalar, prime), prime)
    tensor = {powers: value for powers, value in tensor.items() if value}
    if endpoint == 'filtered':
        alphabet = states(prime, modes)
        lookup = {state: i for i, state in enumerate(alphabet)}
        inverse = endpoint_inverse(prime, modes, kernel)
        repaired = {}
        for powers, value in tensor.items():
            column = lookup[tuple(powers[n+j] for j in selected)]
            for row, state in enumerate(alphabet):
                if inverse[row][column]:
                    target = list(powers)
                    for j, digit in zip(selected, state):
                        target[n+j] = digit
                    repaired[tuple(target)] = add(repaired.get(tuple(target), {}), value, prime, inverse[row][column])
        tensor = {powers: value for powers, value in repaired.items() if value}
    endpoint_value = lambda value: sum(coefficient*(-1)**degree for degree, coefficient in value.items()) % prime
    assert all(endpoint_value(tensor.get(powers, {})) == endpoint_value(gaussian.get(powers, {}))
               for powers in set(tensor)|set(gaussian))
    bound = P*(max(degree for value in gaussian.values() for degree in value)+2*n*modes*(prime-1))
    bound += max(degree+clearing for value in kernel.values() for degree in value)
    assert max(degree for value in tensor.values() for degree in value) <= bound
    return tensor, clearing, gaussian


def profile_data(prime, parent, deficit, initial_modes, terminal_modes, P, J, kernel, blocks=None, endpoint='identity'):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    if not isinstance(J, (int, Fraction)):
        raise ValueError('Use an exact rational aggregate charge for physical profile data.')
    if blocks is None:
        blocks = zero_profiles(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint)
    n, selected, retained = data['n'], data['selected'], data['retained']
    return {'zero_active': tuple(sorted(value+P*(1-J)*block['charge'] for block in blocks for value in block['orders'])),
            'infinity_active': infinity_profiles(prime, parent, deficit, initial_modes, terminal_modes, P, kernel, endpoint),
            'zero_spectators': tuple(P*(2*i+1-J) for i in range(n) if i not in selected)
                +tuple(P*(2*i+1+J) for i in range(n) if i not in retained),
            'infinity_spectators': tuple(P*(2*i+1-2*n) for i in range(2*n) if i not in selected),
            'paired_physical_first': tuple(P*(4*i-2*(n-1)+J) for i in range(n)),
            'paired_physical_second': tuple(P*(4*i-2*(n-1)-J) for i in range(n))}


class Pressure:
    def __init__(self, prime, parent, deficit, initial_modes, terminal_modes, P, blocks, infinity):
        data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
        self.prime, self.n, self.P = prime, data['n'], P
        self.selected, self.retained = data['selected'], data['retained']
        self.zero = np.array([value/P for block in blocks for value in block['orders']])
        self.charges = np.array([block['charge'] for block in blocks for value in block['orders']])
        self.infinity = np.array(infinity)/P

    def grid(self, scales, charges, relative=False):
        scales, charges = np.broadcast_arrays(np.asarray(scales), np.asarray(charges))
        zero = self.zero[:, None]+self.charges[:, None]*(1-charges)
        value = logsumexp(-zero*scales, axis=0)+logsumexp(-self.infinity[:, None]*scales, axis=0)
        digits = np.arange(self.prime)[:, None]
        for j in self.selected:
            value -= logsumexp(-digits*(2*j+1-charges)*scales, axis=0)
            value -= logsumexp(-digits*(2*j+1-2*self.n)*scales, axis=0)
        for j in self.retained:
            value -= logsumexp(-digits*(2*j+1+charges)*scales, axis=0)
        if not relative:
            value += np.array([gaussian_pressure(self.prime, self.n, float(J), float(t))
                               for t, J in zip(scales.flat, charges.flat)]).reshape(scales.shape)
        return value
