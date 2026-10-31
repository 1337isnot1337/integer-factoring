#!/usr/bin/env python3
"""Binary identity occupation filters with arbitrary nonzero F4 endpoints.

F4 is F2[omega]/(omega^2+omega+1), encoded by 0,1,2,3.
The geometric tensor bases lie in F2; scalar filter coefficients lie in
F4. A diagonal physical-Y map restores the crossing endpoint for free.
Every native Smith call packs every coefficient in both bit planes.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import product

from coupled_identity_mode_filters import zero_charge_blocks
from coupled_identity_operator_filters import (
    _geometry, active_components, original_components,
    profile_data as operator_profile_data,
)
from native_matrix_smith_f4 import MatrixOracle


MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
INV = (0, 1, 3, 2)


def add(first, second, coefficient=1):
    result = dict(first)
    for degree, value in second.items():
        result[degree] = result.get(degree, 0) ^ MUL[coefficient][value]
        if not result[degree]:
            del result[degree]
    return result


def multiply(first, second):
    result = {}
    for i, a in first.items():
        for j, b in second.items():
            result[i + j] = result.get(i + j, 0) ^ MUL[a][b]
            if not result[i + j]:
                del result[i + j]
    return result


def evaluate_one(polynomial):
    result = 0
    for value in polynomial.values():
        result ^= value
    return result


def normalized_filters(modes, alpha):
    states = tuple(product(range(2), repeat=modes))
    if set(alpha) != set(states):
        raise ValueError('Supply one F4 filter for every binary occupation state.')
    result = {}
    for state in states:
        polynomial = alpha[state]
        if any(not isinstance(degree, int) or not isinstance(value, int)
               or not 0 <= value < 4 for degree, value in polynomial.items()):
            raise ValueError('Use integer Laurent powers and F4 coefficients encoded by 0,1,2,3.')
        value = {degree: coefficient for degree, coefficient in polynomial.items() if coefficient}
        if not evaluate_one(value):
            raise ValueError('Every occupation endpoint must be nonzero.')
        result[state] = value
    return result


def monomial_filters(modes, exponents, coefficients):
    states = tuple(product(range(2), repeat=modes))
    if len(exponents) != len(states) or len(coefficients) != len(states):
        raise ValueError('Supply one exponent and coefficient per binary state.')
    return normalized_filters(modes, {state: {-exponent: coefficient}
                                     for state, exponent, coefficient
                                     in zip(states, exponents, coefficients)})


def active_tensor(n, ramification, modes, alpha, geometry='initial'):
    alpha = normalized_filters(modes, alpha)
    components = active_components(2, n, ramification, modes, geometry)
    tensor = {}
    for state, polynomial in alpha.items():
        scalar = {2 * degree: value for degree, value in polynomial.items()}
        for powers, value in components[state, state].items():
            tensor[powers] = add(tensor.get(powers, {}), multiply(value, scalar))
    return {powers: value for powers, value in tensor.items() if value}


def laurent_orders(matrix):
    """Finite exact Smith orders of a full-rank square F4 Laurent matrix."""
    shifts, packed, maximum_bits = [], [], 1
    for row in matrix:
        row_shifts, row_packed = [], []
        for value in row:
            shift = min(value) if value else 0
            planes = tuple(sum(1 << (degree - shift) for degree, coefficient in value.items()
                               if coefficient & (1 << plane)) for plane in range(2))
            maximum_bits = max(maximum_bits, *(plane.bit_length() for plane in planes))
            row_shifts.append(shift)
            row_packed.append(planes)
        shifts.append(row_shifts)
        packed.append(row_packed)
    precision = max(512, 1 << (maximum_bits - 1).bit_length())
    while precision <= 8192:
        try:
            orders = MatrixOracle(packed, len(matrix), shifts, precision).profile()['orders']
            if any(value != round(value) for value in orders):
                raise ArithmeticError('An integer Laurent matrix returned a noninteger Smith order.')
            return tuple(int(value) for value in orders), precision
        except ArithmeticError:
            precision *= 2
    raise ArithmeticError('The complete F4 source or its pivots exceed the supported Smith window.')


def zero_profiles(n, ramification, modes, alpha, geometry='initial'):
    tensor = active_tensor(n, ramification, modes, alpha, geometry)
    result = []
    for charge, matrix in zero_charge_blocks(2, modes, tensor):
        orders, precision = laurent_orders(matrix)
        result.append({'charge': charge, 'orders': orders,
                       'exact_jet_precision': precision, 'coefficient_field': 'F4'})
    return tuple(result)


def reduced_zero_profiles(n, ramification, modes, alpha, geometry='initial'):
    """The same charge profiles through the smaller exact Weyl matrix."""
    from coupled_identity_weyl_profiles import reduced_components, weyl_data

    alpha = normalized_filters(modes, alpha)
    states = tuple(alpha)
    data = weyl_data(2, n, modes, geometry)
    components = reduced_components(2, n, modes, geometry)
    matrix = [[{} for _ in states] for _ in states]
    for state, scalar in alpha.items():
        for (row, column), polynomial in components[state, state].items():
            value = {ramification * degree: coefficient for degree, coefficient in polynomial.items()}
            matrix[row][column] = add(matrix[row][column], multiply(value, scalar))
    for row, state in enumerate(states):
        weight = ramification * sum(j * digit for j, digit in zip(data['selected'], state))
        matrix[row] = [{degree + weight: value for degree, value in polynomial.items()}
                       for polynomial in matrix[row]]
    reduced, precision = [], 0
    for occupation in range(modes + 1):
        indices = [i for i, state in enumerate(states) if sum(state) == occupation]
        orders, used = laurent_orders([[matrix[i][j] for j in indices] for i in indices])
        reduced.append((occupation, tuple(2 * value for value in orders)))
        precision = max(precision, used)
    result = defaultdict(list)
    for auxiliary in states:
        shift = -2 * ramification * sum(j * digit for j, digit in zip(data['selected'], auxiliary))
        for occupation, orders in reduced:
            result[occupation - sum(auxiliary)].extend(value + shift for value in orders)
    return tuple({'charge': charge, 'orders': tuple(sorted(orders)),
                  'exact_jet_precision': precision, 'coefficient_field': 'F4',
                  'backend': 'Weyl reduced exact matrix'}
                 for charge, orders in sorted(result.items()))


def compiled_monomial_profiles(n, ramification, modes, exponents, geometry='initial'):
    """Compile one exponent table; vary its nonzero F4 colors by bit XORs.

    Each state coefficient multiplies an F2 polynomial. Therefore its
    two bit planes can be assembled without any truncated series or
    repeated polynomial convolution. This is the same complete Weyl
    matrix used by reduced_zero_profiles.
    """
    from coupled_identity_weyl_profiles import reduced_components, weyl_data

    states = tuple(product(range(2), repeat=modes))
    monomial_filters(modes, exponents, (1,) * len(states))
    data = weyl_data(2, n, modes, geometry)
    components = reduced_components(2, n, modes, geometry)
    entries = {}
    for state_index, (state, exponent) in enumerate(zip(states, exponents)):
        for (row, column), polynomial in components[state, state].items():
            weight = ramification * sum(j * digit for j, digit in zip(data['selected'], states[row]))
            entries.setdefault((row, column), {})[state_index] = tuple(
                ramification * degree - exponent + weight
                for degree, coefficient in polynomial.items() if coefficient)
    compiled = []
    for occupation in range(modes + 1):
        indices = [i for i, state in enumerate(states) if sum(state) == occupation]
        shifts, sources, maximum_bits = [], [], 1
        for row in indices:
            row_shifts, row_sources = [], []
            for column in indices:
                value = entries.get((row, column), {})
                minimum = min((degree for polynomial in value.values() for degree in polynomial), default=0)
                bits = tuple((state, sum(1 << (degree - minimum) for degree in polynomial))
                             for state, polynomial in value.items())
                maximum_bits = max([maximum_bits] + [value.bit_length() for _, value in bits])
                row_shifts.append(minimum)
                row_sources.append(bits)
            shifts.append(row_shifts)
            sources.append(row_sources)
        compiled.append((occupation, shifts, sources, max(512, 1 << (maximum_bits - 1).bit_length())))
    auxiliary_shifts = tuple((sum(state), -2 * ramification * sum(
        j * digit for j, digit in zip(data['selected'], state))) for state in states)

    def evaluate(coefficients):
        if len(coefficients) != len(states) or any(value not in (1, 2, 3) for value in coefficients):
            raise ValueError('Use one nonzero F4 coefficient per occupation state.')
        result, maximum_precision = defaultdict(list), 0
        for occupation, shifts, sources, initial_precision in compiled:
            matrix = []
            for row in sources:
                values = []
                for entry in row:
                    first = second = 0
                    for state, bits in entry:
                        if coefficients[state] & 1:
                            first ^= bits
                        if coefficients[state] & 2:
                            second ^= bits
                    values.append((first, second))
                matrix.append(values)
            precision = initial_precision
            while precision <= 8192:
                try:
                    profile = MatrixOracle(matrix, len(matrix), shifts, precision).profile()['orders']
                    if any(value != round(value) for value in profile):
                        raise ArithmeticError('An integer Laurent matrix returned a noninteger order.')
                    orders = tuple(2 * int(value) for value in profile)
                    break
                except ArithmeticError:
                    precision *= 2
            else:
                raise ArithmeticError('The complete F4 matrix exceeds the supported Smith window.')
            maximum_precision = max(maximum_precision, precision)
            for auxiliary_occupation, shift in auxiliary_shifts:
                result[occupation - auxiliary_occupation].extend(value + shift for value in orders)
        return tuple({'charge': charge, 'orders': tuple(sorted(orders)),
                      'exact_jet_precision': maximum_precision, 'coefficient_field': 'F4',
                      'backend': 'Compiled complete Weyl bit planes'}
                     for charge, orders in sorted(result.items()))

    return evaluate


def profile_data(n, ramification, charge, modes, alpha, geometry='initial', blocks=None):
    alpha = normalized_filters(modes, alpha)
    if blocks is None:
        blocks = zero_profiles(n, ramification, modes, alpha, geometry)
    identity = {(state, state): {0: 1} for state in alpha}
    result = operator_profile_data(2, n, ramification, charge, modes,
                                   identity, geometry, blocks)
    selected = _geometry(2, n, modes, geometry)['selected']
    result['infinity_active'] = tuple(sorted(
        sum(ramification * (2 * j + 1 - 2 * n) * digit
            for j, digit in zip(selected, state)) - 2 * max(polynomial)
        for state, polynomial in alpha.items()))
    result['coefficient_field'] = 'F4'
    result['endpoint_values'] = tuple(evaluate_one(value) for value in alpha.values())
    return result


def polynomial_seed(n, ramification, modes, alpha, geometry='initial', repair_endpoint=True):
    """Complete original physical tensor, for independent small controls."""
    alpha = normalized_filters(modes, alpha)
    selected = _geometry(2, n, modes, geometry)['selected']
    components, gaussian = original_components(2, n, modes, geometry)
    clearing = max(0, -min(degree for value in alpha.values() for degree in value))
    tensor = {}
    for state, polynomial in alpha.items():
        scalar = {degree + clearing: value for degree, value in polynomial.items()}
        for powers, value in components[state, state].items():
            ramified = {ramification * degree: coefficient for degree, coefficient in value.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(ramified, scalar))
    tensor = {powers: value for powers, value in tensor.items() if value}
    if repair_endpoint:
        endpoints = {state: evaluate_one(value) for state, value in alpha.items()}
        tensor = {powers: {degree: MUL[coefficient][
            INV[endpoints[tuple(powers[n + j] for j in selected)]]]
                           for degree, coefficient in value.items()}
                  for powers, value in tensor.items()}
    return tensor, clearing, gaussian
