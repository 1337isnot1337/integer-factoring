#!/usr/bin/env python3
"""Terminal binary operator cores at n=2^k−1, with one extra XV mode."""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import comb

import numpy as np

from coupled_identity_mode_filters import add as add2, multiply as multiply2, shift, algebra_multiply
from coupled_identity_operator_filters import _components, states
from coupled_identity_weyl_profiles import adjugate, determinant
from f4_identity_operator_filters import normalized_kernel
from f4_repaired_identity_filters import add, multiply, evaluate_one, laurent_orders
from search_coupled_identity_operator_filters import Pressure


@lru_cache(maxsize=16)
def geometry_data(n, modes, initial_modes=0):
    parent = n+1
    if parent <= 1 or parent & (parent-1) or not 1 <= modes < n:
        raise ValueError('Use n=2^k−1 and between one and n−1 selected modes.')
    if initial_modes not in (0, 1) or initial_modes >= modes:
        raise ValueError('Use zero or one initial mode and at least one terminal mode.')
    auxiliary = modes+1
    terminal_modes = modes-initial_modes
    selected = tuple(range(initial_modes))+tuple(range(n-terminal_modes, n))
    augmented = selected+(n,)
    retained = tuple(range(terminal_modes+1))+tuple(range(n-initial_modes, n))
    discarded = tuple(i for i in range(n) if i not in retained)

    def form(i, j):
        # Parent-length reversal: f_j(i)=tau^i f_(n-j)(i,(1+tau)/tau).
        if j == 0:
            return {0: 1}
        h = n-j
        result = {}
        for a in range(min(i, h)+1):
            for degree in range(a+1):
                coefficient = comb(h, a)*comb(i, a)*comb(a, degree) % 2
                if coefficient:
                    exponent = i-a+degree
                    result[exponent] = result.get(exponent, 0) ^ 1
        return {degree: value for degree, value in result.items() if value}

    forms = tuple(tuple(form(i, j) for j in augmented) for i in range(n))
    retained_matrix = tuple(forms[i] for i in retained)
    unit = determinant(retained_matrix, 2)
    assert unit.get(0) == 1
    inverse_numerator = adjugate(retained_matrix, 2)
    interpolation, margins = [], []
    for i in discarded:
        row = []
        for a, t in enumerate(retained):
            value = {}
            for j in range(auxiliary):
                value = add2(value, multiply2(forms[i][j], inverse_numerator[j][a], 2), 2)
            assert not value or min(value) >= i-t
            margins.append(None if not value else min(value)-max(0, i-t))
            row.append(value)
        interpolation.append(tuple(row))
    active_v = []
    for a, t in enumerate(retained):
        row = []
        for j in range(auxiliary):
            value = multiply2(unit, retained_matrix[a][j], 2)
            for k, i in enumerate(discarded):
                term = multiply2(interpolation[k][a], forms[i][j], 2)
                value = add2(value, shift(term, t-i, 2), 2)
            assert not value or min(value) >= 0
            row.append(value)
        active_v.append(tuple(row))
    assert determinant(active_v, 2).get(0) == 1
    return {'n': n, 'modes': modes, 'initial_modes': initial_modes, 'auxiliary_modes': auxiliary,
            'selected': selected, 'augmented': augmented, 'retained': retained,
            'discarded': discarded, 'unit': unit, 'forms': forms,
            'interpolation_numerators': tuple(interpolation), 'weighted_order_margins': tuple(margins),
            'augmented_x': retained_matrix, 'augmented_v': tuple(active_v),
            'active_x': tuple(tuple(row[:modes]) for row in retained_matrix),
            'active_v': tuple(tuple(row[:modes]) for row in active_v)}


@lru_cache(maxsize=8)
def active_components(n, ramification=1, modes=2, initial_modes=0):
    if ramification <= 0 or ramification % 2 == 0:
        raise ValueError('Use a positive odd ramification.')
    data = geometry_data(n, modes, initial_modes)
    auxiliary = data['auxiliary_modes']
    zero = (0,)*(2*(modes+auxiliary))
    gaussian = {zero: {0: 1}}
    for j, selected in enumerate(data['selected']):
        powers = list(zero)
        powers[auxiliary+j] = powers[auxiliary+modes+j] = 1
        gaussian = algebra_multiply(gaussian, {zero: {0: 1}, tuple(powers): {2*ramification*selected: 1}}, 2)
    for j in range(auxiliary):
        powers = list(zero)
        powers[j] = powers[auxiliary+2*modes+j] = 1
        coefficient = {2*ramification*(degree+data['retained'][j]+1): value for degree, value in data['unit'].items()}
        gaussian = algebra_multiply(gaussian, {zero: {0: 1}, tuple(powers): coefficient}, 2)
    left, right = [], []
    for j, selected in enumerate(data['selected']):
        powers = list(zero)
        powers[auxiliary+j] = 1
        first = {tuple(powers): {2*ramification*selected: 1}}
        powers = list(zero)
        powers[auxiliary+modes+j] = 1
        second = {tuple(powers): {0: 1}}
        for i in range(auxiliary):
            powers = list(zero)
            powers[i] = 1
            value = multiply2({0: 1, 1: 1}, data['active_v'][i][j], 2)
            first[tuple(powers)] = {2*ramification*degree: coefficient for degree, coefficient in value.items()}
            powers = list(zero)
            powers[auxiliary+2*modes+i] = 1
            second[tuple(powers)] = {2*ramification*(degree-selected): coefficient
                                    for degree, coefficient in data['active_x'][i][j].items()}
        left.append(first)
        right.append(second)
    return _components(gaussian, left, right, 2)


def active_tensor(n, ramification, modes, kernel, initial_modes=0):
    kernel, inverse = normalized_kernel(modes, kernel, geometry_data(n, modes, initial_modes)['selected'])
    result = {}
    components = active_components(n, ramification, modes, initial_modes)
    for index, scalar in kernel.items():
        scalar = {2*degree: value for degree, value in scalar.items()}
        for powers, polynomial in components[index].items():
            result[powers] = add(result.get(powers, {}), multiply(polynomial, scalar))
    return {powers: polynomial for powers, polynomial in result.items() if polynomial}


def charge_blocks(modes, tensor):
    auxiliary = modes+1
    selected_states, auxiliary_states = states(2, modes), states(2, auxiliary)
    result = []
    for charge in range(-auxiliary, modes+1):
        labels = [(y, x) for y, x in product(selected_states, auxiliary_states) if sum(y)-sum(x) == charge]
        result.append((charge, [[tensor.get(v+y+u+x, {}) for u, v in labels] for y, x in labels]))
    return result


def zero_profiles(n, ramification, modes, kernel, initial_modes=0):
    tensor = active_tensor(n, ramification, modes, kernel, initial_modes)
    return tuple({'charge': charge, 'orders': laurent_orders(matrix)[0]}
                 for charge, matrix in charge_blocks(modes, tensor))


def infinity_profiles(n, ramification, modes, kernel, initial_modes=0):
    selected = geometry_data(n, modes, initial_modes)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    alphabet = states(2, modes)
    matrix = []
    for row in alphabet:
        entries = []
        for column in alphabet:
            weight = ramification*sum((2*j+1-2*n)*digit for j, digit in zip(selected, column))
            entries.append({weight-2*degree: value for degree, value in kernel.get((row, column), {}).items()})
        matrix.append(entries)
    return laurent_orders(matrix)[0]


def profile_data(n, ramification, charge, modes, kernel, blocks=None, initial_modes=0):
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an exact rational aggregate charge.')
    if ramification <= 0 or ramification % 2 == 0:
        raise ValueError('Use a positive odd ramification.')
    if blocks is None:
        from near_power_terminal_weyl import zero_profiles as reduced_profiles
        blocks = reduced_profiles(n, ramification, modes, kernel, initial_modes)
    data = geometry_data(n, modes, initial_modes)
    selected, retained = data['selected'], data['retained']
    return {
        'prime': 2, 'coefficient_field': 'F4', 'n': n, 'ramification': ramification,
        'modes': modes, 'charge': charge, 'identity_indices': selected,
        'retained_xv_indices': retained,
        'zero_active': tuple(sorted(value+ramification*(1-charge)*block['charge']
                                    for block in blocks for value in block['orders'])),
        'zero_charge_blocks_at_J1': blocks,
        'zero_spectators': tuple(ramification*(2*i+1-charge) for i in range(n) if i not in selected)
            + tuple(ramification*(2*i+1+charge) for i in range(n) if i not in retained),
        'infinity_active': infinity_profiles(n, ramification, modes, kernel, initial_modes),
        'infinity_spectators': tuple(ramification*(2*i+1-2*n) for i in range(2*n) if i not in selected),
        'paired_physical_first': tuple(ramification*(4*i-2*(n-1)+charge) for i in range(n)),
        'paired_physical_second': tuple(ramification*(4*i-2*(n-1)-charge) for i in range(n)),
    }


class TerminalPressure(Pressure):
    def __init__(self, n, ramification, modes, blocks, infinity, initial_modes=0):
        data = geometry_data(n, modes, initial_modes)
        self.n, self.ramification, self.modes = n, ramification, modes
        self.selected, self.retained = data['selected'], data['retained']
        self.zero = np.array([value/ramification for block in blocks for value in block['orders']])
        self.charges = np.array([block['charge'] for block in blocks for value in block['orders']])
        self.infinity = np.array(infinity)/ramification


@lru_cache(maxsize=4)
def original_components(n, modes, initial_modes=0):
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel
    selected = geometry_data(n, modes, initial_modes)['selected']
    matrix = pascal_kernel(range(n), 2*n, 2)
    gaussian = {powers: {degree: value for degree, value in enumerate(polynomial) if value}
                for powers, polynomial in direct_kernel(matrix, 2).items()}
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
    return _components(gaussian, left, right, 2), gaussian


def polynomial_seed(n, ramification, modes, kernel, initial_modes=0):
    selected = geometry_data(n, modes, initial_modes)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    components, gaussian = original_components(n, modes, initial_modes)
    clearing = max(0, -min(degree for value in kernel.values() for degree in value))
    tensor = {}
    for index, scalar in kernel.items():
        scalar = {degree+clearing: value for degree, value in scalar.items()}
        for powers, polynomial in components[index].items():
            value = {ramification*degree: coefficient for degree, coefficient in polynomial.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(value, scalar))
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
    repaired = {powers: value for powers, value in repaired.items() if value}
    assert all(evaluate_one(repaired.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
               for powers in set(repaired)|set(gaussian))
    bound = ramification*(max(degree for value in gaussian.values() for degree in value)+2*n*modes)
    bound += max(degree+clearing for value in kernel.values() for degree in value)
    assert max(degree for value in repaired.values() for degree in value) <= bound
    return repaired, clearing, gaussian
