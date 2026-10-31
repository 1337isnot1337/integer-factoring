#!/usr/bin/env python3
"""Exact four-mode resonant blocks, using two-by-two matrices and scalars."""
from __future__ import annotations

from collections import defaultdict
from itertools import product

import numpy as np
from scipy.special import logsumexp

from coupled_identity_operator_filters import identity_kernel
from f4_repaired_identity_filters import add, multiply, evaluate_one
from prime_state_filter_pressure import gaussian_pressure


STATES = tuple(product(range(2), repeat=4))
FIRST, SECOND = (1, 0, 0, 1), (0, 1, 1, 0)


def validate_selection(n, selected):
    selected = tuple(selected)
    if n <= 0 or n & (n-1):
        raise ValueError('Use a power-of-two interval length.')
    if len(selected) != 4 or selected != tuple(sorted(set(selected))) or min(selected) < 0 or max(selected) >= n:
        raise ValueError('Select four distinct increasing interval indices.')
    if selected[0]+selected[3] != selected[1]+selected[2]:
        raise ValueError('The complementary two-particle states must have equal energy.')
    return selected


def determinant(matrix):
    return add(multiply(matrix[0][0], matrix[1][1]), multiply(matrix[0][1], matrix[1][0]))


def matrix_orders(matrix, infinity=False):
    det = determinant(matrix)
    if not det:
        raise ArithmeticError('A nonzero-endpoint resonance matrix became singular.')
    if infinity:
        first = -2*max(degree for row in matrix for value in row for degree in value)
        total = -2*max(det)
    else:
        first = min(degree for row in matrix for value in row for degree in value)
        total = min(det)
    result = first, total-first
    assert result[0] <= result[1]
    return result


def signature(n, ramification, matrix):
    """Return unweighted zero scalar/pair orders and the infinity pair orders."""
    if ramification <= 0 or ramification % 2 == 0:
        raise ValueError('Use a positive odd ramification.')
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError('The resonant occupation block must be two by two.')
    if any(not isinstance(degree, int) or not isinstance(value, int) or not 0 <= value < 4
           for row in matrix for entry in row for degree, value in entry.items()):
        raise ValueError('Use complete F4 Laurent polynomial entries.')
    matrix = tuple(tuple({degree: value for degree, value in entry.items() if value}
                         for entry in row) for row in matrix)
    if not evaluate_one(determinant(matrix)):
        raise ValueError('The resonant endpoint block must be invertible.')
    m = n*ramification
    c = {0: 1, m: 1}
    delta_a, delta_b = add(matrix[0][0], {0: 1}), add(matrix[1][1], {0: 1})
    vacuum = add({4*m: 1}, multiply({0: 1, 2*m: 1}, add(delta_a, delta_b)))
    single_a = add({3*m: 1}, multiply(c, delta_a))
    single_b = add({3*m: 1}, multiply(c, delta_b))
    pair = tuple(tuple(add(matrix[i][j], {0: 1, 2*m: 1}) if i == j else matrix[i][j]
                       for j in range(2)) for i in range(2))
    assert vacuum and single_a and single_b
    return (min(vacuum), min(single_a), min(single_b),
            *matrix_orders(pair), *matrix_orders(matrix, True))


def kernel(matrix):
    result = identity_kernel(2, 4)
    for i, row in enumerate((FIRST, SECOND)):
        for j, column in enumerate((FIRST, SECOND)):
            result[row, column] = dict(matrix[i][j])
    return {index: value for index, value in result.items() if value}


def reduced_orders(n, ramification, selected, values):
    selected = validate_selection(n, selected)
    vacuum, single_a, single_b, pair_a, pair_b, infinity_a, infinity_b = values
    result = []
    for state in STATES:
        occupation = sum(state)
        weight = ramification*sum(j*digit for j, digit in zip(selected, state))
        if occupation == 0:
            extra = vacuum
        elif occupation == 1:
            extra = single_a if state[0] or state[3] else single_b
        elif state == FIRST:
            extra = pair_a
        elif state == SECOND:
            extra = pair_b
        else:
            extra = (4-occupation)*n*ramification
        result.append(weight+extra)
    return tuple(result)


def profiles(n, ramification, selected, values):
    selected = validate_selection(n, selected)
    reduced = reduced_orders(n, ramification, selected, values)
    zero = defaultdict(list)
    for state, order in zip(STATES, reduced):
        for auxiliary in STATES:
            shift = -2*ramification*sum(j*digit for j, digit in zip(selected, auxiliary))
            zero[sum(state)-sum(auxiliary)].append(2*order+shift)
    infinity = []
    for state in STATES:
        order = ramification*sum((2*j+1-2*n)*digit for j, digit in zip(selected, state))
        if state == FIRST:
            order += values[5]
        elif state == SECOND:
            order += values[6]
        infinity.append(order)
    return (tuple({'charge': charge, 'orders': tuple(sorted(orders))} for charge, orders in sorted(zero.items())),
            tuple(sorted(infinity)))


def pressure(n, ramification, selected, values, scales, charges, gaussian=None):
    """Full pressure, cancelling the common auxiliary partition exactly."""
    selected = validate_selection(n, selected)
    scales, charges = np.broadcast_arrays(np.asarray(scales), np.asarray(charges))
    if scales.ndim != 1:
        raise ValueError('This vectorized helper accepts one-dimensional scale and charge arrays.')
    occupations = np.array([sum(state) for state in STATES])
    weights = np.array([sum(j*digit for j, digit in zip(selected, state)) for state in STATES])
    zero = 2*np.array(reduced_orders(n, ramification, selected, values))/ramification
    baseline_zero = 2*(weights+(4-occupations)*n)
    infinity = 2*weights+(1-2*n)*occupations
    changed_infinity = infinity.astype(float)
    changed_infinity[STATES.index(FIRST)] += values[5]/ramification
    changed_infinity[STATES.index(SECOND)] += values[6]/ramification
    shift = occupations[:, None]*(1-charges)
    change = logsumexp(-(zero[:, None]+shift)*scales, axis=0)
    change -= logsumexp(-(baseline_zero[:, None]+shift)*scales, axis=0)
    change += logsumexp(-changed_infinity[:, None]*scales, axis=0)
    change -= logsumexp(-infinity[:, None]*scales, axis=0)
    if gaussian is None:
        gaussian = np.array([gaussian_pressure(2, n, float(J), float(t)) for t, J in zip(scales, charges)])
    return change+gaussian
