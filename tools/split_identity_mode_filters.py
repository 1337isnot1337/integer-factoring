#!/usr/bin/env python3
"""Certified mixed-end identity cores from exact polynomial interpolation.

Select a initial and b terminal identity directions. Retain b lowest
and a highest auxiliary modes. The finite constructor verifies both
regularity conditions of the existing weighted concentration theorem.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb, factorial

from coupled_identity_mode_filters import add, multiply, shift, algebra_multiply
from prime_state_filter_profiles import validate


@lru_cache(maxsize=64)
def geometry_data(prime, n, initial_modes=1, terminal_modes=2):
    validate(prime, n, 1, (0,) * prime)
    modes = initial_modes + terminal_modes
    if min(initial_modes, terminal_modes) < 0 or not 1 <= modes <= n:
        raise ValueError('Select nonnegative end counts, with total between one and n.')
    selected = tuple(range(initial_modes)) + tuple(range(n-terminal_modes, n))
    retained = tuple(range(terminal_modes)) + tuple(range(n-initial_modes, n))
    discarded = tuple(i for i in range(n) if i not in retained)

    def form(i, j):
        value = {}
        for h in range(min(i, j)+1):
            scalar = comb(i, h) * comb(j, h)
            for degree in range(h+1):
                coefficient = scalar * comb(h, degree) * (-1)**degree % prime
                value[degree] = (value.get(degree, 0)+coefficient) % prime
        return {degree: coefficient for degree, coefficient in value.items() if coefficient}

    forms = tuple(tuple(form(i, j) for j in selected) for i in range(n))
    retained_matrix = tuple(forms[i] for i in retained)
    # Import lazily: the Weyl layer calls this constructor through its
    # geometry interface, while these two polynomial operations are pure.
    from coupled_identity_weyl_profiles import determinant, adjugate
    unit = determinant(retained_matrix, prime)
    if not unit.get(0, 0):
        raise ArithmeticError('The proposed retained minor is not a zero unit.')
    adjoint = adjugate(retained_matrix, prime)
    numerators = []
    margins = []
    for i in discarded:
        row = []
        for a, t in enumerate(retained):
            value = {}
            for j in range(modes):
                value = add(value, multiply(forms[i][j], adjoint[j][a], prime), prime)
            required = max(0, i-t)
            if value and min(value) < required:
                raise ArithmeticError('The proposed interpolation fails weighted regularity.')
            margins.append(None if not value else min(value)-required)
            row.append(value)
        numerators.append(tuple(row))
    active_v = []
    for a, t in enumerate(retained):
        row = []
        for j in range(modes):
            value = multiply(unit, retained_matrix[a][j], prime)
            for k, i in enumerate(discarded):
                correction = multiply(numerators[k][a], forms[i][j], prime)
                value = add(value, shift(correction, t-i, prime), prime)
            assert not value or min(value) >= 0
            row.append(value)
        active_v.append(tuple(row))
    active_v = tuple(active_v)
    assert determinant(active_v, prime).get(0, 0)
    return {'prime': prime, 'n': n, 'modes': modes,
            'geometry': f'split:{initial_modes}:{terminal_modes}',
            'selected': selected, 'retained': retained, 'discarded': discarded,
            'all_forms': forms, 'interpolation_numerators': tuple(numerators),
            'interpolation_denominator_unit': unit, 'weighted_order_margins': tuple(margins),
            'active_x': retained_matrix, 'active_v': active_v, 'gaussian_v_unit': unit}


def active_gaussian(prime, n, ramification, modes, geometry):
    counts = tuple(map(int, geometry.split(':')[1:]))
    if len(counts) != 2 or sum(counts) != modes:
        raise ValueError('The split counts must sum to the number of modes.')
    data = geometry_data(prime, n, *counts)
    zero = (0,) * (4*modes)
    gaussian = {zero: {0: 1}}
    unit = {2*ramification*degree: value for degree, value in data['gaussian_v_unit'].items()}
    for j in range(modes):
        a = 2*ramification*data['selected'][j]
        b = 2*ramification*(data['retained'][j]+1)
        for left, right, order, scalar in ((modes+j, 2*modes+j, a, {0: 1}),
                                           (j, 3*modes+j, b, unit)):
            factor, power = {}, {0: 1}
            for digit in range(prime):
                state = list(zero)
                state[left] = state[right] = digit
                factor[tuple(state)] = shift(power, order*digit, prime, pow(factorial(digit), -1, prime))
                power = multiply(power, scalar, prime)
            gaussian = algebra_multiply(gaussian, factor, prime)
    return gaussian
