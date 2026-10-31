#!/usr/bin/env python3
"""Exact mixed-end identity geometry a fixed deficit below a prime power.

This module supplies the polynomial concentration and square dummy-Weyl
parameter. It does not implement a new pressure search or claim a gain.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from math import comb

from coupled_identity_mode_filters import add, multiply, shift
from coupled_identity_weyl_profiles import adjugate, determinant
from prime_state_filter_profiles import validate


@lru_cache(maxsize=64)
def geometry_data(prime, parent, deficit, initial_modes, terminal_modes):
    validate(prime, parent, 1, (0,)*prime)
    power = parent
    while power > 1 and power % prime == 0:
        power //= prime
    if power != 1 or parent <= 1:
        raise ValueError('The parent interval must be a positive power of the characteristic.')
    if any(not isinstance(value, int) or value < 0
           for value in (deficit, initial_modes, terminal_modes)):
        raise ValueError('Use nonnegative integer deficit and end counts.')
    n, modes = parent-deficit, initial_modes+terminal_modes
    size = modes+deficit
    if not 1 <= modes or size > n:
        raise ValueError('The actual interval must contain all r+d retained auxiliary modes.')
    selected = tuple(range(initial_modes))+tuple(range(n-terminal_modes, n))
    dummy = tuple(range(n, parent))
    augmented = selected+dummy
    retained = tuple(range(terminal_modes+deficit))+tuple(range(n-initial_modes, n))
    discarded = tuple(i for i in range(n) if i not in retained)

    def form(i, j):
        result = {}
        if j < initial_modes:
            for h in range(min(i, j)+1):
                for degree in range(h+1):
                    value = comb(i, h)*comb(j, h)*comb(h, degree)*(-1)**degree % prime
                    result[degree] = (result.get(degree, 0)+value) % prime
        else:
            reverse = parent-1-j
            for h in range(min(i, reverse)+1):
                for degree in range(h+1):
                    value = comb(i, h)*comb(reverse, h)*comb(h, degree)*(-1)**(h+degree) % prime
                    exponent = i-h+degree
                    result[exponent] = (result.get(exponent, 0)+value) % prime
        return {degree: value for degree, value in result.items() if value}

    all_forms = tuple(tuple(form(i, j) for j in augmented) for i in range(parent))
    forms, tail_forms = all_forms[:n], all_forms[n:]
    fx = tuple(forms[i] for i in retained)
    unit = determinant(fx, prime)
    assert unit.get(0, 0)
    inverse_numerator = adjugate(fx, prime)
    interpolation, margins = [], []
    for i in discarded:
        row = []
        for k, t in enumerate(retained):
            value = {}
            for j in range(size):
                value = add(value, multiply(forms[i][j], inverse_numerator[j][k], prime), prime)
            assert not value or min(value) >= max(0, i-t)
            margins.append(None if not value else min(value)-max(0, i-t))
            row.append(value)
        interpolation.append(tuple(row))
    fv = []
    for a, t in enumerate(retained):
        row = []
        for j in range(size):
            value = multiply(unit, fx[a][j], prime)
            for k, i in enumerate(discarded):
                term = multiply(interpolation[k][a], forms[i][j], prime)
                value = add(value, shift(term, t-i, prime), prime)
            assert not value or min(value) >= 0
            row.append(value)
        fv.append(tuple(row))
    fv = tuple(fv)
    assert determinant(fv, prime).get(0, 0)
    return {'prime': prime, 'parent': parent, 'n': n, 'deficit': deficit,
            'modes': modes, 'augmented_modes': size, 'initial_modes': initial_modes,
            'terminal_modes': terminal_modes, 'selected': selected, 'dummy': dummy,
            'augmented': augmented, 'retained': retained, 'discarded': discarded,
            'forms': forms, 'tail_forms': tail_forms,
            'interpolation_numerators': tuple(interpolation),
            'weighted_order_margins': tuple(margins), 'unit': unit,
            'augmented_x': fx, 'augmented_v': fv}


def weyl_data(prime, parent, deficit, initial_modes, terminal_modes):
    data = geometry_data(prime, parent, deficit, initial_modes, terminal_modes)
    fx, fv, size = data['augmented_x'], data['augmented_v'], data['augmented_modes']
    c = {0: 1, 1: prime-1}
    delta = multiply(c, multiply(determinant(fx, prime), determinant(fv, prime), prime), prime)
    ax, av = adjugate(fx, prime), adjugate(fv, prime)
    numerator = []
    for i in range(size):
        row = []
        for j in range(size):
            value = dict(delta) if i == j else {}
            for k, t in enumerate(data['retained']):
                term = multiply(multiply(ax[i][k], av[j][k], prime), data['unit'], prime)
                value = add(value, shift(term, data['augmented'][i]+t+1, prime), prime)
            assert value.get(0, 0) == (delta[0] if i == j else 0)
            if value:
                assert min(value) >= max(0, data['augmented'][i]-data['augmented'][j])
            row.append(value)
        numerator.append(tuple(row))
    return {'denominator_unit': delta, 'inverse_parameter_numerator': tuple(numerator)}


def remove_dummy_factor(blocks, prime, order, minimum_charge, maximum_charge):
    """Divide one known Gaussian digit factor in integer multiplicities."""
    result = {}
    for charge in range(minimum_charge, maximum_charge+1):
        current = Counter(blocks.get(charge, ()))
        for digit in range(1, prime):
            for value, count in result.get(charge-digit, Counter()).items():
                current[value+digit*order] -= count
        if any(count < 0 for count in current.values()):
            raise ArithmeticError('A promised Gaussian digit factor is absent.')
        result[charge] = +current
    residual_start = maximum_charge-(prime-1)+1
    if any(result[charge] for charge in range(residual_start, maximum_charge+1)):
        raise ArithmeticError('The dummy Gaussian division has a nonzero final residual.')
    return {charge: tuple(sorted(values.elements())) for charge, values in result.items()
            if charge < residual_start}
