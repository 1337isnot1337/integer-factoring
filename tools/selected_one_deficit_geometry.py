#!/usr/bin/env python3
"""Certified arbitrary identity selections at n=p^k-1."""
from __future__ import annotations

from functools import lru_cache
from math import comb

from coupled_identity_mode_filters import add, multiply, shift
from coupled_identity_weyl_profiles import determinant, adjugate
from prime_state_filter_profiles import validate


@lru_cache(maxsize=128)
def geometry_data(prime, parent, selected):
    validate(prime, parent, 1, (0,)*prime)
    power = parent
    while power > 1 and power % prime == 0:
        power //= prime
    if power != 1 or parent <= 1:
        raise ValueError('The parent interval must be a power of the characteristic.')
    selected = tuple(selected)
    n, modes = parent-1, len(selected)
    if not selected or selected != tuple(sorted(set(selected))) or min(selected) < 0 or max(selected) >= n:
        raise ValueError('Use distinct increasing actual indices in0..n−1.')
    if modes+1 > n:
        raise ValueError('The one-dummy reduction requires r+1≤n.')
    augmented = selected+(n,)
    prefix = next(j for j in range(parent) if j not in augmented)
    retained = tuple(sorted(n-j for j in augmented if j > prefix))+tuple(range(n-prefix, n))
    assert len(retained) == modes+1 and retained == tuple(sorted(set(retained)))
    discarded = tuple(i for i in range(n) if i not in retained)

    def form(i, j):
        reverse = parent-1-j
        value = {}
        if j <= reverse:
            for h in range(min(i, j)+1):
                for degree in range(h+1):
                    coefficient = comb(i, h)*comb(j, h)*comb(h, degree)*(-1)**degree % prime
                    value[degree] = (value.get(degree, 0)+coefficient) % prime
        else:
            for h in range(min(i, reverse)+1):
                for degree in range(h+1):
                    coefficient = comb(i, h)*comb(reverse, h)*comb(h, degree)*(-1)**(h+degree) % prime
                    exponent = i-h+degree
                    value[exponent] = (value.get(exponent, 0)+coefficient) % prime
        return {degree: coefficient for degree, coefficient in value.items() if coefficient}

    forms = tuple(tuple(form(i, j) for j in augmented) for i in range(n))
    fx = tuple(forms[i] for i in retained)
    unit = determinant(fx, prime)
    assert unit.get(0, 0)
    inverse_numerator = adjugate(fx, prime)
    interpolation, margins = [], []
    for i in discarded:
        row = []
        for a, t in enumerate(retained):
            value = {}
            for j in range(modes+1):
                value = add(value, multiply(forms[i][j], inverse_numerator[j][a], prime), prime)
            assert not value or min(value) >= max(0, i-t)
            margins.append(None if not value else min(value)-max(0, i-t))
            row.append(value)
        interpolation.append(tuple(row))
    fv = []
    for a, t in enumerate(retained):
        row = []
        for j in range(modes+1):
            value = multiply(unit, fx[a][j], prime)
            for k, i in enumerate(discarded):
                value = add(value, shift(multiply(interpolation[k][a], forms[i][j], prime), t-i, prime), prime)
            assert not value or min(value) >= 0
            row.append(value)
        fv.append(tuple(row))
    fv = tuple(fv)
    assert determinant(fv, prime).get(0, 0)
    return {'prime': prime, 'parent': parent, 'n': n, 'deficit': 1, 'modes': modes,
            'augmented_modes': modes+1, 'selected': selected, 'dummy': (n,), 'augmented': augmented,
            'initial_prefix_length': prefix, 'retained': retained, 'discarded': discarded,
            'forms': forms, 'unit': unit, 'interpolation_numerators': tuple(interpolation),
            'weighted_order_margins': tuple(margins), 'augmented_x': fx, 'augmented_v': fv}


def compact_parameter(data):
    prime, parent, selected = data['prime'], data['parent'], data['augmented']
    c, d = {0: 1, 1: prime-1}, {0: 1, parent: prime-1}
    e = add(d, multiply(c, {j: 1 for j in selected}, prime), prime, -1)
    h = min(e)
    assert h == data['initial_prefix_length']
    e0 = shift(e, -h, prime)
    delta = multiply(d, e0, prime)
    numerator = tuple(tuple(add(e0 if i == j else {}, shift(c, parent-h+index, prime), prime)
                            for j in range(len(selected))) for i, index in enumerate(selected))
    return {'denominator_unit': delta, 'inverse_parameter_numerator': numerator,
            'missing_order': h, 'missing_polynomial_unit': e0}


def adjugate_parameter(data):
    prime, fx, fv = data['prime'], data['augmented_x'], data['augmented_v']
    c = {0: 1, 1: prime-1}
    delta = multiply(c, multiply(determinant(fx, prime), determinant(fv, prime), prime), prime)
    ax, av = adjugate(fx, prime), adjugate(fv, prime)
    numerator = []
    for i, selected in enumerate(data['augmented']):
        row = []
        for j in range(data['augmented_modes']):
            value = dict(delta) if i == j else {}
            for k, retained in enumerate(data['retained']):
                term = multiply(multiply(ax[i][k], av[j][k], prime), data['unit'], prime)
                value = add(value, shift(term, selected+retained+1, prime), prime)
            row.append(value)
        numerator.append(tuple(row))
    return {'denominator_unit': delta, 'inverse_parameter_numerator': tuple(numerator)}


def contraction_eigenbasis(data):
    prime, indices = data['prime'], data['augmented']
    minimum = min(indices)
    size = len(indices)
    basis = []
    for i, index in enumerate(indices):
        row = [{index-minimum: 1}]
        row += [({0: 1} if i == j else {0: prime-1} if i == 0 else {}) for j in range(1, size)]
        basis.append(tuple(row))
    basis = tuple(basis)
    unit = determinant(basis, prime)
    assert unit.get(0, 0)
    inverse_numerator = adjugate(basis, prime)
    for matrix in (basis, inverse_numerator):
        for i in range(size):
            for j in range(size):
                assert not matrix[i][j] or min(matrix[i][j]) >= max(0, indices[i]-indices[j])
    return {'basis': basis, 'determinant_unit': unit,
            'inverse_numerator': inverse_numerator,
            'endpoint_determinant': sum(unit.values()) % prime}
