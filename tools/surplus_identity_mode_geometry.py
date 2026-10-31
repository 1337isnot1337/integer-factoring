#!/usr/bin/env python3
"""Exact arbitrary-selection identity geometry above a prime power.

An old reflected Pascal basis leaves a square core on at most r+2d
physical identity directions, with identity operators on the added
actual directions. Terminal selection needs only r+d directions.
There are no formal or dummy modes.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb

from coupled_identity_mode_filters import add, multiply, shift
from coupled_identity_weyl_profiles import adjugate, determinant
from prime_state_filter_profiles import validate


@lru_cache(maxsize=4096)
def form(prime, i, j):
    """Complete f_j(i;1−tau), with no valuation or coefficient cutoff."""
    result = {}
    for h in range(min(i, j)+1):
        scalar = comb(i, h)*comb(j, h) % prime
        if scalar:
            for degree in range(h+1):
                value = scalar*comb(h, degree)*(-1)**degree % prime
                result[degree] = (result.get(degree, 0)+value) % prime
    return {degree: value for degree, value in result.items() if value}


@lru_cache(maxsize=64)
def selected_geometry_data(prime, parent, surplus, selected):
    validate(prime, parent, 1, (0,)*prime)
    power = parent
    while power > 1 and power % prime == 0:
        power //= prime
    if power != 1 or parent <= 1:
        raise ValueError('The parent interval must be a positive power of the characteristic.')
    if not isinstance(surplus, int) or not 1 <= surplus <= parent:
        raise ValueError('Require an integer surplus 1≤d≤N.')
    n = parent+surplus
    selected = tuple(selected)
    if (not selected or len(set(selected)) != len(selected)
            or any(not isinstance(j, int) or not 0 <= j < n for j in selected)):
        raise ValueError('Select distinct actual identity indices in 0,...,n−1.')
    modes = len(selected)
    completed = set(selected)|set(range(surplus))|set(range(parent, n))
    unchanged = tuple(sorted(completed-set(selected)))
    augmented = selected+unchanged
    old_columns = tuple(sorted(j for j in augmented if j < parent))
    retained = tuple(sorted(parent-1-j for j in old_columns))+tuple(range(parent, n))
    assert len(augmented) == len(retained) <= modes+2*surplus
    old_unit = {degree: 1 for degree in range(parent)}
    new_factor = {degree: value for degree, value in ((0, 2 % prime), (parent, prime-1)) if value}
    fx, fv = [], []
    for i in retained:
        xr, vr = [], []
        for j in augmented:
            if i < parent:
                old_column = j if j < parent else j-parent
                present = i == parent-1-old_column
                xr.append({0: 1} if present else {})
                vr.append(dict(old_unit) if present else {})
            else:
                value = form(prime, i-parent, j if j < parent else j-parent)
                if j >= parent:
                    value = multiply(value, new_factor, prime)
                xr.append(value)
                vr.append(value)
        fx.append(tuple(xr))
        fv.append(tuple(vr))
    fx, fv = tuple(fx), tuple(fv)
    dx, dv = determinant(fx, prime), determinant(fv, prime)
    assert dx.get(0, 0) and dv.get(0, 0)
    return {'prime': prime, 'parent': parent, 'n': n, 'surplus': surplus,
            'modes': modes, 'augmented_modes': len(augmented), 'old_augmented': old_columns,
            'selected': selected, 'unchanged_identity_indices': unchanged,
            'augmented': augmented, 'retained': retained,
            'augmented_x': fx, 'augmented_v': fv,
            'determinant_x': dx, 'determinant_v': dv, 'unit': {0: 1}}


def geometry_data(prime, parent, surplus, modes):
    if (not isinstance(modes, int) or not isinstance(surplus, int)
            or not 1 <= surplus <= modes <= parent):
        raise ValueError('The terminal specialization requires 1≤d≤r≤N.')
    return selected_geometry_data(prime, parent, surplus, tuple(range(parent+surplus-modes, parent+surplus)))


@lru_cache(maxsize=64)
def selected_weyl_data(prime, parent, surplus, selected):
    data = selected_geometry_data(prime, parent, surplus, selected)
    fx, fv, size = data['augmented_x'], data['augmented_v'], data['augmented_modes']
    delta = multiply({0: 1, 1: prime-1},
                     multiply(data['determinant_x'], data['determinant_v'], prime), prime)
    ax, av = adjugate(fx, prime), adjugate(fv, prime)
    numerator = []
    for i in range(size):
        row = []
        for j in range(size):
            value = dict(delta) if i == j else {}
            for k, t in enumerate(data['retained']):
                term = multiply(ax[i][k], av[j][k], prime)
                value = add(value, shift(term, data['augmented'][i]+t+1, prime), prime)
            assert value.get(0, 0) == (delta[0] if i == j else 0)
            if value:
                assert min(value) >= max(0, data['augmented'][i]-data['augmented'][j])
            row.append(value)
        numerator.append(tuple(row))
    return {'denominator_unit': delta, 'inverse_parameter_numerator': tuple(numerator)}


def weyl_data(prime, parent, surplus, modes):
    data = geometry_data(prime, parent, surplus, modes)
    return selected_weyl_data(prime, parent, surplus, data['selected'])
