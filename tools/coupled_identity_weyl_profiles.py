#!/usr/bin/env python3
"""Exact Weyl reduction of coupled identity-operator zero profiles.

The original p^(2r) active matrix is replaced by a p^r matrix and
explicit auxiliary digit shifts. All auxiliary rational denominators
are cleared by one proved zero unit. Complete polynomial arrays are
passed to the exact Smith routine; no series truncation is used here.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import factorial

from coupled_identity_mode_filters import add, multiply, shift, algebra_multiply
from coupled_identity_operator_filters import (
    _geometry, _laurent_orders, normalized_kernel, states,
)
from prime_state_filter_profiles import validate


def determinant(matrix, prime):
    size = len(matrix)
    if size == 0:
        return {0: 1}
    if size == 1:
        return dict(matrix[0][0])
    result = {}
    for column, value in enumerate(matrix[0]):
        if not value:
            continue
        minor = [[row[j] for j in range(size) if j != column] for row in matrix[1:]]
        result = add(result, multiply(value, determinant(minor, prime), prime), prime, (-1)**column)
    return result


def adjugate(matrix, prime):
    size = len(matrix)
    return tuple(tuple(shift(determinant(
        [[matrix[a][b] for b in range(size) if b != i]
         for a in range(size) if a != j], prime), 0, prime, (-1)**(i+j))
        for j in range(size)) for i in range(size))


@lru_cache(maxsize=64)
def weyl_data(prime, n, modes=2, geometry='initial'):
    """Return Z=C^-1=N/Delta as complete polynomials in tau=q^P."""
    data = _geometry(prime, n, modes, geometry)
    fx, fv = data['active_x'], data['active_v']
    selected, retained = data['selected'], data['retained']
    unit = data.get('gaussian_v_unit', {0: 1})
    dx, dv = determinant(fx, prime), determinant(fv, prime)
    assert dx.get(0, 0) and dv.get(0, 0)
    delta = multiply({0: 1, 1: -1 % prime}, multiply(dx, dv, prime), prime)
    ax, av = adjugate(fx, prime), adjugate(fv, prime)
    numerator = []
    for i in range(modes):
        row = []
        for j in range(modes):
            value = dict(delta) if i == j else {}
            for k in range(modes):
                term = multiply(multiply(ax[i][k], av[j][k], prime), unit, prime)
                value = add(value, shift(term, selected[i] + retained[k] + 1, prime), prime)
            assert value.get(0, 0) == (delta[0] if i == j else 0)
            if i != j and value:
                assert min(value) >= selected[i] - selected[j]
            row.append(value)
        numerator.append(tuple(row))
    return {'selected': selected, 'retained': retained, 'denominator_unit': delta,
            'inverse_parameter_numerator': tuple(numerator)}


@lru_cache(maxsize=32)
def reduced_components(prime, n, modes=2, geometry='initial'):
    """Sparse entries of Delta^(r(p-1)) Psi_C(E_ab), in tau."""
    data = weyl_data(prime, n, modes, geometry)
    delta, numerator = data['denominator_unit'], data['inverse_parameter_numerator']
    alphabet = states(prime, modes)
    lookup = {state: i for i, state in enumerate(alphabet)}
    maximum = modes * (prime - 1)
    delta_powers = [{0: 1}]
    for exponent in range(maximum):
        delta_powers.append(multiply(delta_powers[-1], delta, prime))
    forms = []
    for j in range(modes):
        form = {}
        for i in range(modes):
            powers = [0] * modes
            powers[i] = 1
            form[tuple(powers)] = shift(numerator[j][i], 0, prime, -1)
        forms.append({powers: value for powers, value in form.items() if value})
    cache = {(0,) * modes: {(0,) * modes: {0: 1}}}

    def monomial(powers):
        if powers not in cache:
            j = next(i for i, digit in enumerate(powers) if digit)
            previous = list(powers)
            previous[j] -= 1
            cache[powers] = algebra_multiply(forms[j], monomial(tuple(previous)), prime)
        return cache[powers]

    result = {}
    for b in alphabet:
        exponents = tuple(prime - 1 - digit for digit in b)
        unit_power = delta_powers[maximum - sum(exponents)]
        polynomial = {powers: multiply(value, unit_power, prime)
                      for powers, value in monomial(exponents).items()}
        for a in alphabet:
            if sum(a) != sum(b):
                continue
            derivatives = tuple(prime - 1 - digit for digit in a)
            scalar = (-1)**modes
            for digit in a:
                scalar *= factorial(digit)
            matrix = {}
            for column, source in enumerate(alphabet):
                for powers, value in polynomial.items():
                    total = tuple(x + y for x, y in zip(source, powers))
                    if any(digit >= prime for digit in total):
                        continue
                    if any(digit < derivative for digit, derivative in zip(total, derivatives)):
                        continue
                    target = tuple(digit - derivative for digit, derivative in zip(total, derivatives))
                    coefficient = scalar
                    for digit, derivative in zip(total, derivatives):
                        coefficient *= factorial(digit) // factorial(digit - derivative)
                    coefficient %= prime
                    if coefficient:
                        entry = (lookup[target], column)
                        matrix[entry] = add(matrix.get(entry, {}), value, prime, coefficient)
            result[a, b] = {entry: value for entry, value in matrix.items() if value}
    return result


def reduced_matrix(prime, n, ramification, modes, kernel, geometry='initial', endpoint='identity'):
    validate(prime, n, ramification, (0,) * prime)
    data = weyl_data(prime, n, modes, geometry)
    kernel = normalized_kernel(prime, modes, kernel, endpoint, data['selected'])
    alphabet = states(prime, modes)
    size = len(alphabet)
    matrix = [[{} for j in range(size)] for i in range(size)]
    components = reduced_components(prime, n, modes, geometry)
    for index, scalar in kernel.items():
        for (row, column), value in components[index].items():
            value = {ramification * degree: coefficient for degree, coefficient in value.items()}
            matrix[row][column] = add(matrix[row][column], multiply(value, scalar, prime), prime)
    for row, state in enumerate(alphabet):
        weight = ramification * sum(j * digit for j, digit in zip(data['selected'], state))
        matrix[row] = [shift(value, weight, prime) for value in matrix[row]]
    return matrix


def zero_profiles(prime, n, ramification, modes, kernel, geometry='initial', endpoint='identity'):
    """Return the same charge-block z-orders as the full active engine."""
    alphabet = states(prime, modes)
    selected = weyl_data(prime, n, modes, geometry)['selected']
    matrix = reduced_matrix(prime, n, ramification, modes, kernel, geometry, endpoint)
    reduced = []
    for occupation in range(modes * (prime - 1) + 1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        if not indices:
            continue
        block = [[matrix[i][j] for j in indices] for i in indices]
        orders = tuple(2 * value for value in _laurent_orders(block, prime))
        reduced.append((occupation, orders))
    result = defaultdict(list)
    for auxiliary in alphabet:
        shift_order = -2 * ramification * sum(j * digit for j, digit in zip(selected, auxiliary))
        for occupation, orders in reduced:
            result[occupation - sum(auxiliary)].extend(order + shift_order for order in orders)
    return tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                 for charge, orders in sorted(result.items()))
