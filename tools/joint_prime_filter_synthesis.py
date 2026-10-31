#!/usr/bin/env python3
"""Exact joint p-state residue synthesis over a fixed finite extension.

The transform has prime-field coefficients, so additive coordinates in
F_(p^(r+1)) suffice for all coefficient equations and profile replays.
The nonzero physical-Y endpoint repair is valid over that fixed field.
"""
from functools import lru_cache
from itertools import product
from math import comb

from polynomial_state_filter_synthesis import _rref, residue_patterns
from prime_power_identity_filters import finite_difference_terms


def extension_add(first, second, prime):
    result, place = 0, 1
    while first or second:
        result += ((first % prime + second % prime) % prime) * place
        first //= prime
        second //= prime
        place *= prime
    return result


def extension_scale(value, scalar, prime):
    scalar %= prime
    result, place = 0, 1
    while value:
        result += (value % prime * scalar % prime) * place
        value //= prime
        place *= prime
    return result


@lru_cache(maxsize=8)
def difference_matrix(prime, modes):
    states = tuple(product(range(prime), repeat=modes))
    lookup = {state: index for index, state in enumerate(states)}
    result = [[{} for _ in states] for _ in states]
    for row, (k, terms) in enumerate(finite_difference_terms(prime, modes).items()):
        for state, power, scalar in terms:
            result[row][lookup[state]] = {degree: scalar * comb(power, degree) * (-1) ** degree % prime
                                          for degree in range(power + 1)
                                          if scalar * comb(power, degree) % prime}
    return tuple(tuple(row) for row in result)


@lru_cache(maxsize=32768)
def local_endpoint_image(prime, modes, degrees, orders):
    size = prime ** modes
    variables = tuple((state, degree) for state, bound in enumerate(degrees)
                      for degree in range(max(0, bound + 1)))
    matrix = difference_matrix(prime, modes)
    equations = []
    for row, order in enumerate(orders):
        maximum = max((degree + power for state, degree in variables for power in matrix[row][state]), default=-1)
        for at in range(min(order, maximum + 1)):
            equations.append([matrix[row][state].get(at - degree, 0) for state, degree in variables])
    reduced, pivots = _rref(equations, len(variables), prime)
    pivot_set = set(pivots)
    image = {}
    for free in range(len(variables)):
        if free in pivot_set:
            continue
        lift = [0] * len(variables)
        lift[free] = 1
        for row, pivot in zip(reduced, pivots):
            lift[pivot] = -row[free] % prime
        endpoint = [0] * size
        for coefficient, (state, degree) in zip(lift, variables):
            endpoint[state] = (endpoint[state] + coefficient) % prime
        while any(endpoint):
            pivot = max(index for index, value in enumerate(endpoint) if value)
            coefficient = endpoint[pivot]
            if pivot in image:
                old, old_lift = image[pivot]
                endpoint = [(a - coefficient * b) % prime for a, b in zip(endpoint, old)]
                lift = [(a - coefficient * b) % prime for a, b in zip(lift, old_lift)]
            else:
                inverse = pow(coefficient, -1, prime)
                image[pivot] = (tuple(value * inverse % prime for value in endpoint),
                                tuple(value * inverse % prime for value in lift))
                break
    return variables, tuple(image.values())


def endpoint_basis(prime, modes, m, degrees, orders, construct=True):
    basis = {}
    for residue, local_degrees, local_orders in residue_patterns(m, degrees, orders):
        variables, image = local_endpoint_image(prime, modes, local_degrees, local_orders)
        for vector, coefficients in image:
            endpoint = list(vector)
            lift = {(state, residue + m * degree): value
                    for (state, degree), value in zip(variables, coefficients) if value} if construct else {}
            while any(endpoint):
                pivot = max(index for index, value in enumerate(endpoint) if value)
                coefficient = endpoint[pivot]
                if pivot in basis:
                    old, old_lift = basis[pivot]
                    endpoint = [(a - coefficient * b) % prime for a, b in zip(endpoint, old)]
                    if construct:
                        for key, value in old_lift.items():
                            lift[key] = (lift.get(key, 0) - coefficient * value) % prime
                            if not lift[key]:
                                del lift[key]
                else:
                    inverse = pow(coefficient, -1, prime)
                    basis[pivot] = (tuple(value * inverse % prime for value in endpoint),
                                    {key: value * inverse % prime for key, value in lift.items()})
                    break
    return tuple(basis.values())


def transform(prime, modes, m, alpha):
    states = tuple(product(range(prime), repeat=modes))
    result = {}
    for state, row in zip(states, difference_matrix(prime, modes)):
        polynomial = {}
        for source, coefficient in zip(states, row):
            for degree, value in alpha[source].items():
                for power, scalar in coefficient.items():
                    at = degree + m * power
                    polynomial[at] = extension_add(polynomial.get(at, 0),
                                                    extension_scale(value, scalar, prime), prime)
        result[state] = {degree: value for degree, value in polynomial.items() if value}
    return result


def synthesize(prime, modes, m, degrees, orders, construct=True):
    size, field_size = prime ** modes, prime ** (modes + 1)
    degrees, orders = tuple(degrees), tuple(orders)
    if (m < 1 or len(degrees) != size or len(orders) != size
            or any(not isinstance(value, int) or value < 0 for value in degrees + orders)):
        raise ValueError('Use positive m and one nonnegative degree/order bound per joint state.')
    basis = endpoint_basis(prime, modes, m, degrees, orders, construct)
    if not all(any(vector[j] for vector, lift in basis) for j in range(size)):
        return None
    if not construct:
        return True
    endpoint = [0] * size
    polynomial = {}
    for vector, lift in basis:
        union = tuple(j for j in range(size) if endpoint[j] or vector[j])
        if all(endpoint[j] for j in union):
            continue
        coefficient = next(value for value in range(field_size) if all(
            extension_add(endpoint[j], extension_scale(value, vector[j], prime), prime) for j in union))
        endpoint = [extension_add(a, extension_scale(coefficient, b, prime), prime)
                    for a, b in zip(endpoint, vector)]
        for key, value in lift.items():
            polynomial[key] = extension_add(polynomial.get(key, 0), extension_scale(coefficient, value, prime), prime)
            if not polynomial[key]:
                del polynomial[key]
        if all(endpoint):
            break
    assert all(endpoint)
    states = tuple(product(range(prime), repeat=modes))
    alpha = {state: {degree: value for (index, degree), value in polynomial.items() if index == i}
             for i, state in enumerate(states)}
    beta = transform(prime, modes, m, alpha)
    actual_degrees = tuple(max(alpha[state]) for state in states)
    actual_orders = tuple(min(beta[state]) for state in states)
    assert all(a <= b for a, b in zip(actual_degrees, degrees))
    assert all(a >= b for a, b in zip(actual_orders, orders))
    return {'prime': prime, 'modes': modes, 'm': m, 'coefficient_field_order': field_size,
            'degree_bounds': degrees, 'order_bounds': orders,
            'endpoint_image_dimension': len(basis), 'endpoint_values': tuple(endpoint),
            'alpha': alpha, 'beta': beta, 'actual_degrees': actual_degrees, 'actual_orders': actual_orders,
            'residue_pattern_count': len(residue_patterns(m, degrees, orders))}
