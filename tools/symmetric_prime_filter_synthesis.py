#!/usr/bin/env python3
"""Exact prime-state symmetric filters in r*(p-1)+1 occupation weights.

The Vandermonde grouping retains all digit binomial coefficients. Here
f_w=alpha_(p-1-k) for |k|=w, and beta_w=sum_h (-1)^h binom(w,h)
(1-x^m)^(w-h) f_h. Endpoint beta_w(1)=(-1)^w f_w(1).
"""

from functools import lru_cache
from math import comb, isqrt

from joint_prime_filter_synthesis import extension_add, extension_scale
from polynomial_state_filter_synthesis import _rref, residue_patterns


@lru_cache(maxsize=16)
def difference_matrix(prime, modes):
    last = modes * (prime - 1)
    return tuple(tuple({degree: (-1) ** (source + degree) * comb(weight, source)
                                 * comb(weight - source, degree) % prime
                        for degree in range(weight - source + 1)
                        if comb(weight, source) * comb(weight - source, degree) % prime}
                       if source <= weight else {}
                       for source in range(last + 1))
                 for weight in range(last + 1))


@lru_cache(maxsize=8192)
def local_endpoint_image(prime, modes, degrees, orders):
    size = modes * (prime - 1) + 1
    variables = tuple((weight, degree) for weight, bound in enumerate(degrees)
                      for degree in range(max(0, bound + 1)))
    matrix = difference_matrix(prime, modes)
    equations = []
    for weight, order in enumerate(orders):
        maximum = max((degree + power for source, degree in variables
                       for power in matrix[weight][source]), default=-1)
        for at in range(min(order, maximum + 1)):
            equations.append([matrix[weight][source].get(at - degree, 0)
                              for source, degree in variables])
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
        for coefficient, (weight, degree) in zip(lift, variables):
            endpoint[weight] = (endpoint[weight] + coefficient) % prime
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


def transform(prime, modes, m, values):
    result = []
    for row in difference_matrix(prime, modes):
        value = {}
        for polynomial, multiplier in zip(values, row):
            for degree, coefficient in polynomial.items():
                for power, scalar in multiplier.items():
                    at = degree + m * power
                    value[at] = extension_add(value.get(at, 0),
                                              extension_scale(coefficient, scalar, prime), prime)
        result.append({degree: coefficient for degree, coefficient in value.items() if coefficient})
    return tuple(result)


def synthesize(prime, modes, m, degrees, orders, construct=True):
    size = modes * (prime - 1) + 1
    degrees, orders = tuple(degrees), tuple(orders)
    if (prime < 2 or any(prime % d == 0 for d in range(2, isqrt(prime) + 1))
            or modes < 1 or m < 1 or len(degrees) != size or len(orders) != size
            or any(not isinstance(value, int) or value < 0 for value in degrees + orders)):
        raise ValueError('Use a prime, positive modes and m, and one nonnegative integer bound per weight.')
    patterns = residue_patterns(m, degrees, orders)
    basis = {}
    for residue, local_degrees, local_orders in patterns:
        variables, image = local_endpoint_image(prime, modes, local_degrees, local_orders)
        for vector, coefficients in image:
            endpoint = list(vector)
            lift = {(weight, residue + m * degree): value
                    for (weight, degree), value in zip(variables, coefficients) if value} if construct else {}
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
    if not all(any(vector[index] for vector, lift in basis.values()) for index in range(size)):
        return None
    if not construct:
        return True
    field_size = prime
    while field_size <= size:
        field_size *= prime
    endpoint = [0] * size
    polynomial = {}
    for vector, lift in basis.values():
        union = tuple(index for index in range(size) if endpoint[index] or vector[index])
        if all(endpoint[index] for index in union):
            continue
        scalar = next(value for value in range(field_size) if all(extension_add(
            endpoint[index], extension_scale(value, vector[index], prime), prime) for index in union))
        endpoint = [extension_add(a, extension_scale(scalar, b, prime), prime)
                    for a, b in zip(endpoint, vector)]
        for key, value in lift.items():
            polynomial[key] = extension_add(polynomial.get(key, 0),
                                            extension_scale(scalar, value, prime), prime)
            if not polynomial[key]:
                del polynomial[key]
        if all(endpoint):
            break
    values = tuple({degree: value for (weight, degree), value in polynomial.items() if weight == index}
                   for index in range(size))
    beta = transform(prime, modes, m, values)
    actual_degrees, actual_orders = tuple(map(max, values)), tuple(map(min, beta))
    assert all(endpoint)
    assert all(actual <= bound for actual, bound in zip(actual_degrees, degrees))
    assert all(actual >= bound for actual, bound in zip(actual_orders, orders))
    return {
        'prime': prime, 'modes': modes, 'm': m, 'coefficient_field_order': field_size,
        'degree_bounds': degrees, 'order_bounds': orders,
        'actual_degrees': actual_degrees, 'actual_orders': actual_orders,
        'endpoint_values': tuple(endpoint), 'endpoint_image_dimension': len(basis),
        'residue_pattern_count': len(patterns), 'f': values, 'beta': beta,
    }
