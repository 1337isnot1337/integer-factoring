#!/usr/bin/env python3
"""Exact residue-compressed synthesis for joint binary diagonal filters.

The tensor-product finite-difference matrix depends only on x^m.  At most
2*2^r+1 residue patterns suffice.  A fixed coefficient extension of order
strictly greater than 2^r makes nonzero endpoint coordinate projections
equivalent to endpoint-repair feasibility.  Its additive coordinates are
enough for all equations and exact beta-profile calculations.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from polynomial_state_filter_synthesis import residue_patterns
from prime_power_identity_filters import finite_difference_terms, beta_polynomials
from weyl_lattice_filter_synthesis import endpoint_image


@lru_cache(maxsize=8)
def difference_columns(modes):
    states = tuple(product(range(2), repeat=modes))
    lookup = {state: index for index, state in enumerate(states)}
    columns = [[] for _ in states]
    for row, (k, terms) in enumerate(finite_difference_terms(2, modes).items()):
        for state, power, scalar in terms:
            shifts = tuple(shift for shift in range(power + 1) if shift & power == shift)
            columns[lookup[state]].append((row, shifts))
    return tuple(tuple(column) for column in columns)


@lru_cache(maxsize=32768)
def local_endpoint_image(modes, degrees, orders):
    states = tuple(product(range(2), repeat=modes))
    lookup = {state: index for index, state in enumerate(states)}
    variables = tuple((index, degree) for index, bound in enumerate(degrees)
                      for degree in range(max(0, bound + 1)))
    masks = [0] * len(states)
    equations = {}
    transform = difference_columns(modes)
    for column, (index, degree) in enumerate(variables):
        bit = 1 << column
        masks[index] |= bit
        for row, shifts in transform[index]:
            for shift in shifts:
                at = degree + shift
                if at < orders[row]:
                    equations[row, at] = equations.get((row, at), 0) ^ bit
    pivots = {}
    for equation in equations.values():
        while equation:
            pivot = equation.bit_length() - 1
            if pivot in pivots:
                equation ^= pivots[pivot]
            else:
                pivots[pivot] = equation
                break
    return variables, endpoint_image(variables, pivots, masks)


def endpoint_basis(modes, m, degrees, orders, construct=True):
    basis = {}
    for residue, local_degrees, local_orders in residue_patterns(m, degrees, orders):
        variables, image = local_endpoint_image(modes, local_degrees, local_orders)
        for endpoint, coefficients in image:
            lift = set()
            if construct:
                for column, (index, degree) in enumerate(variables):
                    if (coefficients >> column) & 1:
                        lift.add((index, residue + m * degree))
            while endpoint:
                pivot = endpoint.bit_length() - 1
                if pivot in basis:
                    old, old_lift = basis[pivot]
                    endpoint ^= old
                    if construct:
                        lift.symmetric_difference_update(old_lift)
                else:
                    basis[pivot] = endpoint, lift
                    break
    return tuple(basis.values())


def synthesize(modes, m, degrees, orders, construct=True, field_size=None):
    size = 2 ** modes
    degrees, orders = tuple(degrees), tuple(orders)
    field_size = 2 * size if field_size is None else field_size
    if (m < 1 or len(degrees) != size or len(orders) != size
            or any(not isinstance(value, int) or value < 0 for value in degrees + orders)):
        raise ValueError('Use positive m and one nonnegative degree/order bound per joint state.')
    if field_size <= size or field_size & (field_size - 1):
        raise ValueError('Use a binary extension with strictly more elements than endpoint coordinates.')
    basis = endpoint_basis(modes, m, degrees, orders, construct)
    support = 0
    for endpoint, lift in basis:
        support |= endpoint
    if support != (1 << size) - 1:
        return None
    if not construct:
        return True
    endpoint = [0] * size
    polynomial = {}
    for vector, lift in basis:
        union = tuple(j for j in range(size) if endpoint[j] or (vector >> j) & 1)
        if all(endpoint[j] for j in union):
            continue
        scalar = next(value for value in range(field_size)
                      if all(endpoint[j] ^ (value if (vector >> j) & 1 else 0) for j in union))
        for j in range(size):
            if (vector >> j) & 1:
                endpoint[j] ^= scalar
        for key in lift:
            polynomial[key] = polynomial.get(key, 0) ^ scalar
            if not polynomial[key]:
                del polynomial[key]
        if all(endpoint):
            break
    assert all(endpoint)
    states = tuple(product(range(2), repeat=modes))
    alpha = {state: {degree: value for (index, degree), value in polynomial.items() if index == i}
             for i, state in enumerate(states)}
    # beta_polynomials needs a prime-power n; the transform itself depends
    # only on m, so replay it here directly for arbitrary positive m.
    beta = {}
    for k, terms in finite_difference_terms(2, modes).items():
        value = {}
        for state, power, scalar in terms:
            for shift in range(power + 1):
                if shift & power != shift:
                    continue
                for degree, coefficient in alpha[state].items():
                    at = degree + m * shift
                    value[at] = value.get(at, 0) ^ coefficient
        beta[k] = {degree: coefficient for degree, coefficient in value.items() if coefficient}
    actual_degrees = tuple(max(alpha[state]) for state in states)
    actual_orders = tuple(min(beta[state]) for state in states)
    assert all(a <= b for a, b in zip(actual_degrees, degrees))
    assert all(a >= b for a, b in zip(actual_orders, orders))
    return {'modes': modes, 'm': m, 'coefficient_field_order': field_size,
            'degree_bounds': degrees, 'order_bounds': orders,
            'endpoint_image_dimension': len(basis), 'endpoint_values': tuple(endpoint),
            'alpha': alpha, 'beta': beta, 'actual_degrees': actual_degrees,
            'actual_orders': actual_orders,
            'residue_pattern_count': len(residue_patterns(m, degrees, orders))}
