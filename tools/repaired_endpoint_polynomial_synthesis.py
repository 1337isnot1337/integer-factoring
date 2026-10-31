#!/usr/bin/env python3
"""Polynomial synthesis with nonzero endpoints and a free diagonal repair.

The endpoint image need not contain the all-ones vector. It contains a
vector with every coordinate nonzero exactly when none of its coordinate
projections is zero. For p states over F_p no field extension is needed:
at most p proper coordinate hyperplanes cannot cover the image.
"""
from __future__ import annotations

import argparse
import json

from polynomial_state_filter_synthesis import (
    _sparse_add, endpoint_image, residue_patterns, transform,
)
from lucas_weighted_pascal import _check_prime


def endpoint_basis(prime, m, degrees, orders, construct=True):
    basis = []
    for residue, local_degrees, local_orders in residue_patterns(m, degrees, orders):
        variables, image = endpoint_image(prime, local_degrees, local_orders)
        for _, endpoint, coefficients in image:
            endpoint = list(endpoint)
            lift = {(j, residue + m * power): value
                    for (j, power), value in zip(variables, coefficients) if value} if construct else {}
            for pivot, vector, old_lift in basis:
                coefficient = endpoint[pivot]
                if coefficient:
                    endpoint = [(a - coefficient * b) % prime for a, b in zip(endpoint, vector)]
                    if construct:
                        lift = _sparse_add(lift, old_lift, -coefficient, prime)
            pivot = next((j for j, value in enumerate(endpoint) if value), None)
            if pivot is None:
                continue
            inverse = pow(endpoint[pivot], -1, prime)
            endpoint = tuple(value * inverse % prime for value in endpoint)
            lift = {key: value * inverse % prime for key, value in lift.items()}
            basis.append((pivot, endpoint, lift))
    return tuple(basis)


def synthesize(prime, m, degrees, orders, construct=True):
    _check_prime(prime)
    degrees, orders = tuple(degrees), tuple(orders)
    if (not isinstance(m, int) or m < 1 or len(degrees) != prime or len(orders) != prime
            or any(not isinstance(value, int) or value < 0 for value in degrees + orders)):
        raise ValueError('Use positive m and one nonnegative degree/order bound per state.')
    basis = endpoint_basis(prime, m, degrees, orders, construct)
    support = tuple(any(row[1][j] for row in basis) for j in range(prime))
    if not all(support):
        return None
    if not construct:
        return True
    # Increase endpoint support without losing an already nonzero entry.
    # If a basis vector is fully supported, use it directly. Otherwise
    # fewer than p scalars are forbidden when adjoining that vector.
    endpoint = [0] * prime
    lift = {}
    for _, vector, old_lift in basis:
        if all(endpoint):
            break
        if all(vector):
            endpoint, lift = list(vector), dict(old_lift)
            break
        union = [j for j in range(prime) if endpoint[j] or vector[j]]
        if all(endpoint[j] for j in union):
            continue
        coefficient = next(value for value in range(prime)
                           if all((endpoint[j] + value * vector[j]) % prime for j in union))
        endpoint = [(a + coefficient * b) % prime for a, b in zip(endpoint, vector)]
        lift = _sparse_add(lift, old_lift, coefficient, prime)
    assert all(endpoint)
    basis_lifts = (lift,)
    alpha_basis = tuple(tuple({degree: value for (digit, degree), value in lift.items() if digit == j}
                              for j in range(prime)) for lift in basis_lifts)
    beta_basis = tuple(transform(alpha, m, prime) for alpha in alpha_basis)
    actual_degrees = tuple(max(degree for alpha in alpha_basis for degree in alpha[j])
                           for j in range(prime))
    actual_orders = tuple(min(degree for beta in beta_basis for degree in beta[k])
                          for k in range(prime))
    endpoint_values = tuple(tuple(sum(alpha[j].values()) % prime for alpha in alpha_basis)
                            for j in range(prime))
    assert all(any(value) for value in endpoint_values)
    assert all(a <= b for a, b in zip(actual_degrees, degrees))
    assert all(a >= b for a, b in zip(actual_orders, orders))
    return {
        'prime': prime, 'm': m, 'degree_bounds': degrees, 'order_bounds': orders,
        'endpoint_image_dimension': len(basis), 'extension_degree': len(alpha_basis),
        'coefficient_field': f'F_{prime}',
        'alpha_basis': alpha_basis, 'beta_basis': beta_basis,
        'endpoint_values_in_extension_basis': endpoint_values,
        'actual_degrees': actual_degrees, 'actual_orders': actual_orders,
        'residue_pattern_count': len(residue_patterns(m, degrees, orders)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=3)
    parser.add_argument('--m', type=int, default=16)
    parser.add_argument('--degrees', default='0,0,1')
    parser.add_argument('--orders', default='1,0,1')
    args = parser.parse_args()
    result = synthesize(args.prime, args.m, tuple(map(int, args.degrees.split(','))),
                        tuple(map(int, args.orders.split(','))))
    print(json.dumps({'feasible_with_nonzero_endpoint_repair': result is not None,
                      'witness': result}, indent=2))


if __name__ == '__main__':
    main()
