#!/usr/bin/env python3
"""Exact polynomial state-filter synthesis with residue compression.

Find alpha_j(q) of prescribed degree bounds, alpha_j(1)=1, whose
finite-difference beta_k have prescribed zero orders. The transform
depends on q only through u=q^m. Residue classes modulo m separate,
and at most 2*p+1 distinct bound patterns are needed, independent of m.

All arithmetic is over F_p. A returned witness is verified directly as
sparse polynomials; feasibility is not inferred from floating pressure.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
from math import comb, factorial

from lucas_weighted_pascal import _check_prime


def _rref(matrix, columns, prime):
    work = [list(row) for row in matrix]
    pivots = []
    rank = 0
    for column in range(columns):
        selected = next((i for i in range(rank, len(work)) if work[i][column] % prime), None)
        if selected is None:
            continue
        work[rank], work[selected] = work[selected], work[rank]
        inverse = pow(work[rank][column] % prime, -1, prime)
        work[rank] = [value * inverse % prime for value in work[rank]]
        for i in range(len(work)):
            if i == rank or not work[i][column]:
                continue
            coefficient = work[i][column]
            work[i] = [(a - coefficient * b) % prime
                       for a, b in zip(work[i], work[rank])]
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return work[:rank], tuple(pivots)


@lru_cache(maxsize=None)
def difference_matrix(prime):
    _check_prime(prime)
    result = [[{} for _ in range(prime)] for _ in range(prime)]
    for k in range(prime):
        for j in range(k + 1):
            alpha = prime - 1 - k + j
            scalar = (-1) ** (k - j) * comb(k, j)
            result[k][alpha] = {
                degree: scalar * comb(j, degree) * (-1) ** degree % prime
                for degree in range(j + 1)
                if scalar * comb(j, degree) % prime
            }
    return tuple(tuple(row) for row in result)


@lru_cache(maxsize=200000)
def endpoint_image(prime, degrees, orders):
    """Endpoint subspace and exact coefficient lifts for one residue pattern."""
    variables = tuple((j, degree) for j, bound in enumerate(degrees)
                      for degree in range(max(0, bound + 1)))
    matrix = difference_matrix(prime)
    constraints = []
    for k, order in enumerate(orders):
        maximum = max((alpha_degree + power
                       for j, alpha_degree in variables for power in matrix[k][j]),
                      default=-1)
        # Coefficients above this degree vanish identically, even when an
        # impossible order target is given as a very large integer.
        for degree in range(min(order, maximum + 1)):
            constraints.append([
                matrix[k][j].get(degree - alpha_degree, 0)
                for j, alpha_degree in variables
            ])
    reduced, pivots = _rref(constraints, len(variables), prime)
    pivot_set = set(pivots)
    image = []
    for free in range(len(variables)):
        if free in pivot_set:
            continue
        lift = [0] * len(variables)
        lift[free] = 1
        for row, pivot in zip(reduced, pivots):
            lift[pivot] = -row[free] % prime
        endpoint = [0] * prime
        for coefficient, (j, degree) in zip(lift, variables):
            endpoint[j] = (endpoint[j] + coefficient) % prime
        for pivot, vector, old_lift in image:
            coefficient = endpoint[pivot]
            if coefficient:
                endpoint = [(a - coefficient * b) % prime for a, b in zip(endpoint, vector)]
                lift = [(a - coefficient * b) % prime for a, b in zip(lift, old_lift)]
        pivot = next((i for i, value in enumerate(endpoint) if value), None)
        if pivot is None:
            continue
        inverse = pow(endpoint[pivot], -1, prime)
        endpoint = tuple(value * inverse % prime for value in endpoint)
        lift = tuple(value * inverse % prime for value in lift)
        image.append((pivot, endpoint, lift))
    return variables, tuple(image)


def residue_patterns(m, degrees, orders):
    cuts = {0, m}
    for degree in degrees:
        cut = degree % m + 1
        if cut < m:
            cuts.add(cut)
    for order in orders:
        cut = order % m
        if cut:
            cuts.add(cut)
    cuts = sorted(cuts)
    result = []
    for residue in cuts[:-1]:
        local_degrees = tuple((degree - residue) // m for degree in degrees)
        local_orders = tuple(max(0, (order - residue + m - 1) // m) for order in orders)
        result.append((residue, local_degrees, local_orders))
    return tuple(result)


def _sparse_add(first, second, coefficient, prime):
    result = dict(first)
    for key, value in second.items():
        result[key] = (result.get(key, 0) + coefficient * value) % prime
        if not result[key]:
            del result[key]
    return result


def transform(alpha, m, prime):
    result = []
    for row in difference_matrix(prime):
        polynomial = {}
        for value, multiplier in zip(alpha, row):
            for degree, coefficient in value.items():
                for power, scalar in multiplier.items():
                    at = degree + m * power
                    polynomial[at] = (polynomial.get(at, 0) + coefficient * scalar) % prime
        result.append({degree: value for degree, value in polynomial.items() if value})
    return tuple(result)


def synthesize(prime, m, degrees, orders, construct=True):
    """Return None if impossible; otherwise return a verified sparse witness.

    With construct=False, return True/None without exporting the lift.
    The endpoint target is the all-ones vector in F_p^p.
    """
    _check_prime(prime)
    degrees = tuple(degrees)
    orders = tuple(orders)
    if not isinstance(m, int) or m < 1:
        raise ValueError('The mode scale m must be a positive integer.')
    if len(degrees) != prime or len(orders) != prime:
        raise ValueError('One degree bound and zero-order target per state are required.')
    if any(not isinstance(value, int) or value < 0 for value in degrees + orders):
        raise ValueError('Degree bounds and zero orders must be nonnegative integers.')
    basis = []
    patterns = residue_patterns(m, degrees, orders)
    for residue, local_degrees, local_orders in patterns:
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
            pivot = next((i for i, value in enumerate(endpoint) if value), None)
            if pivot is None:
                continue
            inverse = pow(endpoint[pivot], -1, prime)
            endpoint = tuple(value * inverse % prime for value in endpoint)
            if construct:
                lift = {key: value * inverse % prime for key, value in lift.items()}
            basis.append((pivot, endpoint, lift))

    target = [1] * prime
    witness = {}
    for pivot, vector, lift in basis:
        coefficient = target[pivot]
        if coefficient:
            target = [(a - coefficient * b) % prime for a, b in zip(target, vector)]
            if construct:
                witness = _sparse_add(witness, lift, coefficient, prime)
    if any(target):
        return None
    if not construct:
        return True
    alpha = tuple({degree: value for (j, degree), value in witness.items() if j == i}
                  for i in range(prime))
    beta = transform(alpha, m, prime)
    assert all(value and max(value) <= bound and sum(value.values()) % prime == 1
               for value, bound in zip(alpha, degrees))
    assert all(value and min(value) >= bound for value, bound in zip(beta, orders))
    return {
        'prime': prime, 'm': m, 'degree_bounds': degrees, 'order_bounds': orders,
        'alpha': alpha, 'beta': beta,
        'actual_degrees': tuple(max(value) for value in alpha),
        'actual_orders': tuple(min(value) for value in beta),
        'residue_pattern_count': len(patterns),
    }


def profile_data(prime, n, ramification, charge, alpha, include_spectators=True):
    """Complete profiles for an endpoint-one polynomial state filter.

    Polynomial exponents are in q=-v; the physical substitution v=-z^2
    doubles them. Infinity uses every polynomial's full degree. Scalar
    common factors are included, not removed without accounting for them.
    """
    from prime_state_filter_profiles import active_profiles

    data = active_profiles(prime, n, ramification, charge,
                           exponents=(0,) * prime,
                           include_spectators=include_spectators)
    alpha = tuple({degree: value % prime for degree, value in polynomial.items()
                   if value % prime} for polynomial in alpha)
    if len(alpha) != prime or any(not value for value in alpha):
        raise ValueError('Supply one nonzero sparse polynomial for every digit.')
    if any(not isinstance(degree, int) or degree < 0
           for value in alpha for degree in value):
        raise ValueError('Polynomial powers in q=-v must be nonnegative integers.')
    if any(sum(value.values()) % prime != 1 for value in alpha):
        raise ValueError('Every digit filter must equal one at the crossing endpoint q=1.')
    beta = transform(alpha, n * ramification, prime)
    orders = tuple(2 * min(value) for value in beta)
    degrees = tuple(max(value) for value in alpha)
    a = ramification * (1 - charge)
    kappa = ramification * (1 - 2 * n)
    data.pop('exponents')
    data.update({
        'alpha': alpha, 'beta': beta, 'alpha_degrees': degrees,
        'beta_orders': orders,
        'zero_active': tuple(sorted(orders[k] + a * (digit - k)
                                    for k in range(prime) for digit in range(prime))),
        'zero_active_first_factor': tuple(a * digit for digit in range(prime)),
        'zero_active_second_factor': tuple(orders[k] - a * k for k in range(prime)),
        'infinity_active': tuple(sorted(kappa * digit - 2 * degree
                                        for digit, degree in enumerate(degrees))),
    })
    return data


def polynomial_seed(prime, n, ramification, alpha):
    """Literal finite polynomial tensor; intended for small exact controls.

    The state-filter multiplier is Exp_p(-Q) times
    sum_j alpha_j(q)*Q^j/j!. This constructor charges the ordinary
    v degree. Its cost is exponential in the interval length.
    Coefficients have raw truncated-exponential normalization in
    physical variable order (V,Y,U,X), as in the predecessor tools.
    """
    from check_all_field_pascal_crossings import add, multiply, pascal_kernel, trim
    from check_truncated_symmetric_crossings import direct_kernel

    profile_data(prime, n, ramification, 0, alpha, include_spectators=False)

    def kernel_product(first, second):
        result = {}
        for powers, x in first.items():
            for other, y in second.items():
                target = tuple(a + b for a, b in zip(powers, other))
                if max(target) >= prime:
                    continue
                result[target] = add(result.get(target, ()), multiply(x, y, prime), prime)
        return {powers: value for powers, value in result.items() if value}

    gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, prime), prime)
    left, right = {}, {}
    for index, value in ((n, (1,)), (0, (1, 1))):
        powers = [0] * (4 * n)
        powers[index] = 1
        left[tuple(powers)] = value
    powers = [0] * (4 * n)
    powers[2 * n] = 1
    right[tuple(powers)] = (1,)
    for b in range(n):
        powers = [0] * (4 * n)
        powers[3 * n + b] = 1
        right[tuple(powers)] = (0,) * (2 * n - 1 - b) + ((-1) ** b % prime,)
    q = kernel_product(left, right)
    current = gaussian
    result = {}
    for k in range(prime):
        scalar = {}
        for j in range(k + 1):
            coefficient = (-1) ** (k - j) * pow(factorial(j) * factorial(k - j), -1, prime)
            scalar = _sparse_add(scalar, alpha[j], coefficient, prime)
        for powers, polynomial in current.items():
            if not scalar:
                continue
            degree = max(scalar) + ramification * (len(polynomial) - 1)
            value = [0] * (degree + 1)
            for exponent, coefficient in scalar.items():
                for i, x in enumerate(polynomial):
                    at = exponent + ramification * i
                    value[at] = (value[at] + (-1) ** exponent * coefficient * x) % prime
            result[powers] = add(result.get(powers, ()), trim(value, prime), prime)
        if k + 1 < prime:
            current = kernel_product(q, current)
    return {powers: value for powers, value in result.items() if value}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=3)
    parser.add_argument('--m', type=int, default=16)
    parser.add_argument('--degrees', default='0,0,0')
    parser.add_argument('--orders', default='0,16,32')
    args = parser.parse_args()
    result = synthesize(args.prime, args.m,
                        tuple(map(int, args.degrees.split(','))),
                        tuple(map(int, args.orders.split(','))))
    print(json.dumps({'feasible': result is not None, 'witness': result}, indent=2))


if __name__ == '__main__':
    main()
