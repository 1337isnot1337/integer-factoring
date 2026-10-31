#!/usr/bin/env python3
"""Exact bounded cores for correlated filters on polynomial identity modes.

For r selected identity modes and n>=r, regular Gaussian-preserving
changes leave 2r active zero modes and r active infinity modes. The
state weights may be genuinely coupled across all r occupations.
Sparse scalar polynomials use x=-v; zero coefficients use v=-z^2.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import comb, factorial

from prime_state_filter_profiles import validate


def add(first, second, prime, coefficient=1):
    result = dict(first)
    for power, value in second.items():
        result[power] = (result.get(power, 0) + coefficient * value) % prime
        if not result[power]:
            del result[power]
    return result


def multiply(first, second, prime):
    result = {}
    for a, x in first.items():
        for b, y in second.items():
            result[a + b] = (result.get(a + b, 0) + x * y) % prime
    return {power: value for power, value in result.items() if value}


def shift(value, amount, prime, coefficient=1):
    return {power + amount: x * coefficient % prime for power, x in value.items()
            if x * coefficient % prime}


def inverse_matrix(matrix, prime):
    size = len(matrix)
    work = [[value % prime for value in row]
            + [int(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i][column])
        work[column], work[pivot] = work[pivot], work[column]
        inverse = pow(work[column][column], -1, prime)
        work[column] = [value * inverse % prime for value in work[column]]
        for i in range(size):
            if i != column and work[i][column]:
                coefficient = work[i][column]
                work[i] = [(a - coefficient * b) % prime
                           for a, b in zip(work[i], work[column])]
    return tuple(tuple(row[size:]) for row in work)


@lru_cache(maxsize=256)
def concentration_data(prime, n, modes=2):
    """Return the exact interpolation and active form matrices in tau=z^(2P).

    F is the last-r-row interpolation matrix, B interpolates its earlier
    rows, and H=W F is the transformed coefficient matrix on V. No
    local-field inverses appear in the returned active coefficients.
    """
    validate(prime, n, 1, (0,) * prime)
    if not isinstance(modes, int) or not 1 <= modes <= n:
        raise ValueError('Select between one and n initial identity modes.')
    top = tuple(range(n - modes, n))
    vandermonde = [[comb(i, j) % prime for j in range(modes)] for i in top]
    inverse = inverse_matrix(vandermonde, prime)
    interpolation = tuple(tuple(sum(comb(i, k) * inverse[k][j]
                                    for k in range(modes)) % prime
                                for j in range(modes)) for i in range(n - modes))

    def f(i, j):
        result = {}
        for h in range(j + 1):
            scalar = comb(j, h) * comb(i, h)
            for power in range(h + 1):
                coefficient = scalar * comb(h, power) * (-1) ** power
                result[power] = (result.get(power, 0) + coefficient) % prime
        return {power: value for power, value in result.items() if value}

    all_forms = tuple(tuple(f(i, j) for j in range(modes)) for i in range(n))
    active_x = tuple(all_forms[i] for i in top)
    active_v = []
    for h, t in enumerate(top):
        row = []
        for j in range(modes):
            value = dict(active_x[h][j])
            for i in range(n - modes):
                value = add(value, shift(all_forms[i][j], t - i, prime),
                            prime, interpolation[i][h])
            row.append(value)
        active_v.append(tuple(row))
    return {'prime': prime, 'n': n, 'modes': modes, 'top': top,
            'interpolation': interpolation, 'all_forms': all_forms,
            'active_x': active_x, 'active_v': tuple(active_v)}


def algebra_multiply(first, second, prime):
    result = {}
    for powers, x in first.items():
        for other, y in second.items():
            target = tuple(a + b for a, b in zip(powers, other))
            if max(target) >= prime:
                continue
            result[target] = add(result.get(target, {}), multiply(x, y, prime), prime)
    return {powers: value for powers, value in result.items() if value}


@lru_cache(maxsize=64)
def active_power_kernels(prime, n, ramification=1, modes=2):
    """G times every product of Q_j powers at the base charge J=1."""
    validate(prime, n, ramification, (0,) * prime)
    data = concentration_data(prime, n, modes)
    one = {0: 1}
    gaussian = {((0,) * (4 * modes)): one}
    zero = (0,) * (4 * modes)
    for j in range(modes):
        a = 2 * ramification * j
        b = 2 * ramification * (data['top'][j] + 1)
        for left, right, order in ((modes + j, 2 * modes + j, a),
                                   (j, 3 * modes + j, b)):
            factor = {}
            for digit in range(prime):
                powers = list(zero)
                powers[left] = powers[right] = digit
                factor[tuple(powers)] = {order * digit: pow(factorial(digit), -1, prime)}
            gaussian = algebra_multiply(gaussian, factor, prime)
    q_modes = []
    for j in range(modes):
        a = 2 * ramification * j
        powers = list(zero)
        powers[modes + j] = 1
        left = {tuple(powers): {a: 1}}
        powers = list(zero)
        powers[2 * modes + j] = 1
        right = {tuple(powers): one}
        for h in range(modes):
            powers = list(zero)
            powers[h] = 1
            coefficient = multiply({0: 1, 1: -1 % prime}, data['active_v'][h][j], prime)
            left[tuple(powers)] = {2 * ramification * power: value
                                   for power, value in coefficient.items()}
            powers = list(zero)
            powers[3 * modes + h] = 1
            right[tuple(powers)] = {2 * ramification * power - a: -value % prime
                                    for power, value in data['active_x'][h][j].items()}
        q_modes.append(algebra_multiply(left, right, prime))
    kernels = {(0,) * modes: gaussian}
    for state in product(range(prime), repeat=modes):
        if not any(state):
            continue
        j = next(index for index, value in enumerate(state) if value)
        previous = list(state)
        previous[j] -= 1
        kernels[state] = algebra_multiply(q_modes[j], kernels[tuple(previous)], prime)
    return kernels


def normalized_filters(prime, modes, alpha):
    states = tuple(product(range(prime), repeat=modes))
    if set(alpha) != set(states):
        raise ValueError('Supply a scalar filter for every joint occupation state.')
    result = {}
    for state in states:
        value = {power: coefficient % prime for power, coefficient in alpha[state].items()
                 if coefficient % prime}
        if any(not isinstance(power, int) for power in value):
            raise ValueError('Scalar Laurent powers in x=-v must be integers.')
        if not value or sum(value.values()) % prime != 1:
            raise ValueError('Every joint occupation filter must equal one at x=1.')
        result[state] = value
    return result


def multiplier_coefficients(prime, modes, alpha):
    alpha = normalized_filters(prime, modes, alpha)
    result = {}
    for powers in product(range(prime), repeat=modes):
        value = {}
        for state in product(*(range(power + 1) for power in powers)):
            coefficient = (-1) ** sum(a - b for a, b in zip(powers, state))
            for a, b in zip(powers, state):
                coefficient *= pow(factorial(b) * factorial(a - b), -1, prime)
            value = add(value, alpha[state], prime, coefficient)
        if value:
            result[powers] = value
    return result


def active_tensor(prime, n, ramification, modes, alpha):
    kernels = active_power_kernels(prime, n, ramification, modes)
    scalar = multiplier_coefficients(prime, modes, alpha)
    result = {}
    for powers, coefficient in scalar.items():
        coefficient = {2 * degree: value for degree, value in coefficient.items()}
        for variables, value in kernels[powers].items():
            result[variables] = add(result.get(variables, {}),
                                    multiply(coefficient, value, prime), prime)
    return {powers: value for powers, value in result.items() if value}


def zero_charge_blocks(prime, modes, tensor):
    states = tuple(product(range(prime), repeat=modes))
    result = []
    for charge in range(-(prime - 1) * modes, (prime - 1) * modes + 1):
        labels = [(y, x) for y, x in product(states, repeat=2)
                  if sum(y) - sum(x) == charge]
        matrix = []
        for y, x in labels:
            row = []
            for u, v in labels:
                row.append(tensor.get(v + y + u + x, {}))
            matrix.append(row)
        result.append((charge, matrix))
    return result


def zero_profiles(prime, n, ramification, modes, alpha):
    """Exact complete active profile, organized by common charge at J=1."""
    from check_all_field_pascal_crossings import smith_orders

    tensor = active_tensor(prime, n, ramification, modes, alpha)
    result = []
    for charge, matrix in zero_charge_blocks(prime, modes, tensor):
        minimum = min(power for row in matrix for value in row for power in value)
        maximum = max(power for row in matrix for value in row for power in value)
        dense = [[tuple(value.get(power + minimum, 0)
                        for power in range(max(value, default=minimum) - minimum + 1))
                  if value else () for value in row] for row in matrix]
        bound = max(512, 2 * len(matrix) * (maximum - minimum + 1))
        orders, precision = smith_orders(dense, prime, maximum_precision=bound)
        result.append({'charge': charge, 'orders': tuple(value + minimum for value in orders),
                       'exact_jet_precision': precision})
    return tuple(result)


def profile_data(prime, n, ramification, charge, modes, alpha, blocks=None):
    validate(prime, n, ramification, (0,) * prime)
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an integer or Fraction for an exact physical charge.')
    alpha = normalized_filters(prime, modes, alpha)
    if blocks is None:
        blocks = zero_profiles(prime, n, ramification, modes, alpha)
    kappa = tuple(ramification * (2 * j + 1 - 2 * n) for j in range(modes))
    result = {
        'prime': prime, 'n': n, 'ramification': ramification,
        'modes': modes, 'charge': charge,
        'zero_active': tuple(sorted(value + ramification * (1 - charge) * block['charge']
                                    for block in blocks for value in block['orders'])),
        'zero_charge_blocks_at_J1': blocks,
        'infinity_active': tuple(sorted(sum(a * b for a, b in zip(kappa, state))
                                        - 2 * max(value) for state, value in alpha.items())),
        'zero_spectators': tuple(ramification * (2 * i + 1 - charge) for i in range(modes, n))
        + tuple(ramification * (2 * i + 1 + charge) for i in range(n - modes)),
        'infinity_spectators': tuple(ramification * (2 * i + 1 - 2 * n)
                                     for i in range(modes, 2 * n)),
        'paired_physical_first': tuple(ramification * (4 * i - 2 * (n - 1) + charge)
                                       for i in range(n)),
        'paired_physical_second': tuple(ramification * (4 * i - 2 * (n - 1) - charge)
                                        for i in range(n)),
    }
    return result


def monomial_filters(prime, modes, exponents):
    states = tuple(product(range(prime), repeat=modes))
    if len(exponents) != len(states):
        raise ValueError('Supply one Laurent exponent per joint occupation state.')
    return {state: {-exponent: 1} for state, exponent in zip(states, exponents)}


def polynomial_seed(prime, n, ramification, modes, alpha, identity_indices=None):
    """Construct and monomial-clear the actual original polynomial tensor.

    This reference path uses every physical mode and is exponential in n.
    Return (tensor, clearing), with raw coefficient normalization. The
    clearing x^h adds +2h at zero and -2h at infinity, and is one at x=1.
    """
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel

    validate(prime, n, ramification, (0,) * prime)
    concentration_data(prime, n, modes)
    if identity_indices is None:
        identity_indices = tuple(range(modes))
    identity_indices = tuple(identity_indices)
    if (len(identity_indices) != modes or len(set(identity_indices)) != modes
            or any(not isinstance(index, int) or not 0 <= index < n for index in identity_indices)):
        raise ValueError('Select distinct physical identity-mode indices.')
    alpha = normalized_filters(prime, modes, alpha)
    clearing = max(0, -min(power for value in alpha.values() for power in value))
    matrix = pascal_kernel(tuple(range(n)), 2 * n, prime)
    gaussian = {powers: {degree: value for degree, value in enumerate(polynomial) if value}
                for powers, polynomial in direct_kernel(matrix, prime).items()}
    zero = (0,) * (4 * n)
    q_modes = []
    for j in identity_indices:
        powers = list(zero)
        powers[n + j] = 1
        left = {tuple(powers): {0: 1}}
        powers = list(zero)
        powers[2 * n + j] = 1
        right = {tuple(powers): {0: 1}}
        for i in range(n):
            powers = list(zero)
            powers[i] = 1
            value = {degree: coefficient for degree, coefficient in enumerate(matrix[i][j])
                     if coefficient}
            if value:
                left[tuple(powers)] = value
            powers = list(zero)
            powers[3 * n + i] = 1
            value = {degree: coefficient for degree, coefficient in enumerate(matrix[n + j][n + i])
                     if coefficient}
            if value:
                right[tuple(powers)] = value
        q_modes.append(algebra_multiply(left, right, prime))
    scalar = multiplier_coefficients(prime, modes, alpha)
    kernels = {(0,) * modes: gaussian}
    for state in product(range(prime), repeat=modes):
        if not any(state):
            continue
        j = next(index for index, value in enumerate(state) if value)
        previous = list(state)
        previous[j] -= 1
        kernels[state] = algebra_multiply(q_modes[j], kernels[tuple(previous)], prime)
    tensor = {}
    for state, coefficient in scalar.items():
        coefficient = {clearing + degree: value * (-1) ** (clearing + degree) % prime
                       for degree, value in coefficient.items()}
        for powers, value in kernels[state].items():
            value = {ramification * degree: x for degree, x in value.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(value, coefficient, prime), prime)
    tensor = {powers: tuple(value.get(degree, 0) for degree in range(max(value) + 1))
              for powers, value in tensor.items() if value}
    return tensor, clearing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=2)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--charge', type=Fraction, default=Fraction(1))
    parser.add_argument('--exponents', default='0,0,0,1')
    args = parser.parse_args()
    exponents = tuple(map(int, args.exponents.split(',')))
    alpha = monomial_filters(args.prime, args.modes, exponents)
    data = profile_data(args.prime, args.n, args.ramification, args.charge, args.modes, alpha)
    print(json.dumps(data, indent=2, default=str))


if __name__ == '__main__':
    main()
