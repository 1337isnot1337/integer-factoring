#!/usr/bin/env python3
"""Terminal and mixed-extreme polynomial identity filters at prime-power n.

The original Q_j remain polynomial. Only zero-cut basis changes use
zero units; any denominator in an active form is cleared by a regular
change of the active V variables before the exact Smith calculation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import comb, factorial

from coupled_identity_mode_filters import (
    add, algebra_multiply, inverse_matrix, monomial_filters, multiply,
    multiplier_coefficients, normalized_filters, shift, zero_charge_blocks,
)
from prime_state_filter_profiles import validate


def is_power(value, prime):
    while value > 1 and value % prime == 0:
        value //= prime
    return value == 1


def reverse_form(i, offset, prime):
    """f_(n-1-offset)(i,c)=tau^i f_offset(i,-c/tau), n=p^k."""
    result = {}
    for h in range(min(i, offset) + 1):
        scalar = comb(offset, h) * comb(i, h)
        for power in range(h + 1):
            exponent = i - h + power
            coefficient = scalar * comb(h, power) * (-1) ** (h + power)
            result[exponent] = (result.get(exponent, 0) + coefficient) % prime
    return {power: value for power, value in result.items() if value}


@lru_cache(maxsize=128)
def geometry_data(prime, n, modes=2, geometry='mixed'):
    validate(prime, n, 1, (0,) * prime)
    if not is_power(n, prime):
        raise ValueError('The terminal Lucas identities require n to be a power of the characteristic.')
    if not isinstance(modes, int) or not 1 <= modes <= n:
        raise ValueError('Select between one and n identity modes.')
    if geometry == 'terminal':
        selected = tuple(range(n - modes, n))
        retained = tuple(range(modes))
        discarded = tuple(range(modes, n))
        forms = tuple(tuple(reverse_form(i, n - 1 - j, prime) for j in selected)
                      for i in range(n))
        vandermonde = [[comb(i, h) % prime for h in range(modes)] for i in retained]
        inverse = inverse_matrix(vandermonde, prime)
        interpolation = tuple(tuple(sum(comb(i, h) * inverse[h][j]
                                        for h in range(modes)) % prime
                                    for j in range(modes)) for i in discarded)
        f = tuple(forms[i] for i in retained)
        h = []
        for a in range(modes):
            row = []
            for j in range(modes):
                value = dict(f[a][j])
                for k, i in enumerate(discarded):
                    value = add(value, forms[i][j], prime, interpolation[k][a])
                row.append(value)
            h.append(tuple(row))
        unit = {0: 1}
        return {'prime': prime, 'n': n, 'modes': modes, 'geometry': geometry,
                'selected': selected, 'retained': retained, 'discarded': discarded,
                'all_forms': forms, 'interpolation_constants': interpolation,
                'active_x': f, 'active_v': tuple(h), 'gaussian_v_unit': unit}
    if geometry != 'mixed' or modes != 2 or n < 2:
        raise ValueError('Mixed geometry uses exactly two modes and n>=2.')
    t = n - 1
    selected = retained = (0, t)
    discarded = tuple(range(1, t))
    forms = tuple(({0: 1}, {i: 1}) for i in range(n))
    f = (forms[0], forms[t])
    unit = {0: 1, t: -1 % prime}
    numerator = tuple(({i: 1, t: -1 % prime}, {0: 1, i: -1 % prime})
                      for i in discarded)
    h = []
    for a, retained_index in enumerate(retained):
        row = []
        for j in range(modes):
            value = multiply(unit, f[a][j], prime)
            for k, i in enumerate(discarded):
                term = multiply(numerator[k][a], forms[i][j], prime)
                value = add(value, shift(term, retained_index - i, prime), prime)
            assert all(power >= 0 for power in value)
            row.append(value)
        h.append(tuple(row))
    return {'prime': prime, 'n': n, 'modes': modes, 'geometry': geometry,
            'selected': selected, 'retained': retained, 'discarded': discarded,
            'all_forms': forms, 'interpolation_numerators': numerator,
            'active_x': f, 'active_v': tuple(h), 'gaussian_v_unit': unit}


@lru_cache(maxsize=32)
def active_power_kernels(prime, n, ramification=1, modes=2, geometry='mixed'):
    validate(prime, n, ramification, (0,) * prime)
    data = geometry_data(prime, n, modes, geometry)
    one = {0: 1}
    zero = (0,) * (4 * modes)
    gaussian = {zero: one}
    for j in range(modes):
        a = 2 * ramification * data['selected'][j]
        b = 2 * ramification * (data['retained'][j] + 1)
        for left, right, order, unit in (
                (modes + j, 2 * modes + j, a, one),
                (j, 3 * modes + j, b,
                 {2 * ramification * power: value
                  for power, value in data['gaussian_v_unit'].items()})):
            factor = {}
            coefficient = one
            for digit in range(prime):
                powers = list(zero)
                powers[left] = powers[right] = digit
                factor[tuple(powers)] = shift(coefficient, order * digit, prime,
                                              pow(factorial(digit), -1, prime))
                coefficient = multiply(coefficient, unit, prime)
            gaussian = algebra_multiply(gaussian, factor, prime)
    q_modes = []
    for j in range(modes):
        a = 2 * ramification * data['selected'][j]
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


def active_tensor(prime, n, ramification, modes, alpha, geometry='mixed'):
    kernels = active_power_kernels(prime, n, ramification, modes, geometry)
    scalar = multiplier_coefficients(prime, modes, alpha)
    result = {}
    for powers, coefficient in scalar.items():
        coefficient = {2 * degree: value for degree, value in coefficient.items()}
        for variables, value in kernels[powers].items():
            result[variables] = add(result.get(variables, {}),
                                    multiply(coefficient, value, prime), prime)
    return {powers: value for powers, value in result.items() if value}


def zero_profiles(prime, n, ramification, modes, alpha, geometry='mixed'):
    from check_all_field_pascal_crossings import smith_orders

    tensor = active_tensor(prime, n, ramification, modes, alpha, geometry)
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


def profile_data(prime, n, ramification, charge, modes, alpha, geometry='mixed', blocks=None):
    validate(prime, n, ramification, (0,) * prime)
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an integer or Fraction for an exact physical charge.')
    data = geometry_data(prime, n, modes, geometry)
    alpha = normalized_filters(prime, modes, alpha)
    if blocks is None:
        blocks = zero_profiles(prime, n, ramification, modes, alpha, geometry)
    kappa = tuple(ramification * (2 * j + 1 - 2 * n) for j in data['selected'])
    return {
        'prime': prime, 'n': n, 'ramification': ramification,
        'modes': modes, 'charge': charge, 'geometry': geometry,
        'identity_indices': data['selected'], 'retained_xv_indices': data['retained'],
        'zero_active': tuple(sorted(value + ramification * (1 - charge) * block['charge']
                                    for block in blocks for value in block['orders'])),
        'zero_charge_blocks_at_J1': blocks,
        'infinity_active': tuple(sorted(sum(a * b for a, b in zip(kappa, state))
                                        - 2 * max(value) for state, value in alpha.items())),
        'zero_spectators': tuple(ramification * (2 * i + 1 - charge)
                                 for i in range(n) if i not in data['selected'])
        + tuple(ramification * (2 * i + 1 + charge)
                for i in range(n) if i not in data['retained']),
        'infinity_spectators': tuple(ramification * (2 * i + 1 - 2 * n)
                                     for i in range(2 * n) if i not in data['selected']),
        'paired_physical_first': tuple(ramification * (4 * i - 2 * (n - 1) + charge)
                                       for i in range(n)),
        'paired_physical_second': tuple(ramification * (4 * i - 2 * (n - 1) - charge)
                                        for i in range(n)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=2)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--charge', type=Fraction, default=Fraction(1))
    parser.add_argument('--geometry', choices=('terminal', 'mixed'), default='mixed')
    parser.add_argument('--exponents', default='0,0,0,1')
    args = parser.parse_args()
    alpha = monomial_filters(args.prime, args.modes, tuple(map(int, args.exponents.split(','))))
    data = profile_data(args.prime, args.n, args.ramification, args.charge,
                        args.modes, alpha, args.geometry)
    print(json.dumps(data, indent=2, default=str))


if __name__ == '__main__':
    main()
