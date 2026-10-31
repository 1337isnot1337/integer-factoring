#!/usr/bin/env python3
"""Exact profiles for state filters with arbitrary nonzero endpoint values.

At w=-1, Q_j=Y_j(U_j+C_jX), so occupation filtering multiplies each
endpoint coefficient by the value attached to its physical Y label.
Its inverse is a constant diagonal physical-Y repair. This costs no
degree or weighted Smith order, and restores the original endpoint.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb

from coupled_identity_extreme_filters import is_power
from prime_state_filter_profiles import validate


def beta_polynomials(prime, n, ramification, exponents, coefficients):
    validate(prime, n, ramification, exponents)
    if len(coefficients) != prime or any(not isinstance(value, int) or not value % prime
                                        for value in coefficients):
        raise ValueError('Supply one nonzero prime-field coefficient per digit.')
    result = []
    for k in range(prime):
        value = {}
        for j in range(k + 1):
            index = prime - 1 - k + j
            scalar = coefficients[index] * (-1) ** (k - j) * comb(k, j)
            for h in range(j + 1):
                degree = n * ramification * h - exponents[index]
                coefficient = scalar * (-1) ** h * comb(j, h)
                value[degree] = (value.get(degree, 0) + coefficient) % prime
        value = {degree: coefficient for degree, coefficient in value.items() if coefficient}
        if not value:
            raise ArithmeticError('A finite difference cannot vanish with a nonzero endpoint value.')
        result.append(value)
    return tuple(result)


def profile_data(prime, n, ramification, charge, exponents, coefficients,
                 geometry='initial', include_spectators=True):
    beta = beta_polynomials(prime, n, ramification, exponents, coefficients)
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an integer or Fraction for the exact aggregate charge.')
    if geometry == 'initial':
        selected, retained = 0, n - 1
    elif geometry == 'terminal' and is_power(n, prime):
        selected, retained = n - 1, 0
    else:
        raise ValueError('Terminal scalar geometry requires a prime-power interval length.')
    a = ramification * (2 * selected + 1 - charge)
    kappa = ramification * (2 * selected + 1 - 2 * n)
    orders = tuple(2 * min(value) for value in beta)
    result = {
        'prime': prime, 'n': n, 'ramification': ramification, 'charge': charge,
        'geometry': geometry, 'identity_index': selected, 'retained_xv_index': retained,
        'exponents': tuple(exponents), 'endpoint_values': tuple(value % prime for value in coefficients),
        'beta_orders': orders,
        'zero_active': tuple(sorted(orders[k] + a * (digit - k)
                                    for digit in range(prime) for k in range(prime))),
        'infinity_active': tuple(sorted(kappa * digit + 2 * e
                                        for digit, e in enumerate(exponents))),
        'endpoint_repair': 'Constant diagonal physical-Y inverse; no degree or profile shift.',
    }
    if include_spectators:
        result.update({
            'zero_spectators': tuple(ramification * (2 * i + 1 - charge)
                                     for i in range(n) if i != selected)
            + tuple(ramification * (2 * i + 1 + charge)
                    for i in range(n) if i != retained),
            'infinity_spectators': tuple(ramification * (2 * i + 1 - 2 * n)
                                         for i in range(2 * n) if i != selected),
            'paired_physical_first': tuple(ramification * (4 * i - 2 * (n - 1) + charge)
                                            for i in range(n)),
            'paired_physical_second': tuple(ramification * (4 * i - 2 * (n - 1) - charge)
                                             for i in range(n)),
        })
    return result


def scalar_filters(exponents, coefficients):
    return {(digit,): {-exponent: coefficient}
            for digit, (exponent, coefficient) in enumerate(zip(exponents, coefficients))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=3)
    parser.add_argument('--n', type=int, default=9)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--charge', type=Fraction, default=Fraction(0))
    parser.add_argument('--geometry', choices=('initial', 'terminal'), default='initial')
    parser.add_argument('--exponents', default='0,1,1')
    parser.add_argument('--coefficients', default='1,2,1')
    args = parser.parse_args()
    data = profile_data(args.prime, args.n, args.ramification, args.charge,
                        tuple(map(int, args.exponents.split(','))),
                        tuple(map(int, args.coefficients.split(','))), args.geometry)
    print(json.dumps(data, indent=2, default=str))


if __name__ == '__main__':
    main()
