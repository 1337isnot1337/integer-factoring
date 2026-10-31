#!/usr/bin/env python3
"""Exact all-size profiles for prime-characteristic identity-block filters.

The proof is in research/prime-identity-filter-diagonalization.md.
At zero, the p-by-p active state arrays are regularly equivalent to
diagonal finite differences of the p input digit filters. Computing a
profile therefore requires no truncated series or matrix elimination.

All exponent collisions, including collisions between different powers
of u=z^(2*n*P), are combined in the finite field before taking orders.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb, isqrt


def validate(prime: int, n: int, ramification: int, exponents):
    if not isinstance(prime, int) or prime < 2 or any(
            prime % d == 0 for d in range(2, isqrt(prime) + 1)):
        raise ValueError('The characteristic must be prime.')
    if not isinstance(n, int) or n < 1:
        raise ValueError('The interval length must be a positive integer.')
    if not isinstance(ramification, int) or ramification < 1 or not ramification % 2:
        raise ValueError('Use a positive odd ramification for the crossing endpoint.')
    if len(exponents) != prime or any(not isinstance(e, int) for e in exponents):
        raise ValueError('Supply one integer filter exponent for every digit.')


def beta_polynomials(prime: int, n: int, ramification: int, exponents):
    """Return exact Laurent polynomials beta_k after monomial substitution.

    alpha_d=z^(-2*e_d), C=1-z^(2*n*P), and
    beta_k=sum_j (-1)^(k-j) binom(k,j) C^j alpha_(p-1-k+j).
    """
    validate(prime, n, ramification, exponents)
    result = []
    for k in range(prime):
        terms = {}
        for j in range(k + 1):
            e = exponents[prime - 1 - k + j]
            outer = (-1)**(k - j) * comb(k, j)
            for h in range(j + 1):
                degree = 2 * n * ramification * h - 2 * e
                coefficient = outer * (-1)**h * comb(j, h)
                terms[degree] = (terms.get(degree, 0) + coefficient) % prime
        terms = {degree: value for degree, value in terms.items() if value}
        if not terms:
            raise ArithmeticError('A finite difference vanished despite its nonzero endpoint.')
        result.append(terms)
    return tuple(result)


def active_profiles(prime: int, n: int, ramification: int = 1,
                    charge: int | Fraction = 1, exponents=None,
                    include_spectators: bool = True):
    """Return exact active orders and optional one-particle spectators.

    The filters are alpha_d=(-v)^(-e_d), after w=v^P. Multiplication
    by a common scalar monomial shifts the two cut orders oppositely;
    no such clearing shift is included in this return value.

    `charge` is J in the zero physical gauge X += P*J, Y -= P*J.
    Each spectator order expands to the p digit orders 0,s,...,(p-1)s.
    """
    if exponents is None:
        exponents = (0,) * prime
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an integer or Fraction for the exact physical gauge.')
    differences = beta_polynomials(prime, n, ramification, exponents)
    orders = tuple(min(value) for value in differences)
    a = ramification * (1 - charge)
    kappa = ramification * (1 - 2 * n)
    zero = tuple(sorted(orders[k] + a * (digit - k)
                        for k in range(prime) for digit in range(prime)))
    result = {
        'prime': prime,
        'n': n,
        'ramification': ramification,
        'charge': charge,
        'exponents': tuple(exponents),
        'beta_orders': orders,
        'zero_active': zero,
        'zero_active_first_factor': tuple(a * digit for digit in range(prime)),
        'zero_active_second_factor': tuple(orders[k] - a * k for k in range(prime)),
        'infinity_active': tuple(sorted(kappa * digit + 2 * e
                                        for digit, e in enumerate(exponents))),
    }
    if include_spectators:
        result.update({
            'zero_spectators': tuple(ramification * (2 * i + 1 - charge)
                                     for i in range(1, n))
            + tuple(ramification * (2 * i + 1 + charge)
                    for i in range(n - 1)),
            'infinity_spectators': tuple(ramification * (2 * i + 1 - 2 * n)
                                         for i in range(1, 2 * n)),
            'paired_physical_first': tuple(ramification * (4 * i - 2 * (n - 1) + charge)
                                            for i in range(n)),
            'paired_physical_second': tuple(ramification * (4 * i - 2 * (n - 1) - charge)
                                             for i in range(n)),
        })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=5)
    parser.add_argument('--n', type=int, default=16)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--charge', type=Fraction, default=Fraction(1))
    parser.add_argument('--exponents', default='0,1,2,2,1')
    parser.add_argument('--omit-spectators', action='store_true')
    args = parser.parse_args()
    exponents = tuple(int(value) for value in args.exponents.split(','))
    result = active_profiles(args.prime, args.n, args.ramification,
                             args.charge, exponents,
                             include_spectators=not args.omit_spectators)
    print(json.dumps(result, indent=2, default=str))


if __name__ == '__main__':
    main()
