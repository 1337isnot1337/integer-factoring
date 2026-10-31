#!/usr/bin/env python3
"""All-index joint identity filters at prime-power Pascal interval lengths.

A reflected full Pascal change simultaneously separates every Q_j.
Consequently arbitrary joint diagonal filters have an exact tensor-product
finite-difference profile, independent of the chosen indices until the
explicit physical heights are added.  No Smith calculation is required.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb

from coupled_identity_extreme_filters import is_power
from prime_state_filter_profiles import validate


def validate_selection(prime, n, ramification, selected):
    validate(prime, n, ramification, (0,) * prime)
    if not is_power(n, prime):
        raise ValueError('The simultaneous reflected-Pascal theorem requires n=p^k.')
    selected = tuple(selected)
    if (not selected or len(set(selected)) != len(selected)
            or any(not isinstance(j, int) or not 0 <= j < n for j in selected)):
        raise ValueError('Select distinct identity indices between zero and n−1.')
    return selected


def field_add(first, second, prime):
    return first ^ second if prime == 2 else (first + second) % prime


def prime_scale(value, scalar, prime):
    return value if prime == 2 and scalar & 1 else (0 if prime == 2 else value * scalar % prime)


def normalized_alpha(prime, modes, alpha, field_size=None):
    field_size = prime if field_size is None else field_size
    if field_size != prime and not (prime == 2 and isinstance(field_size, int)
                                   and field_size > 1 and field_size & (field_size - 1) == 0):
        raise ValueError('Use Fp, or a fixed binary extension represented in an additive basis.')
    result = {}
    for state in product(range(prime), repeat=modes):
        polynomial = alpha.get(state, {})
        value = {degree: coefficient for degree, coefficient in polynomial.items() if coefficient}
        if any(not isinstance(degree, int) or not 0 <= coefficient < field_size
               for degree, coefficient in value.items()):
            raise ValueError('Use complete integer Laurent exponents and valid field coefficients.')
        endpoint = 0
        for coefficient in value.values():
            endpoint = field_add(endpoint, coefficient, prime)
        if not endpoint:
            raise ValueError('Every physical-Y endpoint weight must be nonzero.')
        result[state] = value
    return result


@lru_cache(maxsize=16)
def finite_difference_terms(prime, modes):
    """For each k, return (source state, power of C, prime-field scalar)."""
    result = {}
    for k in product(range(prime), repeat=modes):
        terms = []
        for j in product(*(range(value + 1) for value in k)):
            scalar = (-1) ** (sum(k) - sum(j))
            for a, b in zip(k, j):
                scalar *= comb(a, b)
            scalar %= prime
            if scalar:
                state = tuple(prime - 1 - a + b for a, b in zip(k, j))
                terms.append((state, sum(j), scalar))
        result[k] = tuple(terms)
    return result


def beta_polynomials(prime, n, ramification, alpha, field_size=None):
    modes = len(next(iter(alpha)))
    validate_selection(prime, n, ramification, tuple(range(modes)))
    alpha = normalized_alpha(prime, modes, alpha, field_size)
    m = n * ramification
    c_powers = tuple(tuple((m * h, comb(k, h) * (-1) ** h % prime)
                           for h in range(k + 1) if comb(k, h) % prime)
                     for k in range(modes * (prime - 1) + 1))
    result = {}
    for k, terms in finite_difference_terms(prime, modes).items():
        polynomial = {}
        for state, power, scalar in terms:
            for shift, coefficient in c_powers[power]:
                coefficient = scalar * coefficient % prime
                for degree, value in alpha[state].items():
                    exponent = shift + degree
                    polynomial[exponent] = field_add(polynomial.get(exponent, 0),
                                                     prime_scale(value, coefficient, prime), prime)
        polynomial = {degree: value for degree, value in polynomial.items() if value}
        assert polynomial  # beta_k(1)=(-1)^|k| alpha_(p-1-k)(1) is nonzero.
        result[k] = polynomial
    return result


def beta_orders(prime, n, ramification, alpha, field_size=None):
    return tuple(min(value) for value in beta_polynomials(prime, n, ramification, alpha, field_size).values())


def zero_blocks(prime, n, ramification, selected, orders):
    selected = validate_selection(prime, n, ramification, selected)
    states = tuple(product(range(prime), repeat=len(selected)))
    if len(orders) != len(states):
        raise ValueError('Give one beta valuation for every joint digit state.')
    blocks = defaultdict(list)
    for k, value in zip(states, orders):
        for auxiliary in states:
            charge = sum(auxiliary) - sum(k)
            order = 2 * value + 2 * ramification * sum(
                j * (a - b) for j, a, b in zip(selected, auxiliary, k))
            blocks[charge].append(order)
    return tuple({'charge': charge, 'orders': tuple(sorted(values))}
                 for charge, values in sorted(blocks.items()))


def profile_data(prime, n, ramification, charge, selected, alpha, field_size=None):
    selected = validate_selection(prime, n, ramification, selected)
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an integer or Fraction for exact physical charges.')
    modes = len(selected)
    alpha = normalized_alpha(prime, modes, alpha, field_size)
    betas = beta_polynomials(prime, n, ramification, alpha, field_size)
    orders = tuple(min(value) for value in betas.values())
    blocks = zero_blocks(prime, n, ramification, selected, orders)
    retained = tuple(n - 1 - j for j in selected)
    return {
        'prime': prime, 'n': n, 'ramification': ramification, 'modes': modes,
        'coefficient_field': field_size or prime, 'charge': charge,
        'identity_indices': selected, 'retained_xv_indices': retained,
        'beta_orders': orders,
        'zero_active': tuple(sorted(value + ramification * (1 - charge) * block['charge']
                                    for block in blocks for value in block['orders'])),
        'zero_charge_blocks_at_J1': blocks,
        'infinity_active': tuple(sorted(ramification * sum((2 * j + 1 - 2 * n) * digit
                                                           for j, digit in zip(selected, state))
                                        - 2 * max(polynomial)
                                        for state, polynomial in alpha.items())),
        'zero_spectators': tuple(ramification * (2 * j + 1 - charge)
                                 for j in range(n) if j not in selected)
        + tuple(ramification * (2 * j + 1 + charge)
                for j in range(n) if j not in retained),
        'infinity_spectators': tuple(ramification * (2 * j + 1 - 2 * n)
                                     for j in range(2 * n) if j not in selected),
        'paired_physical_first': tuple(ramification * (4 * j - 2 * (n - 1) + charge)
                                       for j in range(n)),
        'paired_physical_second': tuple(ramification * (4 * j - 2 * (n - 1) - charge)
                                        for j in range(n)),
    }
