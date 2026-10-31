#!/usr/bin/env python3
"""Exact fixed-size profiles for identity-block qutrit state filters.

The universal zero tensor has only charge blocks of orders 1, 2, 3, 2, 1.
Its multivariate minors are computed once over F3[u,A,B]. Substitution of
monomial parameters groups equal exponents before taking a valuation;
there is no series cutoff and no generic-coefficient assumption.

Here u=z^(2*n*P), A=z^(-2*e1), B=z^(-2*e2). The charge-block shift is
P*(1-J)*charge. The infinity active orders are 0, P*(1-2*n)+2*e1,
and 2*P*(1-2*n)+2*e2. The all-size endpoint and regular-basis proof is
in research/truncated-symmetric-identity-filters.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import factorial
from typing import Mapping

Polynomial = dict[tuple[int, int, int], int]
Algebra = dict[tuple[int, int, int, int], Polynomial]
ONE = {(0, 0, 0): 1}


def polynomial_add(first: Mapping, second: Mapping, coefficient: int = 1):
    result = dict(first)
    for powers, value in second.items():
        result[powers] = (result.get(powers, 0) + coefficient * value) % 3
        if not result[powers]:
            del result[powers]
    return result


def polynomial_multiply(first: Mapping, second: Mapping):
    result = {}
    for powers, a in first.items():
        for other, b in second.items():
            target = tuple(i + j for i, j in zip(powers, other))
            result[target] = (result.get(target, 0) + a * b) % 3
    return {powers: value for powers, value in result.items() if value}


def algebra_multiply(first: Mapping, second: Mapping):
    result = {}
    for powers, a in first.items():
        for other, b in second.items():
            target = tuple(i + j for i, j in zip(powers, other))
            if max(target) >= 3:
                continue
            result[target] = polynomial_add(
                result.get(target, {}), polynomial_multiply(a, b))
    return {powers: value for powers, value in result.items() if value}


def polynomial_determinant(matrix):
    result = {}
    for permutation in permutations(range(len(matrix))):
        term = ONE
        for i, j in enumerate(permutation):
            term = polynomial_multiply(term, matrix[i][j])
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(len(matrix))
                         for j in range(i + 1, len(matrix)))
        result = polynomial_add(result, term, -1 if inversions % 2 else 1)
    return result


@lru_cache(maxsize=1)
def universal_blocks():
    """Return exact F3[u,A,B] matrices in canonical charge order."""
    gaussian = {
        (s, r, r, s): {(s, 0, 0): pow(factorial(r) * factorial(s), -1, 3)}
        for r, s in product(range(3), repeat=2)
    }
    q = {
        (0, 1, 1, 0): ONE,
        (0, 1, 0, 1): {(0, 0, 0): 2},
        (1, 0, 1, 0): {(0, 0, 0): 1, (1, 0, 0): 2},
        (1, 0, 0, 1): {(0, 0, 0): 2, (1, 0, 0): 1},
    }
    qg = algebra_multiply(q, gaussian)
    qqg = algebra_multiply(q, qg)
    tensor = dict(gaussian)
    linear = {(0, 1, 0): 1, (0, 0, 0): 2}
    quadratic = {(0, 0, 1): 2, (0, 1, 0): 2, (0, 0, 0): 2}
    for coefficient, summands in ((linear, qg), (quadratic, qqg)):
        for powers, value in summands.items():
            tensor[powers] = polynomial_add(
                tensor.get(powers, {}), polynomial_multiply(coefficient, value))

    blocks = []
    for charge in range(-2, 3):
        rows = [(x, y) for x, y in product(range(3), repeat=2) if y - x == charge]
        columns = [(v, u) for v, u in product(range(3), repeat=2) if u - v == charge]
        matrix = []
        for x, y in rows:
            row = []
            for v, u in columns:
                factor = factorial(u) * factorial(x) % 3
                row.append({powers: value * factor % 3
                            for powers, value in tensor.get((v, y, u, x), {}).items()
                            if value * factor % 3})
            matrix.append(tuple(row))
        blocks.append((charge, tuple(rows), tuple(columns), tuple(matrix)))
    return tuple(blocks)


@lru_cache(maxsize=1)
def universal_minors():
    result = []
    for charge, rows, columns, matrix in universal_blocks():
        levels = []
        for size in range(1, len(matrix) + 1):
            level = []
            for chosen_rows in combinations(range(len(matrix)), size):
                for chosen_columns in combinations(range(len(matrix)), size):
                    value = polynomial_determinant(
                        [[matrix[i][j] for j in chosen_columns] for i in chosen_rows])
                    if value:
                        level.append(value)
            levels.append(tuple(level))
        result.append((charge, tuple(levels)))
    return tuple(result)


def monomial_specialization(value, weights):
    """Aggregate exact cancellations after assigning monomial weights."""
    result = {}
    for powers, coefficient in value.items():
        exponent = sum(a * b for a, b in zip(powers, weights))
        result[exponent] = (result.get(exponent, 0) + coefficient) % 3
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def _validate(n, ramification, e1, e2):
    if not isinstance(n, int) or n < 1:
        raise ValueError('The interval length must be a positive integer.')
    if not isinstance(ramification, int) or ramification < 1 or not ramification % 2:
        raise ValueError('The characteristic-three endpoint requires a positive odd ramification.')
    if not isinstance(e1, int) or not isinstance(e2, int):
        raise ValueError('Supply integer exponents for the two digit filters.')


def scalar_valuations(n: int, ramification: int, e1: int, e2: int):
    """Orders of S=alpha2, D=(1-u)alpha2-alpha1, and
    F=1+(1-u)alpha1+(1-u)^2 alpha2 for monomial digit filters.
    """
    _validate(n, ramification, e1, e2)
    size = 2 * n * ramification
    s = -2 * e2
    d = size - 2 * e1 if e1 == e2 else -2 * max(e1, e2)
    f = 2 * size if e1 == e2 == 0 else -2 * max(0, e1, e2)
    return s, d, f


def active_orders_from_minors(n: int, ramification: int,
                             charge: int | Fraction, e1: int, e2: int):
    """Independent fixed-minor route, retaining exact cancellations."""
    _validate(n, ramification, e1, e2)
    a = ramification * (1 - charge)
    weights = (2 * n * ramification, -2 * e1, -2 * e2)
    zero = []
    for block_charge, levels in universal_minors():
        previous = 0
        for level in levels:
            valuations = []
            for value in level:
                specialized = monomial_specialization(value, weights)
                if specialized:
                    valuations.append(min(specialized))
            if not valuations:
                raise ArithmeticError('A promised nonzero charge-block minor vanished.')
            current = min(valuations)
            zero.append(current - previous + a * block_charge)
            previous = current
    return tuple(sorted(zero))


def active_profiles(n: int, ramification: int = 1,
                    charge: int | Fraction = 1, e1: int = 0, e2: int = 1,
                    include_spectators: bool = True):
    """Return the complete exact zero and infinity active profiles.

    Zero spectators have orders P*(2*i+1-J), i=1,...,n-1, and
    P*(2*i+1+J), i=0,...,n-2. Infinity spectators have orders
    P*(2*i+1-2*n), i=1,...,2*n-1. Each spectator is a qutrit mode.
    """
    _validate(n, ramification, e1, e2)
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an integer or Fraction for an exact physical gauge.')
    p = ramification
    a = p * (1 - charge)
    s, d, f = scalar_valuations(n, p, e1, e2)
    zero = (s + 2 * a, s + a, d + a, s, d, f,
            d - a, f - a, f - 2 * a)
    kappa = p * (1 - 2 * n)
    result = {
        'zero_active': tuple(sorted(zero)),
        'infinity_active': tuple(sorted((0, kappa + 2 * e1, 2 * kappa + 2 * e2))),
        'scalar_valuations': (s, d, f),
    }
    if not include_spectators:
        return result
    result.update({
        'zero_spectators': tuple(p * (2 * i + 1 - charge) for i in range(1, n))
        + tuple(p * (2 * i + 1 + charge) for i in range(n - 1)),
        'infinity_spectators': tuple(p * (2 * i + 1 - 2 * n)
                                    for i in range(1, 2 * n)),
        'paired_physical_first': tuple(p * (4 * i - 2 * (n - 1) + charge)
                                       for i in range(n)),
        'paired_physical_second': tuple(p * (4 * i - 2 * (n - 1) - charge)
                                        for i in range(n)),
    })
    return result


def polynomial_seed(n: int, ramification: int = 1, e1: int = 0, e2: int = 1):
    """Construct the actual finite polynomial tensor for a small interval.

    Coefficients use the raw truncated exponential convention, in variable
    order (V,Y,U,X). Multiply coefficients by the U and X occupation
    factorials for the physical normalization. The returned clearing
    order shifts every zero profile order by +2*q and every infinity
    order by -2*q after v=-z^2; its scalar endpoint is already normalized.
    Full construction is exponential in n and is a finite reference path,
    not the fixed-size profile algorithm.
    """
    _validate(n, ramification, e1, e2)
    from check_all_field_pascal_crossings import (
        add, multiply, pascal_kernel, scale, trim,
    )
    from check_truncated_symmetric_crossings import direct_kernel

    def kernel_product(first, second):
        result = {}
        for powers, x in first.items():
            for other, y in second.items():
                target = tuple(a + b for a, b in zip(powers, other))
                if max(target) >= 3:
                    continue
                result[target] = add(result.get(target, ()), multiply(x, y, 3), 3)
        return {powers: value for powers, value in result.items() if value}

    gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, 3), 3)
    left = {}
    powers = [0] * (4 * n)
    powers[n] = 1
    left[tuple(powers)] = (1,)
    powers = [0] * (4 * n)
    powers[0] = 1
    left[tuple(powers)] = (1, 1)
    right = {}
    powers = [0] * (4 * n)
    powers[2 * n] = 1
    right[tuple(powers)] = (1,)
    for b in range(n):
        powers = [0] * (4 * n)
        powers[3 * n + b] = 1
        right[tuple(powers)] = (0,) * (2 * n - 1 - b) + ((-1) ** b % 3,)
    qg = kernel_product(kernel_product(left, right), gaussian)
    qqg = kernel_product(kernel_product(left, right), qg)
    clearing = max(0, e1, e2)
    parity_sign = lambda exponent: -1 if exponent % 2 else 1
    sign = parity_sign(clearing)

    def scalar_terms(terms):
        result = {}
        for exponent, value in terms:
            result[exponent] = (result.get(exponent, 0) + value) % 3
        return {exponent: value for exponent, value in result.items() if value}

    coefficients = (
        {0: 1},
        scalar_terms(((-e1, parity_sign(e1)), (0, -1))),
        scalar_terms(((-e2, 2 * parity_sign(e2)),
                      (-e1, -parity_sign(e1)), (0, 2))),
    )
    result = {}
    for kernel, scalar in zip((gaussian, qg, qqg), coefficients):
        for powers, polynomial in kernel.items():
            if not scalar:
                continue
            degree = clearing + max(scalar) + ramification * (len(polynomial) - 1)
            value = [0] * (degree + 1)
            for exponent, coefficient in scalar.items():
                for i, x in enumerate(polynomial):
                    at = clearing + exponent + ramification * i
                    value[at] = (value[at] + sign * coefficient * x) % 3
            result[powers] = add(result.get(powers, ()), trim(value, 3), 3)
    return {powers: value for powers, value in result.items() if value}, clearing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=16)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--charge', type=Fraction, default=Fraction(1))
    parser.add_argument('--e1', type=int, default=0)
    parser.add_argument('--e2', type=int, default=1)
    args = parser.parse_args()
    result = active_profiles(args.n, args.ramification, args.charge, args.e1, args.e2)
    print(json.dumps(result, indent=2, default=str))


if __name__ == '__main__':
    main()
