#!/usr/bin/env python3
"""Exact controls for slow endpoint corrections and affine Smith profiles.

Small bivariate determinants are expanded directly. Their predicted
weighted profiles are compared with the separate power-series Smith
implementation. No floating computation is used in a certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from math import comb, factorial
from pathlib import Path
import random

from check_all_field_pascal_crossings import (
    add, multiply, pascal_kernel, ramify_negative, scale, trim,
    weighted_orders,
)
from check_truncated_symmetric_crossings import direct_kernel


def bivariate_add(a, b, prime, sign=1):
    result = dict(a)
    for exponent, value in b.items():
        result[exponent] = (result.get(exponent, 0) + sign * value) % prime
        if not result[exponent]:
            del result[exponent]
    return result


def bivariate_multiply(a, b, prime):
    result = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = (i + k, j + l)
            result[key] = (result.get(key, 0) + x * y) % prime
    return {key: value for key, value in result.items() if value}


def bivariate_determinant(matrix, prime):
    size = len(matrix)
    answer = {}
    for permutation in permutations(range(size)):
        term = {(0, 0): 1}
        for i, j in enumerate(permutation):
            term = bivariate_multiply(term, matrix[i][j], prime)
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(size) for j in range(i + 1, size))
        answer = bivariate_add(answer, term, prime,
                               -1 if inversions % 2 else 1)
    return answer


def substitute(value, ramification, prime):
    """Substitute w=v^P, b=1+v, independently of the band formulas."""
    if not value:
        return ()
    degree = max(ramification * a + b for a, b in value)
    result = [0] * (degree + 1)
    for (a, b), coefficient in value.items():
        for j in range(b + 1):
            result[ramification * a + j] = (
                result[ramification * a + j] + coefficient * comb(b, j)
            ) % prime
    return trim(result, prime)


def minor_bands(value, prime, infinity):
    edge = (max if infinity else min)(a for a, _ in value)
    coefficients = {b: coefficient for (a, b), coefficient in value.items()
                    if a == edge}
    if infinity:
        return -edge, -max(coefficients)
    shifted = [0] * (max(coefficients) + 1)
    for b, coefficient in coefficients.items():
        for j in range(b + 1):
            shifted[j] = (shifted[j] + coefficient * comb(b, j)) % prime
    return edge, next(i for i, value in enumerate(shifted) if value)


def all_minors(matrix, prime):
    size = len(matrix)
    result = []
    for k in range(1, size + 1):
        level = []
        for rows in combinations(range(size), k):
            for columns in combinations(range(size), k):
                value = bivariate_determinant(
                    [[matrix[i][j] for j in columns] for i in rows], prime)
                if value:
                    level.append((rows, columns, value))
        result.append(level)
    return result


def generic_prefixes(minors, row_weights, column_weights, prime, infinity):
    answer = [(Fraction(0), 0)]
    for level in minors:
        candidates = []
        for rows, columns, value in level:
            leading, secondary = minor_bands(value, prime, infinity)
            leading += sum(row_weights[i] for i in rows)
            leading += sum(column_weights[j] for j in columns)
            candidates.append((leading, secondary))
        answer.append(min(candidates))
    return answer


def ramify_polynomial(value, denominator):
    if not value:
        return ()
    result = [0] * (denominator * (len(value) - 1) + 1)
    for i, coefficient in enumerate(value):
        result[denominator * i] = coefficient
    return tuple(result)


def matrix_controls():
    rng = random.Random(1337)
    cases = minor_equalities = complete_profiles = 0
    maximum_precision = 0
    for prime, size, b_degree, denominator in product(
            (2, 3, 5), range(1, 5), (1, 2), (1, 2)):
        while True:
            matrix = [[{
                (a, b): coefficient
                for a, b in product(range(3), range(b_degree + 1))
                if (coefficient := rng.randrange(prime))
            } for _ in range(size)] for _ in range(size)]
            minors = all_minors(matrix, prime)
            if minors[-1]:
                break
        rows = [Fraction(rng.randrange(-2, 3), denominator)
                for _ in range(size)]
        columns = [Fraction(rng.randrange(-2, 3), denominator)
                   for _ in range(size)]
        threshold = denominator * size * b_degree
        prefixes = [generic_prefixes(minors, rows, columns, prime, infinity)
                    for infinity in (False, True)]
        for ramification in (threshold + 1, threshold + 3):
            specialized = [[substitute(value, ramification, prime)
                            for value in row] for row in matrix]
            for infinity, prefix in zip((False, True), prefixes):
                for level in minors:
                    for _, _, value in level:
                        polynomial = substitute(value, ramification, prime)
                        actual = (1 - len(polynomial) if infinity else
                                  next(i for i, x in enumerate(polynomial) if x))
                        leading, secondary = minor_bands(value, prime, infinity)
                        assert actual == ramification * leading + secondary
                        minor_equalities += 1
                expected = tuple(
                    ramification * (prefix[k][0] - prefix[k - 1][0])
                    + prefix[k][1] - prefix[k - 1][1]
                    for k in range(1, size + 1))
                integer_matrix = [[ramify_polynomial(value, denominator)
                                   for value in row] for row in specialized]
                integer_rows = [int(ramification * denominator * x) for x in rows]
                integer_columns = [int(ramification * denominator * x)
                                   for x in columns]
                actual, precision = weighted_orders(
                    integer_matrix, integer_rows, integer_columns, prime,
                    infinity=infinity)
                assert tuple(Fraction(x, denominator) for x in actual) == expected
                maximum_precision = max(maximum_precision, precision)
                complete_profiles += 1
        cases += 1

    # At the non-strict degree threshold, substitution can annihilate a
    # nonzero generic minor: w-(b-1)^2 becomes zero for w=v^2,b=1+v.
    bad_minor = {(1, 0): 1, (0, 2): 2, (0, 1): 2, (0, 0): 2}
    assert bad_minor and not substitute(bad_minor, 2, 3)
    return {
        'seed': 1337,
        'random_matrices': cases,
        'separated_minor_equalities': minor_equalities,
        'complete_weighted_profile_equalities': complete_profiles,
        'maximum_series_precision': maximum_precision,
        'strict_threshold_negative_control': True,
    }


def qutrit_tensor(ramification):
    prime = 3
    original = direct_kernel(pascal_kernel((0,), 2, prime), prime)
    kernel = {powers: ramify_polynomial(value, ramification)
              for powers, value in original.items()}
    answer = {powers: scale((0,) + value, -1, prime)
              for powers, value in kernel.items()}
    monomial = (2, 0, 2, 0)
    for powers, value in kernel.items():
        target = tuple(a + b for a, b in zip(powers, monomial))
        if max(target) >= prime:
            continue
        correction = scale(multiply(value, (1, 1), prime), -1, prime)
        answer[target] = add(answer.get(target, ()), correction, prime)
    return original, {powers: value for powers, value in answer.items() if value}


def qutrit_profiles(kernel, ramification):
    tensor = {}
    for labels, value in kernel.items():
        multiplier = factorial(labels[2]) * factorial(labels[3])
        tensor[labels] = ramify_negative(scale(value, multiplier, 3), 3)
    heights = (0, 1, 0, -3, 1, -1, 0, 2)
    profiles = []
    precisions = []
    for endpoint, order in enumerate(((3, 1, 0, 2), (3, 2, 0, 1))):
        matrix = [[()] * 9 for _ in range(9)]
        for labels, value in tensor.items():
            row = 3 * labels[order[0]] + labels[order[1]]
            column = 3 * labels[order[2]] + labels[order[3]]
            matrix[row][column] = value
        rows = [ramification * (heights[4 * endpoint + order[0]] * a
                                + heights[4 * endpoint + order[1]] * b)
                for a, b in product(range(3), repeat=2)]
        columns = [ramification * (heights[4 * endpoint + order[2]] * a
                                   + heights[4 * endpoint + order[3]] * b)
                   for a, b in product(range(3), repeat=2)]
        profile, precision = weighted_orders(
            matrix, rows, columns, 3, infinity=bool(endpoint))
        profiles.append(profile)
        precisions.append(precision)
    return profiles, precisions


def qutrit_controls():
    records = []
    endpoint_equalities = 0
    for ramification in range(1, 22, 2):
        original, kernel = qutrit_tensor(ramification)
        for powers in product(range(3), repeat=4):
            old = sum(coefficient * (-1) ** k for k, coefficient
                      in enumerate(original.get(powers, ()))) % 3
            new = sum(coefficient * (-1) ** k for k, coefficient
                      in enumerate(kernel.get(powers, ()))) % 3
            assert old == new
            endpoint_equalities += 1
        profiles, precisions = qutrit_profiles(kernel, ramification)
        p = ramification
        expected_zero = (0, p + 2, p + 2, 2 * p + 2, 2 * p + 2,
                         2 * p + 2, 3 * p + 2, 3 * p + 2, 4 * p + 4)
        expected_infinity = tuple(sorted((-2 * p - 2, -p - 2, -p - 2,
                                          -2, -2, 0, p - 2, p, 2 * p)))
        assert profiles == [expected_zero, expected_infinity]
        records.append({'P': p, 'zero': profiles[0], 'infinity': profiles[1],
                        'series_precisions': precisions})

    # Rational cost certificate; the common physical energy is zero.
    p = 3
    evaluation = Fraction(11, 10)
    x = evaluation ** -p
    r = evaluation ** 2
    letter = 1 + x + x * x
    gaussian = letter ** 4 / (x * x)
    zero = r + 2 * x + 3 * x * x + 2 * x ** 3 + x ** 4 / r
    infinity = x ** -2 + 2 / x + 2 + 1 / r + x + x / r + x * x / r
    ratio = zero * infinity / gaussian
    expected_ratio = Fraction(99612598272277477054122841,
                              101096444552703336598632100)
    assert ratio == expected_ratio and ratio < 1
    profile = records[1]
    direct = (sum(evaluation ** -order for order in profile['zero'])
              * sum(evaluation ** -order for order in profile['infinity']))
    assert direct == zero * infinity

    zero_leading = (0, 1, 1, 2, 2, 2, 3, 3, 4)
    zero_secondary = (0, 2, 2, 2, 2, 2, 2, 2, 4)
    infinity_leading = (-2, -1, -1, 0, 0, 0, 1, 1, 2)
    infinity_secondary = (-2, -2, -2, -2, -2, 0, -2, 0, 0)
    for x in (Fraction(1, 3), Fraction(1, 2), Fraction(3, 4), Fraction(9, 10)):
        weighted = sum(b * x ** a for a, b in zip(zero_leading, zero_secondary))
        weighted /= sum(x ** a for a in zero_leading)
        other = sum(b * x ** a for a, b in
                    zip(infinity_leading, infinity_secondary))
        other /= sum(x ** a for a in infinity_leading)
        formula = 2 * (-1 + x * x + x ** 3 + 2 * x ** 4) / (1 + x + x * x) ** 2
        assert weighted + other == formula

    # Integer polynomial identity proving gamma(x*)=-1/2 when
    # 3*x*^2+x*-1=0. Coefficient lists are in increasing degree.
    first = (-1, 1, 3)
    second = (3, 1, 3)
    multiplied = [0] * 5
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            multiplied[i + j] += a * b
    assert multiplied == [-3, 2, 7, 6, 9]
    numerator_at_zero = 3
    numerator_derivative_at_zero = -13
    letter_at_zero = 3
    letter_derivative_at_zero = -3
    gamma_at_zero = Fraction(2 * numerator_at_zero, letter_at_zero ** 2)
    gamma_derivative_at_zero = Fraction(
        2 * (numerator_derivative_at_zero * letter_at_zero
             - 2 * numerator_at_zero * letter_derivative_at_zero),
        letter_at_zero ** 3)
    assert gamma_at_zero == Fraction(2, 3)
    assert gamma_derivative_at_zero == Fraction(-14, 9)

    return {
        'endpoint_coefficient_equalities': endpoint_equalities,
        'affine_threshold': 18,
        'profiles': records,
        'beyond_threshold_interpolation_values': [19, 21],
        'cost_certificate': {'P': 3, 'evaluation': '11/10',
                             'ratio_numerator': ratio.numerator,
                             'ratio_denominator': ratio.denominator,
                             'strictly_less_than_one': True},
        'secondary_expectation_identity_controls': 4,
        'gaussian_optimum_secondary_expectation': '-1/2',
        'small_scale_secondary_expectation': ['2/3', '-14/9'],
        'factoring_exponent_improvement_claimed': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'pass', 'matrix_controls': matrix_controls(),
              'qutrit_controls': qutrit_controls()}
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
