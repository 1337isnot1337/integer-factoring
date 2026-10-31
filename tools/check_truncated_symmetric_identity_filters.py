#!/usr/bin/env python3
"""Independent exact controls for the qutrit identity-filter family."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from math import factorial
from pathlib import Path

from check_all_field_pascal_crossings import (
    pascal_kernel, ramify_negative, scale, weighted_orders,
)
from check_truncated_symmetric_crossings import basis, direct_kernel
from truncated_symmetric_identity_filters import (
    ONE, active_orders_from_minors, active_profiles, polynomial_add as add,
    polynomial_determinant as determinant, polynomial_multiply as multiply,
    polynomial_seed, universal_blocks,
)


def matrix_multiply(left, right):
    return [[sum_polynomials(multiply(left[i][k], right[k][j])
                             for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def sum_polynomials(values):
    result = {}
    for value in values:
        result = add(result, value)
    return result


def diagonal(values):
    return [[value if i == j else {} for j in range(len(values))]
            for i, value in enumerate(values)]


def block_factorizations():
    u = {(1, 0, 0): 1}
    r = {(0, 1, 0): 1}
    s = {(0, 0, 1): 1}
    c = add(ONE, u, 2)
    c2 = multiply(c, c)
    negative = lambda value: add({}, value, 2)
    d = add(multiply(c, s), r, 2)
    f = sum_polynomials((ONE, multiply(c, r), multiply(c2, s)))
    plus_left = [[negative(ONE), c], [negative(ONE), add(c, ONE, 2)]]
    plus_right = [[ONE, negative(c)], [ONE, add(ONE, c, 2)]]
    minus = [[ONE, negative(c)], [ONE, add(ONE, c, 2)]]
    middle_left = [
        [ONE, c, c2],
        [ONE, add(c, ONE), add(c2, c, 2)],
        [negative(ONE), add(ONE, c, 2), negative(sum_polynomials((c2, c, ONE)))],
    ]
    middle_right = [
        [ONE, negative(c), negative(c2)],
        [ONE, negative(add(c, ONE)), add(c, c2, 2)],
        [ONE, add(ONE, c, 2), negative(sum_polynomials((c2, c, ONE)))],
    ]
    assert determinant(plus_left) == determinant(plus_right) == ONE
    assert determinant(minus) == ONE
    assert determinant(middle_left) == negative(ONE)
    assert determinant(middle_right) == ONE
    independent = {
        2: [[add({}, s, 2)]],
        1: matrix_multiply(matrix_multiply(plus_left, diagonal((d, s))), plus_right),
        0: matrix_multiply(matrix_multiply(middle_left, diagonal((f, d, s))),
                           middle_right),
        -1: matrix_multiply(matrix_multiply(minus, diagonal((f, d))), minus),
        -2: [[add({}, f, 2)]],
    }
    expected_determinants = {
        2: add({}, s, 2),
        1: multiply(s, d),
        0: negative(multiply(multiply(s, d), f)),
        -1: multiply(d, f),
        -2: add({}, f, 2),
    }
    entries = 0
    for charge, rows, columns, matrix in universal_blocks():
        raw = independent[charge]
        assert determinant(raw) == expected_determinants[charge]
        for i, (x, y) in enumerate(rows):
            for j, (v, occupation_u) in enumerate(columns):
                normalized = add({}, raw[i][j], factorial(occupation_u) * factorial(x))
                assert normalized == matrix[i][j]
                entries += 1
    return {'independent_polynomial_entry_identities': entries,
            'unit_determinant_factorizations': 5,
            'exact_charge_block_determinants': 5}


def finite_difference_profiles():
    controls = 0
    for n, p, charge, e1, e2 in product(
            (1, 2, 5), (1, 3), (Fraction(-2), Fraction(0), Fraction(1), Fraction(5, 2)),
            range(-3, 4), range(-3, 4)):
        expected = active_orders_from_minors(n, p, charge, e1, e2)
        actual = active_profiles(n, p, charge, e1, e2, include_spectators=False)
        assert actual['zero_active'] == expected
        controls += 1
    huge = []
    for n, p, e1, e2 in (
            (2 ** 209, 3 ** 41, 0, 0),
            (2 ** 209, 3 ** 41, 7, 7),
            (3 ** 128, 3 ** 31, -3 ** 159, 3 ** 159),
            (10 ** 100, 1, -10 ** 100, -10 ** 100)):
        actual = active_profiles(n, p, Fraction(7, 3), e1, e2, include_spectators=False)
        assert actual['zero_active'] == active_orders_from_minors(
            n, p, Fraction(7, 3), e1, e2)
        huge.append({'n_bits': n.bit_length(), 'ramification_bits': p.bit_length(),
                     'zero_active': list(map(str, actual['zero_active']))})
    return {'closed_form_against_all_minors': controls,
            'large_integer_controls': huge}


def expand_qutrit_modes(orders):
    result = [0]
    for order in orders:
        result = [value + digit * order for value in result for digit in range(3)]
    return result


def direct_profiles(kernel, n, p, charge):
    digits = basis(n, 3)
    lookup = {word: i for i, word in enumerate(digits)}
    dimension = 3 ** n
    tensor = {}
    for powers, coefficient in kernel.items():
        labels = tuple(lookup[powers[i * n:(i + 1) * n]] for i in range(4))
        factor = 1
        for digit in powers[2 * n:]:
            factor *= factorial(digit)
        tensor[labels] = ramify_negative(scale(coefficient, factor, 3), 3)
    heights = (
        [p * charge] * n,
        [p * (2 * a + 1) for a in range(n)],
        [-p * charge] * n,
        [p * (2 * a - 4 * n + 1) for a in range(n)],
        [p * (2 * a - n + 2) for a in range(n)],
        [-p * n] * n,
        [p * (2 * a - n + 1) for a in range(n)],
        [p * (3 * n - 1)] * n,
    )
    answer = []
    precisions = []
    for endpoint, order in enumerate(((3, 1, 0, 2), (3, 2, 0, 1))):
        full = [[()] * (dimension * dimension) for _ in range(dimension * dimension)]
        for labels, value in tensor.items():
            i = labels[order[0]] * dimension + labels[order[1]]
            j = labels[order[2]] * dimension + labels[order[3]]
            full[i][j] = value
        rows, columns, charges = [], [], []
        for a, b in product(range(dimension), repeat=2):
            charges.append(sum(digits[a]) + (1 if endpoint else -1) * sum(digits[b]))
            rows.append(sum(heights[4 * endpoint + order[k]][i] * digits[letter][i]
                            for k, letter in enumerate((a, b)) for i in range(n)))
            columns.append(sum(heights[4 * endpoint + order[k + 2]][i] * digits[letter][i]
                               for k, letter in enumerate((a, b)) for i in range(n)))
        orders = []
        for sector in sorted(set(charges)):
            indices = [i for i, value in enumerate(charges) if value == sector]
            profile, precision = weighted_orders(
                [[full[i][j] for j in indices] for i in indices],
                [rows[i] for i in indices], [columns[j] for j in indices],
                3, infinity=bool(endpoint))
            orders.extend(profile)
            precisions.append(precision)
        answer.append(tuple(sorted(orders)))
    return answer, max(precisions)


def full_tensor_controls():
    cases = []
    endpoint_entries = 0
    choices = [(1, 1, e1, e2, j)
               for e1, e2 in ((0, 0), (0, 1), (1, 1), (1, 2), (-1, -2), (1, -1))
               for j in (0, 1)]
    choices += [(2, 1, e1, e2, j)
                for e1, e2 in ((0, 0), (0, 1), (1, 1), (-1, 2)) for j in (0, 1)]
    choices += [(2, 3, 0, 1, 1), (2, 3, -1, 2, 1)]
    cache = {}
    for n, p, e1, e2, charge in choices:
        key = (n, p, e1, e2)
        if key not in cache:
            cache[key] = polynomial_seed(n, p, e1, e2)
        kernel, clearing = cache[key]
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, 3), 3)
        evaluate = lambda value: sum(coefficient * (-1) ** i
                                     for i, coefficient in enumerate(value)) % 3
        for powers in product(range(3), repeat=4 * n):
            assert evaluate(kernel.get(powers, ())) == evaluate(gaussian.get(powers, ()))
            endpoint_entries += 1
        actual, precision = direct_profiles(kernel, n, p, charge)
        data = active_profiles(n, p, charge, e1, e2)
        zero_spectators = expand_qutrit_modes(data['zero_spectators'])
        infinity_spectators = expand_qutrit_modes(data['infinity_spectators'])
        expected = [
            tuple(sorted(a + b + 2 * clearing for a in data['zero_active']
                         for b in zero_spectators)),
            tuple(sorted(a + b - 2 * clearing for a in data['infinity_active']
                         for b in infinity_spectators)),
        ]
        assert actual == expected
        cases.append({'n': n, 'P': p, 'J': charge, 'e1': e1, 'e2': e2,
                      'clearing_order': clearing, 'matrix_order': 3 ** (2 * n),
                      'maximum_series_precision': precision})
    return {'endpoint_coefficient_equalities': endpoint_entries,
            'complete_two_cut_controls': cases}


def domination_controls():
    count = 0
    for n, p, evaluation, e1, e2 in product(
            range(1, 5), (1, 3), (Fraction(11, 10), Fraction(2)),
            range(-2, 3), range(-2, 3)):
        for charge in range(1 - 2 * n, 2 * n + 1):
            data = active_profiles(n, p, charge, e1, e2, include_spectators=False)
            actual = (sum(evaluation ** -order for order in data['zero_active'])
                      * sum(evaluation ** -order for order in data['infinity_active']))
            a = p * (1 - charge)
            b = p * (2 * n - 1 + charge)
            kappa = p * (1 - 2 * n)
            background = (sum(evaluation ** (-a * j) for j in range(3))
                          * sum(evaluation ** (-b * j) for j in range(3))
                          * sum(evaluation ** (-kappa * j) for j in range(3)))
            assert actual >= background
            assert (actual == background) == (e1 == e2 == 0)
            count += 1
    return {'exact_rational_comparisons': count,
            'charge_hypothesis': '1-2*n <= J <= 2*n',
            'comparison': 'Gaussian background at the same charge and scale'}


def failed_positive_mode_controls():
    from check_all_field_pascal_crossings import (
        add as polynomial_add, binomial_power, multiply as polynomial_multiply,
    )
    records = []
    for n in (2, 3):
        c_power = binomial_power(n, 3)
        geometric = ()
        for j in range(n):
            geometric = polynomial_add(geometric, binomial_power(j, 3), 3)
        schur_numerator = polynomial_add(
            (0,) * (2 * n - 1) + c_power,
            (0,) * (2 * n) + geometric, 3, -1)
        assert schur_numerator == (0,) * (2 * n - 1) + (1,)
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, 3), 3)
        powers = [0] * (4 * n)
        powers[0] = powers[2 * n] = 2
        assert gaussian[tuple(powers)] == scale(binomial_power(2, 3), 2, 3)
        # The retained one-copy diagonal coefficient has reduced added
        # numerator w^(4n-3)/2 and denominator (1+w)^(2n-3).
        numerator = (0,) * (4 * n - 3) + (2,)
        at_endpoint = sum(value * (-1) ** i for i, value in enumerate(numerator)) % 3
        assert at_endpoint != 0 and 2 * n - 3 > 0
        weighted_order = 2 * (4 * n - 3) + 2 * (1 + 1 - 4 * n)
        assert weighted_order == -2
        records.append({'n': n, 'one_copy_endpoint_pole_order': 2 * n - 3,
                        'zero_gauged_added_entry_order': weighted_order})
    return {'exact_schur_and_coefficient_controls': records,
            'scope': 'Unrepaired full-alphabet positive-one Smith summand'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'pass', 'block_factorizations': block_factorizations(),
              'profile_controls': finite_difference_profiles(),
              'full_tensor_controls': full_tensor_controls(),
              'domination_controls': domination_controls(),
              'failed_positive_mode_controls': failed_positive_mode_controls(),
              'new_factoring_exponent_claimed': False}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
