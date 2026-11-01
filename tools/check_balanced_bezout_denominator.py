#!/usr/bin/env python3
"""Exact parity divided differences and balanced Bezout denominators."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import random
from check_interval_bezout_precision import multiply, scaled_lagrange_inverse, normalized_hensel_inverse

ROOT = Path(__file__).resolve().parents[1]


def order(value):
    value = Fraction(value)
    if not value:
        return None
    return ((value.numerator & -value.numerator).bit_length()
            - (value.denominator & -value.denominator).bit_length())


def factorial_order(n):
    return n - n.bit_count()


def isometry(nodes):
    return all(order(nodes[i] - nodes[j]) == order(i - j)
               for i in range(len(nodes)) for j in range(i))


def parity(value):
    value = Fraction(value)
    assert value.denominator % 2
    return value.numerator % 2


def divided_differences(nodes, values=None):
    values = [Fraction(parity(x)) for x in nodes] if values is None else list(map(Fraction, values))
    result = [values[0]]
    for degree in range(1, len(nodes)):
        values = [(values[i + 1] - values[i]) / (nodes[i + degree] - nodes[i])
                  for i in range(len(values) - 1)]
        result.append(values[0])
    return result


def newton_polynomial(nodes, coefficients):
    result, basis = [Fraction(0)] * len(nodes), [Fraction(1)]
    for degree, coefficient in enumerate(coefficients):
        for i, value in enumerate(basis):
            result[i] += coefficient * value
        if degree + 1 < len(nodes):
            basis = multiply(basis, [-nodes[degree], Fraction(1)])
    return result


def evaluate(polynomial, x):
    result = Fraction(0)
    for coefficient in reversed(polynomial):
        result = result * x + coefficient
    return result


def rescale_polynomial(polynomial, h):
    result = [Fraction(0)] * len(polynomial)
    for degree, coefficient in enumerate(polynomial):
        for i in range(degree + 1):
            result[i] += coefficient * Fraction(comb(degree, i) * (-1) ** (degree - i), 1 << (h * degree))
    return result


def integer_isometry(length, rng):
    bits = max(1, (length - 1).bit_length())
    mask = rng.randrange(1 << bits)
    values = [(i ^ mask) + 2 * (i ^ mask) ** 2 + (rng.randrange(8) << bits)
              for i in range(length)]
    assert isometry(values)
    return values


def controls():
    rng = random.Random(2026110201)
    dd_checks = interpolation_checks = norm_checks = bezout_checks = moment_checks = unit_checks = 0
    examples = []
    for length in range(2, 25):
        raw = integer_isometry(length, rng)
        candidates = [list(map(Fraction, range(length))), [Fraction(3 * x + 2, 2 * x + 1) for x in raw]]
        for nodes in candidates:
            assert isometry(nodes)
            dd = divided_differences(nodes)
            for k in range(1, length):
                assert order(dd[k]) == k - 1 - factorial_order(k)
                dd_checks += 1
            polynomial = newton_polynomial(nodes, dd)
            for node in nodes:
                assert evaluate(polynomial, node) == parity(node)
                interpolation_checks += 1
            for h in (1, 2, 4):
                expected = (h - 1) * (length - 1) + 1 + factorial_order(length - 1)
                actual = max(-order(x) for x in rescale_polynomial(polynomial, h) if x)
                assert actual == expected
                norm_checks += 1
            if length <= 12:
                k = length - 1
                complete = [Fraction(1)] + [Fraction(0)] * 5
                for node in nodes:
                    for j in range(1, 6):
                        complete[j] += node * complete[j - 1]
                normalized = (-1) ** (k + 1) * factorial(k) * dd[k] / (1 << (k - 1))
                partial = Fraction(0)
                for j in range(6):
                    degree = k + j
                    elementary = [1] + [0] * j
                    for value in range(degree):
                        for t in range(j, 0, -1):
                            elementary[t] += value * elementary[t - 1]
                    numerator = sum((-1) ** (j - t) * elementary[j - t] * complete[t]
                                    for t in range(j + 1))
                    direct_values = [prod(node - t for t in range(degree)) / factorial(degree)
                                     for node in nodes]
                    direct = divided_differences(nodes, direct_values)[-1]
                    assert direct == numerator / factorial(degree)
                    quotient = numerator / prod(range(k + 1, k + j + 1))
                    assert not quotient or order(quotient) >= 0
                    partial += (-1) ** j * (1 << j) * quotient
                    error = normalized - partial
                    assert not error or order(error) >= j + 1
                    moment_checks += 1
                    unit_checks += 1
    for m in list(range(1, 26)) + [33]:
        for irregular in (False, True):
            aa0 = integer_isometry(m, rng) if irregular else list(range(m))
            bb0 = integer_isometry(m, rng) if irregular else list(range(m))
            aa, bb = [1 + 4 * x for x in aa0], [1 + 4 * x for x in bb0]
            n = 3 + 4 * rng.randrange(20)
            nodes = [node for a, b in zip(aa, bb)
                     for node in (Fraction(a - 1, 2), (Fraction(n, b) - 1) / 2)]
            assert isometry(nodes)
            expected = m + factorial_order(m - 1)
            first = scaled_lagrange_inverse(n, aa, bb)
            second = normalized_hensel_inverse(n, aa, bb)
            assert first['dyadic_denominator'] == second['dyadic_denominator'] == expected
            dd = divided_differences(nodes)
            assert all(order(dd[k]) == k - 1 - factorial_order(k) for k in range(1, 2 * m))
            dd_checks += 2 * m - 1
            bezout_checks += 1
            if not irregular:
                examples.append({'m': m, 'n': n, 'rows': aa, 'columns': bb,
                                 'dyadic_denominator': expected})
    bad_a, bad_b = [1, 5, 9, 13], [1, 5, 9, 17]
    bad_nodes = [node for a, b in zip(bad_a, bad_b)
                 for node in (Fraction(a - 1, 2), (Fraction(3, b) - 1) / 2)]
    assert not isometry(bad_nodes)
    assert scaled_lagrange_inverse(3, bad_a, bad_b)['dyadic_denominator'] == 6
    m = (1 << 200) + 13
    implicit = {'m': m, 'n': 3, 'a_progression': [1, 4, m], 'b_progression': [1, 4, m],
                'dyadic_denominator': m + factorial_order(m - 1),
                'equivalent_formula': 1 + factorial_order(2 * m - 1),
                'ordinary_resultant_valuation': m * m, 'enumerated_nodes': 0,
                'scope': 'Closed-form theorem application, not expanded interpolation.'}
    assert implicit['dyadic_denominator'] == implicit['equivalent_formula']
    return {'exact_divided_difference_checks': dd_checks, 'exact_interpolation_values': interpolation_checks,
            'affine_coefficient_norm_checks': norm_checks, 'balanced_bezout_inverse_checks': bezout_checks,
            'low_surplus_symmetric_identities': moment_checks, 'normalized_unit_truncation_checks': unit_checks,
            'examples': examples, 'implicit_progression_family': implicit,
            'isometry_assumption_counterexample': {'rows': bad_a, 'columns': bad_b, 'n': 3,
                                                 'actual': 6, 'inapplicable_formula': 5}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'claim': 'Exact parity divided differences on dyadic isometries and balanced Bezout denominator',
              'controls': controls(), 'factoring_exponent_improved': False,
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'check_interval_bezout_precision.py')}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result['controls'].items() if isinstance(value, int)}, indent=2))


if __name__ == '__main__':
    main()
