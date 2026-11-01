#!/usr/bin/env python3
"""Certify an all-scale pressure gain for a 32-letter coherent crossing.

The two full native profiles are recomputed for both displayed maps. The
comparison uses an integer polynomial identity. A separate outward-integer
interval calculation gives a rational strict rate bound. No floating point
is used by this checker, and no post-half factoring bound is asserted.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from check_large_interval_pressure import BinaryIntervals
from coherent_hot_compounds import height, profiles


ROOT = Path(__file__).resolve().parents[1]
BOUNDARY = (24, 36, 66, 129)


def displayed_maps(coherent):
    lower = sorted((w for w in range(256) if height(w) < 16),
                   key=lambda w: (height(w), w))
    assert len(lower) == 28
    base = [1 << word for word in lower + [15]]
    right, left = ([3, 5, 9], [13, 11, 8]) if coherent else ([1, 2, 4], [1, 2, 4])

    def lift(value):
        return sum(1 << word for index, word in enumerate(BOUNDARY)
                   if value >> index & 1)

    return base + list(map(lift, right)), base + list(map(lift, left))


def counting_polynomial(orders):
    count = Counter(orders)
    return [count[index] for index in range(max(orders) + 1)]


def multiply(first, second):
    answer = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            answer[i + j] += a * b
    return answer


def comparison(coordinate, coherent):
    assert coordinate['heights'] == coherent['heights']
    c0 = counting_polynomial(coordinate['zero_orders_t'])
    ci = counting_polynomial(coordinate['infinity_orders_t'])
    z = counting_polynomial(coherent['zero_orders_t'])
    i = counting_polynomial(coherent['infinity_orders_t'])
    left, right = multiply(c0, ci), multiply(z, i)
    difference = [(left[j] if j < len(left) else 0)
                  - (right[j] if j < len(right) else 0)
                  for j in range(max(len(left), len(right)))]
    while difference[-1] == 0:
        difference.pop()
    assert difference[:32] == [0] * 32
    assert difference[32] == 1
    running, quotient = 0, []
    for value in difference[32:-1]:
        running += value
        quotient.append(running)
    assert running + difference[-1] == 0
    assert [0] * 32 + multiply([1, -1], quotient) == difference
    assert len(quotient) == 148 and quotient[0] == 1
    assert all(value >= 0 for value in quotient)
    return {'coordinate_zero_polynomial': c0, 'coordinate_infinity_polynomial': ci,
            'coherent_zero_polynomial': z, 'coherent_infinity_polynomial': i,
            'pressure_difference_coefficients': difference,
            'factorization': 'r^32 * (1-r) * Q(r)',
            'Q_coefficients': quotient, 'Q_degree': len(quotient) - 1,
            'Q_coefficient_sum': sum(quotient),
            'Q_nonnegative_coefficients': True, 'Q_constant_coefficient': 1,
            'strict_pressure_gain_interval': '0 < r < 1',
            'physical_height_multisets_identical': True}


def polynomial_interval(arithmetic, polynomial, value):
    answer = None
    for coefficient in reversed(polynomial):
        if answer is not None:
            answer = arithmetic.multiply(answer, value)
        if coefficient:
            term = arithmetic.rational(coefficient, 1)
            answer = term if answer is None else arithmetic.add(answer, term)
    assert answer is not None
    return answer


def ratio_interval(record, numerator, denominator, root_numerator, root_denominator):
    arithmetic = BinaryIntervals(256)
    root = arithmetic.rational(root_numerator, root_denominator)
    r = arithmetic.power(root, numerator)
    physical = arithmetic.power(root, denominator)
    z = polynomial_interval(arithmetic, counting_polynomial(record['zero_orders_t']), r)
    i = polynomial_interval(arithmetic, counting_polynomial(record['infinity_orders_t']), r)
    a = polynomial_interval(arithmetic, counting_polynomial(record['heights']), physical)
    left = arithmetic.power(arithmetic.multiply(z, i), denominator)
    right = arithmetic.power(a, 4 * numerator)
    return arithmetic.multiply(left, arithmetic.inverse(right))


def interval_controls():
    arithmetic = BinaryIntervals(128)
    checks = 0
    for polynomial in ([1, 0, 2, 1], [0, 3, 0, 0, 1], [7], [1, 3, 5, 7, 9]):
        for root in (Fraction(2, 3), Fraction(7, 11), Fraction(13, 12)):
            interval = polynomial_interval(arithmetic, polynomial,
                                           arithmetic.rational(root.numerator, root.denominator))
            lower, upper = arithmetic.fraction_bounds(interval)
            exact = sum(coefficient * root ** power for power, coefficient in enumerate(polynomial))
            assert lower <= exact <= upper
            checks += 1
    record = {'zero_orders_t': [0, 1, 4], 'infinity_orders_t': [0, 2, 3],
              'heights': [0, 1, 3]}
    for numerator, denominator in ((2, 3), (5, 7)):
        root = Fraction(7, 9)
        z = sum(root ** (numerator * value) for value in record['zero_orders_t'])
        i = sum(root ** (numerator * value) for value in record['infinity_orders_t'])
        a = sum(root ** (denominator * value) for value in record['heights'])
        exact = (z * i) ** denominator / a ** (4 * numerator)
        lower, upper = BinaryIntervals.fraction_bounds(
            ratio_interval(record, numerator, denominator, 7, 9))
        assert lower <= exact <= upper
        checks += 1
    return checks


def encode_interval(interval):
    return {'lower_mantissa': str(interval[0]), 'upper_mantissa': str(interval[1]),
            'binary_exponent': interval[2]}


def rate_certificate(coordinate, coherent):
    # With r=(5405/5406)^1781 and rho=1781/2501, every Gibbs weight
    # is rational: r^(h/rho)=(5405/5406)^(2501*h).
    numerator, denominator = 1781, 2501
    coherent_interval = ratio_interval(coherent, numerator, denominator, 5405, 5406)
    coordinate_interval = ratio_interval(coordinate, numerator, denominator, 5405, 5406)
    lower, upper = BinaryIntervals.fraction_bounds(coherent_interval)
    coordinate_lower, _ = BinaryIntervals.fraction_bounds(coordinate_interval)
    assert 0 < lower <= upper < 1 < coordinate_lower
    gap = 1 - upper
    rho = Fraction(numerator, denominator)
    return {'rho_strict_upper_bound': str(rho),
            'counting_variable': 'r=(5405/5406)^1781',
            'physical_type_formula': 'p_j=(5405/5406)^(2501*h_j) / sum_k (5405/5406)^(2501*h_k)',
            'certified_inequality': '(Z0(r)*Zinfinity(r))^2501 < A((5405/5406)^2501)^7124',
            'coherent_ratio_interval': encode_interval(coherent_interval),
            'coordinate_ratio_interval': encode_interval(coordinate_interval),
            'strict_ratio_gap_lower_bound': [str(gap.numerator), str(gap.denominator)],
            'coordinate_fails_this_rate_at_this_scale': True,
            'mantissa_precision': 256,
            'post_half_target_reached': rho * rho < Fraction(1, 2)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = [profiles(*displayed_maps(coherent), 8) for coherent in (False, True)]
    coordinate, coherent = records
    report = {'status': 'pass', 'scope': __doc__,
              'complete_profiles_recomputed': 4,
              'complete_smith_orders': sum(len(record[key]) for record in records
                                           for key in ('zero_orders_t', 'infinity_orders_t')),
              'independent_interval_controls': interval_controls(),
              'physical_maps': [{'embedding': r['embedding'], 'retraction': r['retraction'],
                                  'heights': r['heights'], 'charges': r['charges'],
                                  'zero_orders_t': r['zero_orders_t'],
                                  'infinity_orders_t': r['infinity_orders_t']} for r in records],
              'all_scale_pressure_comparison': comparison(coordinate, coherent),
              'rational_rate_certificate': rate_certificate(coordinate, coherent),
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'coherent_hot_compounds.py',
                                             'coherent_hot_compounds.cpp',
                                             'selected_gaussian_coefficients.cpp',
                                             'native_rectangular_smith_f4.cpp',
                                             'check_large_interval_pressure.py')},
              'floating_arithmetic_used_for_certification': False,
              'global_optimality_claimed': False,
              'factoring_exponent_improved': False}
    encoded = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('physical_maps', 'all_scale_pressure_comparison', 'source_sha256')},
                     indent=2))


if __name__ == '__main__':
    main()
