#!/usr/bin/env python3
"""Exact rational controls for the charged one-mode Laurent envelopes."""
import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

from check_one_mode_binary_laurent_filters import endpoint, inspect, polynomial
from charged_gaussian_floor import Dyadic


def probabilities(n, P, j, J):
    # t=P*log(2), so every probability and eta=4 are rational.
    lam = 2 * n - 1 - 2 * j
    p1 = 1 / (1 + Fraction(2) ** ((lam + J) * P))
    q1 = 1 / (1 + Fraction(2) ** (lam * P))
    return 1 - p1, p1, 1 - q1, q1


def envelopes(n, P, j, J, exponent):
    p0, p1, q0, q1 = probabilities(n, P, j, J)
    H, z = Fraction(4) ** (n * P), Fraction(4) ** exponent
    A = (p0 + p1 * z) * (q0 / z + q1)
    B = (p0 + p1 * z) * (q0 + q1 * H / z)
    return H, z, A, B, q0


def literal_example(n, P, j, J, first, second, expected):
    orders, degrees, branch, scalar = inspect(first, second, n * P)
    p0, p1, q0, q1 = probabilities(n, P, j, J)
    eta = Fraction(4)
    value = (p0 * eta ** -orders[0] + p1 * eta ** (n * P - orders[1]))
    value *= q0 * eta ** degrees[0] + q1 * eta ** degrees[1]
    assert value == expected and q0 < value < 1
    return {'n': n, 'P': P, 'j': j, 'charge': J, 'scale': 'P*log(2)',
            'f0': sorted(first.items()), 'f1': sorted(second.items()),
            'beta_orders': orders, 'full_degrees': degrees,
            'ratio': str(value), 'gaussian_comparison_only': True}


def exterior_kernel_cover():
    """Prove 2f(t)+f(3t)-2sqrt(2)f(sqrt(2)t)>0 on [1/32,5/2]."""
    arithmetic = Dyadic(80)
    pending = [('', Fraction(1, 32), Fraction(5, 2))]
    cover = []
    evaluated = 0
    while pending:
        path, left, right = pending.pop()
        scale = (arithmetic.rational(left.numerator, left.denominator)[0],
                 arithmetic.rational(right.numerator, right.denominator)[1])
        positive = arithmetic.add(
            arithmetic.multiply(arithmetic.rational(2), arithmetic.softplus(scale)),
            arithmetic.softplus(arithmetic.multiply(arithmetic.rational(3), scale)))
        negative = arithmetic.multiply(
            arithmetic.multiply(arithmetic.rational(2), arithmetic.root_two),
            arithmetic.softplus(arithmetic.multiply(arithmetic.root_two, scale)))
        lower = arithmetic.subtract(positive, negative)[0]
        evaluated += 1
        if lower > 0:
            cover.append([path, str(lower)])
        else:
            assert len(path) < 20
            middle = (left + right) / 2
            pending.extend(((path + '1', middle, right), (path + '0', left, middle)))
    assert Fraction(2, 21) - Fraction(5, 64) > 0
    assert 1 + 1 + Fraction(1, 2) + Fraction(1, 6) > Fraction(5, 2)
    return {'status': 'Every interval has a positive integer outward lower bound.',
            'domain': ['1/32', '5/2'], 'dyadic_precision': 80,
            'intervals_evaluated': evaluated, 'positive_cover_intervals': len(cover),
            'minimum_lower_numerator': str(min(int(value) for path, value in cover)),
            'common_denominator': str(arithmetic.scale), 'cover': cover}


def coefficient_barrier_controls():
    pairs = barrier_cases = 0
    for field, width, shift in ((2, 7, -3), (4, 4, -2)):
        catalogue = tuple(value for coefficients in product(range(field), repeat=width)
                          if (value := polynomial(coefficients, shift)) and endpoint(value))
        for m in (1, 2, 4):
            for first in catalogue:
                for second in catalogue:
                    orders, degrees, branch, scalar = inspect(first, second, m)
                    if orders[1] > orders[0] + m and orders[1] > degrees[1]:
                        assert max(degrees) >= orders[0] + m
                        barrier_cases += 1
                    pairs += 1
    assert barrier_cases > 0
    return {'complete_polynomial_pairs': pairs,
            'genuine_high_cancellation_barrier_cases': barrier_cases}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'charged-one-mode-filters.json')
    args = parser.parse_args()
    band_count = cap_count = 0
    for n in (1, 2, 4, 8, 32):
        for P in (1, 3):
            for j in sorted({0, n // 2, n - 1}):
                lam = 2 * n - 1 - 2 * j
                low, high = -2 * n, 2 * n
                charges = sorted({low, low + 1, -1, 0, 1, n, high - 1, high})
                exponents = sorted({-2 * n * P, -n * P, -2, -1, 0, 1, 2,
                                    n * P - 1, n * P, n * P + 1, 2 * n * P})
                for J in charges:
                    assert low <= J <= high
                    for exponent in exponents:
                        H, z, A, B, q0 = envelopes(n, P, j, J, exponent)
                        if z < 1:
                            assert B > 1
                            p0, p1, q0, q1 = probabilities(n, P, j, J)
                            if J >= -2 * lam:
                                assert A > 1
                            else:
                                first_degree = (p0 + p1 * z) * (q0 * max(H, 1 / z) + q1)
                                second_degree = (p0 + p1 * z) * (q0 / z + q1 * H)
                                assert first_degree > 1 and second_degree > 1
                        elif z > H:
                            assert B > 1
                        elif z == H:
                            assert B > 1
                        else:
                            p0, p1, _, _ = probabilities(n, P, j, J)
                            assert p0 + p1 * z >= 1
                        band_count += 1
                for J in (-100 * n, low - 1, high + 1, 100 * n):
                    for exponent in exponents:
                        H, z, A, B, q0 = envelopes(n, P, j, J, exponent)
                        if z <= 1:
                            assert A > q0
                        assert B > q0
                        cap_count += 1
    examples = [
        literal_example(2, 1, 0, 8, {0: 1}, {-1: 1}, Fraction(1936, 2049)),
        literal_example(2, 1, 1, -12, {0: 1, 2: 1, 4: 1}, {0: 1}, Fraction(513, 683)),
    ]
    record = {'status': 'All exact physical-probability envelope controls pass.',
              'strict_charge_band': '-2*n <= J <= 2*n, uniformly over the selected index',
              'band_controls': band_count, 'outside_charge_loss_cap_controls': cap_count,
              'complete_laurent_relative_improvement_examples': examples,
              'coefficient_at_xm_barrier': coefficient_barrier_controls(),
              'strong_binary_exterior_kernel': exterior_kernel_cover(),
              'scope': 'The uniform theorem is proved analytically; these finite controls use exact rational arithmetic and include both charge boundaries.'}
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    summary = {key: value for key, value in record.items()
               if key not in ('complete_laurent_relative_improvement_examples',
                              'strong_binary_exterior_kernel')}
    summary['strong_binary_exterior_kernel'] = {
        key: value for key, value in record['strong_binary_exterior_kernel'].items()
        if key != 'cover'}
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
