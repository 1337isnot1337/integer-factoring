#!/usr/bin/env python3
"""Exact polynomial controls for arbitrary nonresonant Boolean patterns."""

import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random


def evaluate_anf(coefficients, state):
    return sum((coefficients >> monomial) & 1 for monomial in range(state + 1)
               if monomial & state == monomial) % 2


def multiply(first, second):
    result = 0
    while first:
        low = first & -first
        result ^= second << (low.bit_length() - 1)
        first ^= low
    return result


def full_transform(values, modes, m):
    powers = [sum(1 << (m * shift) for shift in range(power + 1) if comb(power, shift) % 2)
              for power in range(modes + 1)]
    result = []
    for state in range(1 << modes):
        value = 0
        for source in range(1 << modes):
            if source & state == source:
                value ^= multiply(powers[state.bit_count() - source.bit_count()], values[source])
        assert value and value.bit_count() % 2 == 1
        result.append(value)
    return result


def orders(values):
    return tuple((value & -value).bit_length() - 1 for value in values)


def rational_loss_control(n, P, modes, degree, L, truth, beta_orders, base, charge):
    m = n * P
    physical = [base ** (P * sum((2 * n - 1 - 2 * bit) for bit in range(modes)
                                if (state >> bit) & 1)) for state in range(1 << modes)]
    shifted = [weight * base ** (P * charge * state.bit_count())
               for state, weight in enumerate(physical)]
    p = sum(weight for weight, active in zip(physical, truth) if active) / sum(physical)
    zero = sum(weight * base ** (2 * (order - m * state.bit_count()))
               for state, (weight, order) in enumerate(zip(shifted, beta_orders))) / sum(shifted)
    infinity = sum(weight * base ** (-2 * L * active)
                   for weight, active in zip(physical, truth)) / sum(physical)
    ratio = zero * infinity
    bound = p * base ** (-2 * m * degree)
    assert ratio > bound > (1 + base ** (2 * m)) ** (-n)
    if base ** (-2 * m) * (1 - base ** P) >= 1:
        assert ratio > 1


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root / 'evidence' /
                        'nonresonant-boolean-filter-floor.json')
    args = parser.parse_args()
    generator = random.Random(721091)
    patterns = []
    polynomial_cases = state_controls = near_power_controls = rational_controls = 0
    for modes in (1, 2, 3, 4):
        masks = (list(range(2, 1 << (1 << modes))) if modes <= 3 else
                 sorted({generator.randrange(2, 1 << (1 << modes)) for _ in range(128)}))
        n = 1
        while n < modes:
            n *= 2
        patterns.append({'modes': modes, 'patterns': len(masks),
                         'exhaustive': modes <= 3})
        for pattern_index, coefficients in enumerate(masks):
            degree = max(monomial.bit_count() for monomial in range(1 << modes)
                         if (coefficients >> monomial) & 1)
            top = sum(1 << monomial for monomial in range(1 << modes)
                      if monomial.bit_count() == degree and (coefficients >> monomial) & 1)
            truth = tuple(evaluate_anf(coefficients, state) for state in range(1 << modes))
            for P in (1, 3):
                m = n * P
                gaps = {1, m - 1, m + 1, m * degree - 1, m * degree + 1, m * (degree + 1)}
                for gap in sorted(value for value in gaps if value > 0):
                    L = m * degree + gap
                    if L % m == 0 and L <= 2 * m * degree:
                        continue
                    values = tuple(1 ^ ((1 << (m * degree)) ^ (1 << L) if active else 0)
                                   for active in truth)
                    beta = full_transform(values, modes, m)
                    beta_orders = orders(beta)
                    for state, order in enumerate(beta_orders):
                        natural = m * state.bit_count()
                        assert natural <= order <= natural + gap
                        assert values[state].bit_length() - 1 == L * truth[state]
                        assert values[state].bit_count() % 2 == 1
                        state_controls += 1
                        if gap < m:
                            assert order == natural + gap * evaluate_anf(top, state)
                            near_power_controls += 1
                    polynomial_cases += 1
                    if pattern_index % 31 == 0:
                        for base in (Fraction(4, 5), Fraction(49, 50)):
                            for charge in (-2 * n, 0, 2 * n):
                                rational_loss_control(n, P, modes, degree, L, truth,
                                                      beta_orders, base, charge)
                                rational_controls += 1
    negative_values = (1, 1 ^ (1 << 4) ^ (1 << 6), 1, 1)
    negative_beta = full_transform(negative_values, 2, 2)
    assert orders(negative_beta) == (0, 2, 2, 8)
    assert negative_beta[3] == 1 << 8
    assert orders(negative_beta)[3] > 2 * 2 + 6 - 2 * 2
    result = {
        'status': 'Every nonresonant exact order, full degree, endpoint, and rational loss control passed.',
        'scope': 'Arbitrary Boolean g in the stated two-power family outside the finite resonant window; no factoring exponent improvement.',
        'pattern_families': patterns, 'literal_polynomial_cases': polynomial_cases,
        'complete_state_order_degree_endpoint_controls': state_controls,
        'exact_near_power_top_component_orders': near_power_controls,
        'rational_loss_and_large_scale_controls': rational_controls,
        'resonant_negative_control': {
            'n': 2, 'P': 1, 'm': 2, 'degree': 2, 'L': 6,
            'g_anf_monomial_masks': [1, 3], 'full_degrees': [0, 6, 0, 0],
            'beta_orders': list(orders(negative_beta)),
            'state_11_beta': 'x^8', 'false_order_upper_bound': 6,
            'scope': 'The generic upper order cap fails; no negative full pressure is claimed.'
        },
        'proof_dependencies': ['homogeneous-boolean-filter-bounds.md: biased support lemma and threshold',
                               'global-homogeneous-boolean-filter-floor.md: Gaussian support floor'],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
