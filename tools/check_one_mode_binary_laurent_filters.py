#!/usr/bin/env python3
"""Exact finite controls for the all-degree one-mode binary filter theorem."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path


def polynomial(coefficients, shift):
    return {i + shift: value for i, value in enumerate(coefficients) if value}


def endpoint(value):
    result = 0
    for coefficient in value.values():
        result ^= coefficient
    return result


def transform(first, second, m):
    result = dict(second)
    for exponent, coefficient in first.items():
        for position in (exponent, exponent + m):
            result[position] = result.get(position, 0) ^ coefficient
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def inspect(first, second, m):
    beta = transform(first, second, m)
    assert beta and endpoint(beta) == endpoint(second)
    b0, b1 = min(first), min(beta)
    d0, d1 = max(first), max(second)
    assert d0 >= b0
    if b1 <= d1:
        branch = 'diagonal'
        assert m - b1 + d1 >= m > 0
    else:
        branch = 'crossed'
        assert min(second) == b0 and b1 <= d0 + m
        assert d1 >= b0 and d0 + m >= b1
    equality_data = b1 > d1 and b1 == b0 + m and d0 == d1 == b0
    scalar = first == second and len(first) == 1
    assert equality_data == scalar
    return (b0, b1), (d0, d1), branch, scalar


def ratio(orders, degrees, m, q, p0):
    p1 = 1 - p0
    u0, u1 = q ** -orders[0], q ** (m - orders[1])
    v0, v1 = (q ** degree for degree in degrees)
    return (p0 * u0 + p1 * u1) * (p0 * v0 + p1 * v1)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'one-mode-binary-laurent-filters.json')
    args = parser.parse_args()
    rows = []
    examples = {}
    rational_checks = 0
    for field, width, shift in ((2, 7, -3), (4, 4, -2)):
        catalogue = tuple(value for coefficients in product(range(field), repeat=width)
                          if (value := polynomial(coefficients, shift)) and endpoint(value))
        counts = Counter()
        profile_classes = set()
        rational_profiles = set()
        for m in (1, 2, 7, 32):
            for first in catalogue:
                for second in catalogue:
                    orders, degrees, branch, scalar = inspect(first, second, m)
                    counts[branch] += 1
                    counts['scalar_equalities'] += scalar
                    key = (orders, degrees, m)
                    profile_classes.add(key)
                    if key not in rational_profiles:
                        rational_profiles.add(key)
                        for q in (Fraction(3, 2), Fraction(2)):
                            for p0 in (Fraction(3, 5), Fraction(2, 3), Fraction(4, 5)):
                                cost = ratio(orders, degrees, m, q, p0)
                                assert cost >= 1 and (cost == 1) == scalar
                                rational_checks += 1
                    label = (field, branch, scalar)
                    if label not in examples:
                        examples[label] = {
                            'field_order': field, 'm': m,
                            'f0': sorted(first.items()), 'f1': sorted(second.items()),
                            'beta_orders': orders, 'full_degrees': degrees,
                            'branch': branch, 'common_scalar_monomial': scalar,
                        }
        rows.append({'field_order': field, 'support': [shift, shift + width - 1],
                     'endpoint_nonzero_polynomials': len(catalogue),
                     'ramification_products_m': [1, 2, 7, 32],
                     'complete_pairs': counts['diagonal'] + counts['crossed'],
                     'diagonal_cases': counts['diagonal'],
                     'crossed_cases': counts['crossed'],
                     'common_scalar_monomial_cases': counts['scalar_equalities'],
                     'distinct_degree_order_profiles': len(profile_classes)})
    record = {'status': 'All exact degree/order implications and rational cost comparisons pass.',
              'scope': 'Finite independent controls for the uniform Laurent-polynomial theorem; no finite search is used to prove its universal quantifiers.',
              'catalogues': rows, 'exact_rational_cost_comparisons': rational_checks,
              'examples': list(examples.values())}
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({key: value for key, value in record.items() if key != 'examples'}))


if __name__ == '__main__':
    main()
