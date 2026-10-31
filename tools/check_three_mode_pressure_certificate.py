#!/usr/bin/env python3
"""Exact rational pressure certificate for the raw three-mode crossing.

The checker uses integer polynomial arithmetic, exhaustive minors, and
outward-rounded dyadic intervals. No floating-point comparison certifies
the pressure inequality.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from check_schur_downset_crossing import sparse_kernel
from check_weighted_smith_oracle import enumerate_minors
from weighted_smith_oracle import WeightedSmithOracle


ROOT = Path(__file__).resolve().parent.parent


class DyadicIntervals:
    def __init__(self, precision):
        self.precision = precision
        self.scale = 1 << precision

    @staticmethod
    def ceil_div(a, b):
        return (a + b - 1) // b

    def rational(self, numerator, denominator):
        assert numerator > 0 and denominator > 0
        scaled = numerator * self.scale
        return scaled // denominator, self.ceil_div(scaled, denominator)

    def multiply(self, a, b):
        assert min(a + b) >= 0
        return (a[0] * b[0] // self.scale,
                self.ceil_div(a[1] * b[1], self.scale))

    def inverse(self, a):
        assert 0 < a[0] <= a[1]
        square = self.scale * self.scale
        return square // a[1], self.ceil_div(square, a[0])

    def power(self, a, exponent):
        if exponent < 0:
            a, exponent = self.inverse(a), -exponent
        result = (self.scale, self.scale)
        while exponent:
            if exponent & 1:
                result = self.multiply(result, a)
            exponent >>= 1
            if exponent:
                a = self.multiply(a, a)
        return result

    def add_one(self, a):
        return a[0] + self.scale, a[1] + self.scale

    def product(self, terms):
        result = (self.scale, self.scale)
        for term in terms:
            result = self.multiply(result, term)
        return result


def check_interval_arithmetic():
    arithmetic = DyadicIntervals(48)
    checked = 0
    for value in (Fraction(3, 2), Fraction(7, 9), Fraction(33, 32), Fraction(101, 100)):
        interval = arithmetic.rational(value.numerator, value.denominator)
        for exponent in (-31, -4, -1, 0, 1, 7, 29):
            lower, upper = arithmetic.power(interval, exponent)
            exact = value ** exponent
            assert Fraction(lower, arithmetic.scale) <= exact <= Fraction(upper, arithmetic.scale)
            checked += 1
    root = Fraction(33, 32)
    interval = arithmetic.rational(root.numerator, root.denominator)
    terms = [arithmetic.add_one(arithmetic.power(interval, exponent)) for exponent in (-7, -1, 0, 3, 9)]
    product = arithmetic.power(arithmetic.product(terms), 5)
    exact = Fraction(1)
    for exponent in (-7, -1, 0, 3, 9):
        exact *= 1 + root ** exponent
    exact **= 5
    assert Fraction(product[0], arithmetic.scale) <= exact <= Fraction(product[1], arithmetic.scale)
    return checked + 1


def exact_profiles(certificate, independent=True):
    indices = tuple(certificate['indices'])
    ambient = certificate['ambient']
    ramification = certificate['ramification']
    heights = certificate['heights']
    assert indices == (0, 1, 2) and ambient == 4
    assert ramification > 0 and len(heights) == 2
    n = len(indices)
    assert all(len(end) == 8 * n and all(isinstance(h, int) for h in end) for end in heights)
    kernel = sparse_kernel(indices, ambient)
    all_orders, profile_blocks = [], []
    regularity_checks = minor_checks = 0
    minimum_internal_slack = None
    for endpoint, cut in ((0, 0), (1, 1)):
        groups = ((1, 3), (2, 0)) if cut == 0 else ((1, 0), (2, 3))
        hs = heights[endpoint]
        for group in groups:
            for a in group:
                for b in group:
                    for i in range(n):
                        for j in range(n):
                            value = kernel[a * n + i][b * n + j]
                            if not value:
                                continue
                            valuation = 1 - value.bit_length() if endpoint else (value & -value).bit_length() - 1
                            slack = ramification * valuation + hs[a * n + i] + hs[4 * n + b * n + j]
                            assert slack >= 0
                            minimum_internal_slack = slack if minimum_internal_slack is None else min(minimum_internal_slack, slack)
                            regularity_checks += 1
        for direction in (0, 1):
            rows = tuple(g * n + i for g in groups[direction] for i in range(n))
            columns = tuple(g * n + i for g in groups[1 - direction] for i in range(n))
            matrix = [[kernel[i][j] for j in columns] for i in rows]
            rw = [hs[i] for i in rows]
            cw = [hs[4 * n + i] for i in columns]
            oracle = WeightedSmithOracle(matrix, infinity=bool(endpoint))
            orders, _ = oracle.profile(rw, cw, ramification)
            if independent:
                minima = {}
                for size, row_mask, column_mask, value in enumerate_minors(matrix):
                    weight = ramification * oracle.valuation(value)
                    weight += sum(rw[i] for i in range(2 * n) if row_mask >> i & 1)
                    weight += sum(cw[j] for j in range(2 * n) if column_mask >> j & 1)
                    minima[size] = min(minima.get(size, weight), weight)
                    minor_checks += 1
                assert all(sum(orders[:size]) == minima[size] for size in range(1, 2 * n + 1))
            all_orders.extend(orders)
            profile_blocks.append({'endpoint': 'infinity' if endpoint else 'zero',
                                   'cut': cut, 'direction': direction, 'orders': list(orders)})
    weights_a = [sum(h[2 * n + i] + h[3 * n + i] + h[4 * n + i] + h[5 * n + i]
                     for h in heights) for i in range(n)]
    weights_b = [sum(h[i] + h[n + i] + h[6 * n + i] + h[7 * n + i]
                     for h in heights) for i in range(n)]
    return all_orders, weights_a + weights_b, {
        'profile_blocks': profile_blocks, 'physical_mode_weights_A': weights_a,
        'physical_mode_weights_B': weights_b, 'internal_regularity_checks': regularity_checks,
        'minimum_internal_slack': minimum_internal_slack, 'independent_nonzero_minor_checks': minor_checks}


def pressure_interval(certificate, orders, weights):
    numerator = certificate['target_R_numerator']
    denominator = certificate['target_R_denominator']
    root_numerator = certificate['root_numerator']
    root_denominator = certificate['root_denominator']
    assert numerator > 0 and denominator > 0 and root_numerator > root_denominator > 0
    arithmetic = DyadicIntervals(certificate['dyadic_precision'])
    root = arithmetic.rational(root_numerator, root_denominator)
    left_base = arithmetic.product(arithmetic.add_one(arithmetic.power(root, -numerator * order))
                                   for order in orders)
    right_base = arithmetic.product(arithmetic.add_one(arithmetic.power(root, -denominator * weight))
                                    for weight in weights)
    lhs = arithmetic.power(left_base, denominator)
    rhs = arithmetic.power(right_base, numerator)
    return lhs, rhs


def check(certificate):
    interval_controls = check_interval_arithmetic()
    orders, weights, details = exact_profiles(certificate)
    lhs, rhs = pressure_interval(certificate, orders, weights)
    assert lhs[1] < rhs[0], 'outward rational intervals did not prove the pressure inequality'
    numerator = certificate['target_R_numerator']
    denominator = certificate['target_R_denominator']
    rho = Fraction(numerator, 4 * denominator)
    return {'status': 'pass', 'certified_rho_upper': str(rho), 'strict_inequality': True,
            'post_half_target_reached': rho * rho < Fraction(1, 2),
            'scope': 'A strict improvement of the previous Gaussian pressure control; no improved factoring exponent.',
            **details, 'dyadic_precision': certificate['dyadic_precision'],
            'independent_fraction_interval_controls': interval_controls,
            'lhs_interval_numerators': list(map(str, lhs)),
            'rhs_interval_numerators': list(map(str, rhs)),
            'strict_integer_gap': str(rhs[0] - lhs[1]),
            'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                              for name in ('tools/check_three_mode_pressure_certificate.py',
                                           'tools/weighted_smith_oracle.py',
                                           'tools/check_schur_downset_crossing.py',
                                           'tools/check_fixed_mode_crossings.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path,
                        default=ROOT / 'evidence/three-mode-pressure-certificate.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check(json.loads(args.certificate.read_text()))
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
