#!/usr/bin/env python3
"""Complete exact profiles for all 31 single-level quadratic-inverse paths."""
from fractions import Fraction
import json
from math import exp, log
from pathlib import Path
import subprocess
import tempfile

from check_quadratic_crossings import ROOT
from check_polynomial_local_smith import native, nine_letter, prefix_rank


def main():
    records = []
    with tempfile.TemporaryDirectory(prefix='factoring-nine-letter-') as directory:
        executable = Path(directory)/'local-smith'
        subprocess.run(['clang++', '-O3', '-std=c++17', '-Wall', '-Wextra', '-Werror',
                        str(ROOT/'tools/polynomial_local_smith.cpp'), '-o', str(executable)], check=True)
        for level in range(2, 33):
            degree, row, column = nine_letter(level)
            matrices = (row, column)
            profiles = native(executable, matrices, 128)
            for i in range(2):
                precision = 128
                while not profiles[i]['complete']:
                    precision *= 2
                    assert precision <= 8192
                    profiles[i] = native(executable, [matrices[i]], precision)[0]
            mu, nu = [result['orders'] for result in profiles]
            assert len(mu) == len(nu) == 81
            assert sum(nu) == 81*degree and sum(mu) <= 81*degree
            assert min(mu) == min(nu) == 0 and max(nu) <= degree+2
            assert sum(mu[:9])+sum(nu[:9]) <= 9*degree
            assert mu[1] <= 2*(max(nu)-degree)
            first, second = mu.count(1), mu.count(2)
            assert second >= max(0, 81-first*first-7*first-10)
            for matrix, profile in zip(matrices, (mu, nu)):
                for cutoff in (1, 2, 3, 4):
                    assert prefix_rank(matrix, cutoff) == sum(max(0, cutoff-value) for value in profile)
            scale = Fraction(6, 5)
            small, large = 81*second/scale**4, second*scale**(degree-2)
            assert small**7 > 9**20 and large**7 > 9**20
            def pressure(t):
                return degree*t+log(sum(exp(-a*t) for a in mu))+log(sum(exp(-a*t) for a in nu))
            lo, hi = 0, 2
            for _ in range(90):
                a, b = (2*lo+hi)/3, (lo+2*hi)/3
                if pressure(a) < pressure(b):
                    hi = b
                else:
                    lo = a
            t = (lo+hi)/2
            records.append({'one_weight_level': level, 'degree': degree,
                            'infinity_orders': mu, 'zero_orders': nu,
                            'precision': [result['precision'] for result in profiles],
                            'independent_prefix_checks': 8,
                            'all_scale_rho_lower_bound': '5/7',
                            'approximate_rho': pressure(t)/(4*log(9)), 'approximate_scale': exp(t)})
    assert len(records) == 31
    distinct = len({(tuple(row['infinity_orders']), tuple(row['zero_orders'])) for row in records})
    minimum_second = min(row['infinity_orders'].count(2) for row in records)
    assert distinct == 14 and minimum_second == 47
    print(json.dumps({'scope': 'Complete single-level weighted-flag schedule; no factoring improvement',
                      'paths': 31, 'distinct_profiles': distinct,
                      'independent_prefix_checks': 31*8,
                      'minimum_order_two_multiplicity': minimum_second,
                      'all_scale_rho_lower_bound': '5/7', 'profiles': records}, separators=(',', ':')))


if __name__ == '__main__':
    main()
