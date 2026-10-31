#!/usr/bin/env python3
"""Exact energy identities for the all-digit exterior-charge Gaussian floor."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure


def exterior_pressure(digits, n, charge, scale):
    if digits < 2 or n < 1 or abs(charge) < 2*n or scale <= 0:
        raise ValueError('Use q≥2, n≥1, t>0, and |J|≥2n.')
    a = abs(charge)-2*n
    def f(value):
        return float(logsumexp(-value*np.arange(digits)))
    return math.fsum(2*f((2*i+1)*scale)+f((a+4*i+1)*scale)
                     +f((a+4*i+3)*scale)-2*math.sqrt(2)*f((a+4*i+2)*scale/math.sqrt(2))
                     for i in range(n))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    identities, midpoint_checks = 0, 0
    for n in range(1, 65):
        for a in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1),
                  Fraction(7, 3), Fraction(10), Fraction(2*n+1)):
            J = 2*n+a
            zero = [2*i+1+sign*J for i in range(n) for sign in (-1, 1)]
            infinity = [2*i+1-2*n for i in range(2*n)]
            physical = [4*i-2*n+2+sign*J for i in range(n) for sign in (-1, 1)]
            assert sum(-min(value, 0) for value in zero+infinity) == n*J
            assert sum(-min(value, 0) for value in physical) == n*J
            assert Counter(map(abs, zero)) == Counter(a+2*j+1 for j in range(2*n))
            assert Counter(map(abs, infinity)) == Counter([2*i+1 for i in range(n)]*2)
            assert Counter(map(abs, physical)) == Counter([a+4*i+2 for i in range(n)]*2)
            identities += 5
            for i in range(n):
                x = a+4*i+2
                assert (x-1)+(x+1) == 2*x
                assert 2*i+1 <= x/2 and x/2 > 0
                midpoint_checks += 2
    controls, error = 0, 0.
    for digits in (2, 3, 5, 8):
        for n in (1, 2, 3, 8, 32):
            for a in (0., .37, 3., 17.):
                for scale in (.03, .2, 1., 3.):
                    actual = exterior_pressure(digits, n, 2*n+a, scale)
                    for sign in (-1, 1):
                        expected = gaussian_pressure(digits, n, sign*(2*n+a), scale)
                        difference = abs(actual-expected)
                        assert difference < 3e-10
                        error = max(error, difference)
                        controls += 1
    root = Path(__file__).resolve().parent.parent
    source = 'tools/check_exterior_charged_gaussian_floor.py'
    result = {'status': 'Exact energy/ground identities supporting an analytic all-digit exterior-charge floor; numerical controls are not sign certificates.',
              'exact_thermal_and_ground_identities': identities,
              'exact_midpoint_and_monotonicity_relations': midpoint_checks,
              'independent_full_pressure_comparisons': controls,
              'maximum_floating_identity_error': error,
              'analytic_scope': 'Every integer q≥2,n≥1, every t>0, and every real |J|≥2n.',
              'source_sha256': {source: sha256((root/source).read_bytes()).hexdigest()}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
