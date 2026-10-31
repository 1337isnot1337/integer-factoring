#!/usr/bin/env python3
"""Exact signed-energy controls for continuous charge folding."""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure


class FormalPressure:
    """Independent f(E t) and sqrt(2) f(E t/sqrt(2)) coefficient tables.

Negative arguments are removed with f_q(-x)=f_q(x)+(q−1)x.
The ground coefficient stores the multiplier of (q−1)t.
"""
    def __init__(self):
        self.terms = [defaultdict(Fraction), defaultdict(Fraction)]
        self.ground = Fraction(0)

    def term(self, place, energy, coefficient=1):
        energy, coefficient = Fraction(energy), Fraction(coefficient)
        if energy < 0:
            self.ground -= energy*coefficient
        self.terms[place][abs(energy)] += coefficient

    def add(self, other, scalar=1):
        for place in range(2):
            for energy, coefficient in other.terms[place].items():
                self.terms[place][energy] += scalar*coefficient
        self.ground += scalar*other.ground

    def canonical(self):
        return tuple(tuple(sorted((e, c) for e, c in row.items() if c)) for row in self.terms), self.ground


def formal_gaussian(n, charge):
    result = FormalPressure()
    for i in range(n):
        result.term(0, 2*i+1-charge)
        result.term(0, 2*i+1+charge)
        result.term(1, 4*i-2*n+2-charge, -1)
        result.term(1, 4*i-2*n+2+charge, -1)
    for i in range(2*n):
        result.term(0, 2*i+1-2*n)
    return result


def formal_fold(n, k, delta):
    low, high = n-k, n+k
    result = FormalPressure()
    result.add(formal_gaussian(low, delta), Fraction(1, 2))
    result.add(formal_gaussian(high, delta), Fraction(1, 2))
    for length, scalar in ((n, 2), (low, -1), (high, -1)):
        for i in range(length):
            result.term(0, 2*i+1, scalar)
    for ell in range(k):
        energy = 2*(low+1+2*ell)
        for shift in (-1, 1):
            result.term(0, energy+shift+delta, Fraction(1, 2))
            result.term(0, energy+shift-delta, Fraction(-1, 2))
        result.term(1, energy+delta, -1)
        result.term(1, energy-delta, 1)
    return result


def quadratic_add(first, second, scalar=1):
    return tuple(a+scalar*b for a, b in zip(first, second))


def quadratic_nonnegative(value):
    """Exact sign of a+b sqrt(2), with rational a,b."""
    a, b = value
    if a >= 0 and b >= 0:
        return True
    if a <= 0 and b <= 0:
        return a == b == 0
    return 2*b*b >= a*a if b > 0 else a*a >= 2*b*b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    pairs = [(n, k) for n in range(2, 33) for k in range(n-1)]
    pairs += [(64, k) for k in (0, 1, 16, 32, 55, 62)]
    identities, comparisons = 0, 0
    for n, k in pairs:
        low = n-k
        for delta in (Fraction(0), Fraction(1, 3), Fraction(1), Fraction(7, 4), Fraction(2)):
            assert formal_gaussian(n, 2*k+delta).canonical() == formal_fold(n, k, delta).canonical()
            identities += 1
            for ell in range(k):
                energy = 2*(low+1+2*ell)
                minimum = energy-delta
                assert minimum >= 4
                assert Fraction(3, 4)*minimum <= minimum-1
                assert 2*(low+ell)+1 < 2*(n+ell)+1
                comparisons += 3
    assert 27**2 > 2*19**2  # C_9>2, with no floating constant comparison.
    tail_prefixes = 0
    for n in range(3, 65):
        bound = (Fraction(2*n-2), Fraction(-2))
        for fraction in (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)):
            J = tuple(fraction*value for value in bound)
            A, C = (Fraction(0), Fraction(n+1)), (Fraction(2*n), Fraction(0))
            B = (J[1], J[0]/2)  # J/sqrt(2).
            physical = [quadratic_add(A, B, -1)]*2+[quadratic_add(A, B)]*2
            cuts = [quadratic_add(C, J, -1), C, C, quadratic_add(C, J)]
            difference = (Fraction(0), Fraction(0))
            for first, second in zip(physical, cuts):
                difference = quadratic_add(difference, quadratic_add(second, first, -1))
                assert quadratic_nonnegative(difference)
                tail_prefixes += 1
    floating, error = 0, 0.
    for digits in (2, 3, 5):
        def f(x):
            return float(logsumexp(-x*np.arange(digits)))
        for n, k in ((3, 1), (8, 3), (17, 8), (32, 23)):
            low, high = n-k, n+k
            for delta in (0., .3, 1., 2.):
                for scale in (.03, .2, 1., 3.):
                    correction = 2*sum(f((2*i+1)*scale) for i in range(n))
                    correction -= sum(f((2*i+1)*scale) for i in range(low))
                    correction -= sum(f((2*i+1)*scale) for i in range(high))
                    for ell in range(k):
                        energy = 2*(low+1+2*ell)
                        correction += .5*sum(f((energy+s+delta)*scale)-f((energy+s-delta)*scale)
                                              for s in (-1, 1))
                        correction -= math.sqrt(2)*(f((energy+delta)*scale/math.sqrt(2))
                                                   -f((energy-delta)*scale/math.sqrt(2)))
                    actual = .5*(gaussian_pressure(digits, low, delta, scale)
                                 +gaussian_pressure(digits, high, delta, scale))+correction
                    expected = gaussian_pressure(digits, n, 2*k+delta, scale)
                    difference = abs(actual-expected)
                    assert difference < 3e-10
                    error = max(error, difference)
                    floating += 1
    root = Path(__file__).resolve().parent.parent
    source = 'tools/check_folded_charged_gaussian_floor.py'
    result = {'status': 'Exact continuous-charge folding identities and analytic comparison hypotheses; finite comparisons do not prove an unrecorded all-charge floor.',
              'exact_formal_pressure_identities': identities,
              'exact_monotonicity_and_concavity_relations': comparisons,
              'exact_C9_above_two': True,
              'exact_quadratic_field_tail_prefix_inequalities': tail_prefixes,
              'independent_full_pressure_comparisons': floating,
              'maximum_floating_identity_error': error,
              'new_binary_floor': 'Every integer n≥3, t>0, and real |J|≤2n−2−2sqrt(2).',
              'folding_corollary': 'Every integer n≥9, t>0, and real |J|≤2n−16.',
              'source_sha256': {source: sha256((root/source).read_bytes()).hexdigest()}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
