#!/usr/bin/env python3
"""Exact sufficient cold-gap exclusion for rational supporting weights.

The mathematical prefix hypotheses come from a strict polynomial crossing.
This small arithmetic certificate checks the exponential threshold without
floating-point logs or underflow. Failure to certify is inconclusive.
"""
from __future__ import annotations

from fractions import Fraction


def certify(alphabet_x, alphabet_y, a, b, scaled_gap, maximum_terms=256):
    a, b, x = Fraction(a), Fraction(b), Fraction(scaled_gap)
    if (not all(isinstance(q,int) and q >= 2 for q in (alphabet_x,alphabet_y))
            or min(a,b,x) <= 0 or a*a+b*b != 1):
        raise ValueError('Use alphabets at least two, positive rational norm-one weights and gap.')
    if not isinstance(maximum_terms,int) or maximum_terms < 1:
        raise ValueError('Supply a positive integer Taylor ceiling.')
    c, count = a+b, alphabet_x+alphabet_y-2
    exponent, target = x*(c-1)/c, Fraction(2*count)/c
    value = term = Fraction(1)
    degree = 0
    while degree < maximum_terms and value < target:
        degree += 1
        term *= exponent/degree
        value += term
    return {'status':'certified_positive_margin' if value >= target else 'not_certified',
            'alphabet_x':alphabet_x,'alphabet_y':alphabet_y,
            'a':str(a),'b':str(b),'scaled_physical_gap':str(x),
            'exponent':str(exponent),'exponential_target':str(target),
            'Taylor_degree':degree,'Taylor_lower_bound':str(value),
            'scope':'Strict polynomial crossing and a unique physical ground pair are required. '
                    'This excludes this cold parameter point; it is not a global bound improvement.'}
