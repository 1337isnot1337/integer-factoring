#!/usr/bin/env python3
"""Necessary nonlinear graded second-jet caps on square subalphabets."""
from fractions import Fraction as F

from hereditary_adjugate_caps import constraints


def truncated_prefix(orders, k):
    return sum(min(v, 2*orders[1]) for v in orders[:k])


def inequalities(deficit, first, second, energy_x, energy_y):
    q = len(energy_x)
    r, s = tuple(map(F, first)), tuple(map(F, second))
    if (q < 2 or len(r) != q*q or len(s) != q*q or r[0] or s[0]
            or sorted(r) != list(r) or sorted(s) != list(s)):
        raise ValueError('Use two complete sorted normalized profiles.')
    result = []
    for cap in constraints(deficit, energy_x, energy_y):
        for orientation, (a, b) in enumerate(((r, s), (s, r))):
            lhs = sum(c*sum(a[:k]) for k, c in cap['first_prefixes'])
            lhs += truncated_prefix(b, cap['size'])
            result.append(dict(physical_size=cap['physical_size'], size=cap['size'],
                               family=cap['family'].replace('adjugate', 'second jet'),
                               orientation=orientation, lhs=lhs, rhs=cap['rhs'],
                               slack=cap['rhs']-lhs))
    return result
