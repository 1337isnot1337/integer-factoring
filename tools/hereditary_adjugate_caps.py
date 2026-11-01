#!/usr/bin/env python3
"""Graded adjugate caps inherited from every smaller physical alphabet.

This applies the existing graded-adjugate theorem after restriction and
restoration of core minima. These are necessary scalar profile constraints.
"""
from fractions import Fraction as F


def constraints(deficit, energy_x, energy_y):
    ex, ey = (tuple(sorted(map(F, values))) for values in (energy_x, energy_y))
    q, d = len(ex), F(deficit)
    if q < 2 or len(ey) != q or d < 0 or ex[0] != 0 or ey[0] != 0:
        raise ValueError('Use equally sized nonnegative paired lists with minimum zero and d>=0.')
    records = []
    for m in range(2, q+1):
        x, y = ex[:m], ey[:m]
        n = m*m
        constant = (n-1)*d+m*(sum(x)+sum(y))
        cells = sorted(a+b for a in x for b in y)
        for k in range(1, n+1):
            records.append(dict(family='hereditary adjugate cell', physical_size=m, size=k,
                                first_prefixes=[(n-1, k)], second_prefixes=[(2, k-1)],
                                rhs=k*constant-sum(cells[-k:])))
        for l in range(1, m+1):
            records.append(dict(family='hereditary adjugate shared', physical_size=m, size=l,
                                first_prefixes=[(n, l-1), (n-l, 1)], second_prefixes=[(2, l-1)],
                                rhs=l*constant-max(sum(x[-l:])+l*y[-1], sum(y[-l:])+l*x[-1])))
    return records


def inequalities(deficit, first, second, energy_x, energy_y):
    q = len(energy_x)
    r, s = tuple(map(F, first)), tuple(map(F, second))
    if (q < 2 or len(r) != q*q or len(s) != q*q or r != tuple(sorted(r))
            or s != tuple(sorted(s)) or r[0] or s[0]):
        raise ValueError('Use two complete sorted normalized profiles.')
    result = []
    for cap in constraints(deficit, energy_x, energy_y):
        for orientation, (a, b) in enumerate(((r, s), (s, r))):
            lhs = sum(coefficient*sum(a[:index]) for index, coefficient in cap['first_prefixes'])
            lhs += sum(coefficient*sum(b[:index]) for index, coefficient in cap['second_prefixes'])
            result.append(dict(family=cap['family'], physical_size=cap['physical_size'], size=cap['size'],
                               orientation=orientation, lhs=lhs, rhs=cap['rhs'], slack=cap['rhs']-lhs))
    return result
