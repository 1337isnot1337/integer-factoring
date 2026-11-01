#!/usr/bin/env python3
"""Necessary graded adjugate inequalities for complete two-place profiles.

Inputs use the normalization of cell-minor-caps.md. These tests are
necessary constraints, never a construction or sufficient bound certificate.
"""
from __future__ import annotations

from fractions import Fraction as F


def inequalities(deficit, first, second, energy_x, energy_y):
    q = len(energy_x)
    r, s, ex, ey = (tuple(map(F, values)) for values in (first, second, energy_x, energy_y))
    d = F(deficit)
    if (q < 2 or len(ey) != q or len(r) != q*q or len(s) != q*q or d < 0
            or r != tuple(sorted(r)) or s != tuple(sorted(s)) or r[0] or s[0]
            or min(ex) != 0 or min(ey) != 0):
        raise ValueError('Use normalized complete sorted profiles and nonnegative paired energies.')
    ex, ey = tuple(sorted(ex)), tuple(sorted(ey))
    cells = sorted(x+y for x in ex for y in ey)
    n = q*q
    constant = (n-1)*d+q*(sum(ex)+sum(ey))
    records = []
    for orientation, (r, s) in enumerate(((r, s), (s, r))):
        for k in range(1, n+1):
            lhs = (k-1)*s[1]+k*sum(r[:-1])
            rhs = k*constant-sum(cells[-k:])
            records.append({'family': 'adjugate cell', 'orientation': orientation,
                            'size': k, 'lhs': lhs, 'rhs': rhs, 'slack': rhs-lhs})
        for k in range(1, q+1):
            lhs = (k-1)*s[1]+k*sum(r)-sum(r[-k:])
            rhs = k*constant-max(sum(ex[-k:])+k*ey[-1], sum(ey[-k:])+k*ex[-1])
            records.append({'family': 'adjugate shared line', 'orientation': orientation,
                            'size': k, 'lhs': lhs, 'rhs': rhs, 'slack': rhs-lhs})
    return records


def prefix_form(deficit, first, second, energy_x, energy_y, first_minimum):
    """Equivalent inequalities with nonnegative coefficients on absolute prefixes."""
    q, n = len(energy_x), len(first)
    minimum = F(first_minimum)
    mu = [F(v)+minimum for v in first]
    nu = [F(v)-F(deficit)-minimum for v in second]
    physical = q*(sum(map(F, energy_x))+sum(map(F, energy_y)))
    cells = sorted(F(x)+F(y) for x in energy_x for y in energy_y)
    ex, ey = sorted(map(F, energy_x)), sorted(map(F, energy_y))
    out = []
    for orientation, (mu, nu) in enumerate(((mu, nu), (nu, mu))):
        for k in range(1, n+1):
            lhs = k*sum(mu[:-1])+(k-1)*sum(nu[:2])+(k*(n-3)+2)*nu[0]
            rhs = k*physical-sum(cells[-k:])
            out.append({'family': 'adjugate cell', 'orientation': orientation,
                        'size': k, 'slack': rhs-lhs})
        for k in range(1, q+1):
            lhs = (k-1)*sum(mu)+sum(mu[:n-k])+(k-1)*sum(nu[:2])+(k*(n-3)+2)*nu[0]
            rhs = k*physical-max(sum(ex[-k:])+k*ey[-1], sum(ey[-k:])+k*ex[-1])
            out.append({'family': 'adjugate shared line', 'orientation': orientation,
                        'size': k, 'slack': rhs-lhs})
    return out
