#!/usr/bin/env python3
"""Necessary cofactor constraints lifted from physical subalphabets.

These linear conditions are necessary, never sufficient construction tests.
Coefficients refer to normalized complete sorted profiles with both minima
zero and normalized physical paired energies with minima zero.
"""
from __future__ import annotations

from fractions import Fraction as F


def linear_caps(deficit, energy_x, energy_y):
    ex, ey = sorted(map(F, energy_x)), sorted(map(F, energy_y))
    q = len(ex)
    if q < 2 or len(ey) != q or ex[0] or ey[0] or deficit < 0:
        raise ValueError('Equal normalized alphabets and nonnegative deficit required.')
    size = q*q
    result = []

    def put(p, b, family, k, first, second, degree_multiplier, physical_rhs):
        first += [F(0)]*(size-len(first))
        second += [F(0)]*(size-len(second))
        # Absolute-profile coefficients certify that Smith submatrix
        # interlacing can lift each condition to the full tensor.
        absolute_first, absolute_second = first[:], second[:]
        absolute_first[0] += degree_multiplier-sum(first)
        absolute_second[0] += degree_multiplier-sum(second)
        assert min(absolute_first+absolute_second) >= 0
        result.append({'alphabet': (p, b), 'family': family, 'size': k,
                       'first': tuple(first), 'second': tuple(second),
                       'degree_multiplier': degree_multiplier,
                       'physical_rhs': physical_rhs,
                       'rhs': degree_multiplier*F(deficit)+physical_rhs,
                       'absolute_first': tuple(absolute_first),
                       'absolute_second': tuple(absolute_second)})

    for p in range(2, q+1):
        n = p*p
        x, y = ex[:p], ey[:p]
        weight = p*(sum(x)+sum(y))
        cells = sorted(a+b for a in x for b in y)
        for k in range(1, n+1):
            first, second = [F(k)]*(n-1)+[F(0)], [F(0)]*n
            second[1] = k-1
            put(p, p, 'square cofactor cell', k, first, second,
                k*(n-1), k*weight-sum(cells[-k:]))
        for k in range(1, p+1):
            first = [F(k-(i >= n-k)) for i in range(n)]
            second = [F(0)]*n
            second[1] = k-1
            put(p, p, 'square cofactor shared line', k, first, second,
                k*(n-1), k*weight-max(sum(x[-k:])+k*y[-1], sum(y[-k:])+k*x[-1]))
    for p in range(1, q+1):
        for b in range(1, q+1):
            n, h = p*b, min(p, b)**2
            gap = n-1-h
            if p == b or gap <= 0:
                continue
            x, y = ex[:p], ey[:b]
            weight = b*sum(x)+p*sum(y)
            for k in range(1, max(p, b)+1):
                first = [F(k-(i >= n-k)) for i in range(n)]
                second = [F(0)]*n
                second[1] = k*gap
                physical = []
                if k <= p:
                    physical.append(sum(x[-k:])+k*y[-1])
                if k <= b:
                    physical.append(sum(y[-k:])+k*x[-1])
                put(p, b, 'rectangular cofactor shared line', k, first, second,
                    k*(n-1), k*weight-max(physical))
    return result


def inequalities(deficit, first, second, energy_x, energy_y):
    r, s = tuple(map(F, first)), tuple(map(F, second))
    n = len(energy_x)**2
    if len(r) != n or len(s) != n or r[0] or s[0] or sorted(r) != list(r) or sorted(s) != list(s):
        raise ValueError('Normalized complete sorted profiles required.')
    records = []
    for cap in linear_caps(deficit, energy_x, energy_y):
        for orientation, (a, b) in enumerate(((r, s), (s, r))):
            lhs = sum(c*v for c, v in zip(cap['first'], a))
            lhs += sum(c*v for c, v in zip(cap['second'], b))
            records.append({key: cap[key] for key in ('alphabet', 'family', 'size', 'rhs')}
                           | {'orientation': orientation, 'lhs': lhs, 'slack': cap['rhs']-lhs})
    return records
