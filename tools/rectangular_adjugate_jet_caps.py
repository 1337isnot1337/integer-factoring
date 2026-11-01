#!/usr/bin/env python3
"""Stronger rectangular cofactor prefixes from the first surviving jet.

Necessary linear constraints, with nonnegative absolute Smith coefficients.
The predecessor cap builder stays unchanged for historical exact receipts.
"""
from __future__ import annotations

from fractions import Fraction as F


def linear_caps(deficit, energy_x, energy_y):
    ex, ey = sorted(map(F, energy_x)), sorted(map(F, energy_y))
    q, size, d = len(ex), len(ex)**2, F(deficit)
    if q < 2 or len(ey) != q or ex[0] or ey[0] or d < 0:
        raise ValueError('Equal normalized paired alphabets and nonnegative deficit required.')
    records = []
    for p in range(1, q+1):
        for b in range(1, q+1):
            if p == b:
                continue
            n, h = p*b, min(p,b)**2
            nullity, gap = n-h, n-h-1
            x, y = ex[:p], ey[:b]
            weight = b*sum(x)+p*sum(y)
            cells = sorted(a+c for a in x for c in y)
            for family, sizes in (('rectangular jet cell', range(1,n+1)),
                                  ('rectangular jet shared line', range(1,max(p,b)+1))):
                for k in sizes:
                    coefficient = k*gap+max(0,k-nullity)
                    if family == 'rectangular jet cell':
                        first = [F(k)]*(n-1)+[F(0)]
                        physical = k*weight-sum(cells[-k:])
                    else:
                        first = [F(k-(i >= n-k)) for i in range(n)]
                        terms = []
                        if k <= p:
                            terms.append(sum(x[-k:])+k*y[-1])
                        if k <= b:
                            terms.append(sum(y[-k:])+k*x[-1])
                        physical = k*weight-max(terms)
                    first += [F(0)]*(size-n)
                    second = [F(0)]*size
                    second[1] = coefficient
                    multiplier = k*(n-1)
                    absolute_first, absolute_second = first[:], second[:]
                    absolute_first[0] += multiplier-sum(first)
                    absolute_second[0] += multiplier-sum(second)
                    assert min(absolute_first+absolute_second) >= 0
                    records.append({'alphabet':(p,b), 'family':family, 'size':k,
                        'first':tuple(first), 'second':tuple(second),
                        'degree_multiplier':multiplier,'physical_rhs':physical,
                        'rhs':multiplier*d+physical,'absolute_first':tuple(absolute_first),
                        'absolute_second':tuple(absolute_second),
                        'channel_nullity':nullity,'cofactor_gap':gap})
    return records
