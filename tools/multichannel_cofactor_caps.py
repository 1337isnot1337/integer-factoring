#!/usr/bin/env python3
"""Necessary rectangular cofactor caps using all crossed Smith channels."""
from fractions import Fraction as F


def capacity_weights(p, b):
    """n-1 selections, with at most h=min(p,b)^2 from each Smith channel."""
    if min(p, b) < 1:
        raise ValueError('Use positive physical alphabet sizes.')
    n, h = p*b, min(p, b)**2
    full, rest = divmod(n-1, h)
    result = [h]*full+([rest] if rest else [])
    return tuple(result+[0]*(n-len(result)))


def cofactor_order_lower(p, b, orders):
    if len(orders) != p*b or list(orders) != sorted(orders):
        raise ValueError('Use the complete sorted crossed Smith profile.')
    return sum(a*v for a, v in zip(capacity_weights(p, b), orders))


def linear_caps(deficit, energy_x, energy_y):
    ex, ey = sorted(map(F, energy_x)), sorted(map(F, energy_y))
    q, d = len(ex), F(deficit)
    if q < 2 or len(ey) != q or ex[0] or ey[0] or d < 0:
        raise ValueError('Use equal normalized paired alphabets and d>=0.')
    size, records = q*q, []
    for p in range(1, q+1):
        for b in range(1, q+1):
            if p == b:
                continue
            n = p*b
            x, y = ex[:p], ey[:b]
            weight = b*sum(x)+p*sum(y)
            cells = sorted(a+c for a in x for c in y)
            capacity = capacity_weights(p, b)
            for family, sizes in (('multichannel cofactor cell', range(1, n+1)),
                                  ('multichannel cofactor shared line', range(1, max(p, b)+1))):
                for k in sizes:
                    if family.endswith('cell'):
                        first = [F(k)]*(n-1)+[F(0)]
                        physical = k*weight-sum(cells[-k:])
                    else:
                        first = [F(k-(i >= n-k)) for i in range(n)]
                        lines = []
                        if k <= p:
                            lines.append(sum(x[-k:])+k*y[-1])
                        if k <= b:
                            lines.append(sum(y[-k:])+k*x[-1])
                        physical = k*weight-max(lines)
                    first += [F(0)]*(size-n)
                    second = [F(k*v) for v in capacity]+[F(0)]*(size-n)
                    multiplier = k*(n-1)
                    # Both coefficient sums already match the deficit shift.
                    assert sum(first) == sum(second) == multiplier
                    assert min(first+second) >= 0
                    records.append(dict(alphabet=(p, b), family=family, size=k,
                                        first=tuple(first), second=tuple(second),
                                        degree_multiplier=multiplier, physical_rhs=physical,
                                        rhs=multiplier*d+physical))
    return records
