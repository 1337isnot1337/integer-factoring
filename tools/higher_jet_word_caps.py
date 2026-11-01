#!/usr/bin/env python3
"""Necessary graded adjugate prefixes from ordered Smith-channel words."""
from fractions import Fraction as F
from heapq import heappop, heappush
from itertools import accumulate

from hereditary_adjugate_caps import constraints


def word_cost_profile(q, orders):
    """First q^2 ordered-word costs, including empty word, capped at q*s2.

    The first source coordinate is unused. Nonzero raw minima are allowed
    for the shift comparison; the adjugate lemma uses normalized sources.
    Distinct ordered words retain their multiplicities, even at equal cost.
    """
    s = tuple(map(F, orders))
    n = q*q
    if q < 2 or len(s) != n or s != tuple(sorted(s)) or min(s) < 0:
        raise ValueError('Use a complete sorted nonnegative square profile.')
    cap = q*s[1]
    if not cap:
        return (F(0),)*n
    queue, result = [(F(0), 0)], []
    while queue and len(result) < n:
        cost, length = heappop(queue)
        result.append(cost)
        if length < q-1:
            for value in s[1:]:
                if cost+value < cap:
                    heappush(queue, (cost+value, length+1))
    return tuple(result+[cap]*(n-len(result)))


def inequalities(deficit, first, second, energy_x, energy_y):
    q = len(energy_x)
    r, s = tuple(map(F, first)), tuple(map(F, second))
    if (q < 2 or len(r) != q*q or len(s) != q*q or r[0] or s[0]
            or r != tuple(sorted(r)) or s != tuple(sorted(s))):
        raise ValueError('Use two complete sorted normalized profiles.')
    source = (r, s)
    prefixes = [tuple(accumulate((F(0),)+a)) for a in source]
    words = {(m, side): tuple(accumulate((F(0),)+word_cost_profile(m, a[:m*m])))
             for m in range(2, q+1) for side, a in enumerate(source)}
    result = []
    for cap in constraints(deficit, energy_x, energy_y):
        m, k = cap['physical_size'], cap['size']
        for side in (0, 1):
            lhs = sum(c*prefixes[side][i] for i, c in cap['first_prefixes'])
            lhs += words[m, 1-side][k]
            result.append(dict(physical_size=m, size=k, orientation=side,
                               family=cap['family'].replace('adjugate', 'word jet'),
                               lhs=lhs, rhs=cap['rhs'], slack=cap['rhs']-lhs))
    return result
