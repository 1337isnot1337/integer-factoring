#!/usr/bin/env python3
"""Literal six-word reflection tensors and independent F4 split graphs.

This uses the existing split-map identity and contraction implementation.
The images need not preserve occupation, parity or the inherited height.
Every entry and endpoint is computed; Gaussian profile formulas are unused.
"""
from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp, softmax

from physical_split_crossings import F4_PRODUCT, check_split, retract
from reflection_dual_rail_gaussian import laurent_kernel

WORDS = (0, 2, 4, 8, 12, 1)
ZERO = (0,) * 5
FROBENIUS = (0, 1, 3, 2)


def deflated_kernel():
    raw = laurent_kernel(8)
    if any(p % 2 for row in raw for entry in row for p in entry):
        raise ArithmeticError('The reference kernel is not even in w.')
    return tuple(tuple(frozenset(p // 2 for p in entry) for entry in row) for row in raw)


@lru_cache(None)
def ambient_source():
    kernel = deflated_kernel()

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return frozenset((0,))
        if mask.bit_count() % 2:
            return frozenset()
        first = mask & -mask
        i, rest = first.bit_length() - 1, mask ^ first
        remaining, result = rest, set()
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            for a in kernel[i][chosen.bit_length() - 1]:
                for b in pfaffian(rest ^ chosen):
                    result.symmetric_difference_update((a + b,))
        return frozenset(result)

    values = []
    for v, u, y, x in product(WORDS, repeat=4):
        powers = pfaffian(x | y << 4 | u << 8 | v << 12)
        if any(p < -8 for p in powers):
            raise ArithmeticError('Selected source has a pole after the common z^8 shift.')
        value = sum(1 << (p + 8) for p in powers)
        if value.bit_count() % 2 != int(x == v and y == u):
            raise ArithmeticError('An ambient endpoint coefficient failed.')
        values.append((value, 0))
    return {'values': values, 'dimensions': [6] * 4, 'q_x': 6, 'q_u': 6,
            'degree': max(a.bit_length() - 1 for a, _ in values),
            'denominator_degree': 0, 'strict_endpoint': True, 'words': WORDS}


def graph_split(a, b):
    if len(a) != 5 or len(b) != 5 or any(c not in range(4) for c in (*a, *b)):
        raise ValueError('Each graph requires two five-coordinate F4 vectors.')
    right = [[int(i == j) for i in range(5)] + [a[j]] for j in range(5)]
    left = [[int(i == j) ^ F4_PRODUCT[b[i]][a[j]] for j in range(5)] + [b[i]]
            for i in range(5)]
    check_split(left, right)
    return left, right


def tensor(coefficients):
    coefficients = tuple(coefficients)
    if len(coefficients) != 20:
        raise ValueError('Supply a_x,b_x,a_y,b_y in this order, twenty coefficients.')
    lx, px = graph_split(coefficients[:5], coefficients[5:10])
    ly, py = graph_split(coefficients[10:15], coefficients[15:])
    return retract(ambient_source(), lx, px, ly, py)


def frobenius(coefficients):
    return tuple(FROBENIUS[c] for c in coefficients)


def canonical(coefficients):
    c = tuple(coefficients)
    return min(c, frobenius(c))


def sparse_neighborhood():
    vectors = [ZERO] + [tuple(c if k == i else 0 for k in range(5))
                        for i in range(5) for c in (1, 2, 3)]
    raw = set()
    for a, b in product(vectors, repeat=2):
        pair = a + b
        raw.update((pair + ZERO * 2, ZERO * 2 + pair, pair + pair))
    fixed = sum(c == frobenius(c) for c in raw)
    representatives = sorted({canonical(c) for c in raw})
    if (len(raw), fixed, len(representatives)) != (766, 106, 436):
        raise ArithmeticError('The declared sparse graph inventory changed.')
    return representatives


def natural():
    h = (5., 3., 1., 4.)
    return np.asarray((0.,) * 8 + h + h + h + h + (0.,) * 8 + (1.,))


def literal_digest(seed):
    return sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest()


class FixedProfileSupport:
    """Floating scale/angle objective for ONE fixed unit-height profile pair.

    It is valid while only scale and angle change. Height optimization must
    call the complete native oracle again. The returned root is exploratory.
    """
    def __init__(self, profiles, energy_x, energy_y):
        r, s = (np.asarray(v, float) for v in profiles)
        ex, ey = np.asarray(energy_x, float), np.asarray(energy_y, float)
        self.gap = float(min(ex) + min(ey) - min(r) - min(s))
        self.r, self.s = r - min(r), s - min(s)
        self.ex, self.ey = ex - min(ex), ey - min(ey)

    def evaluate(self, variables, gradient=True):
        logt, logtan = map(float, variables)
        t, theta = math.exp(logt), math.atan(math.exp(logtan))
        a, b = math.cos(theta), math.sin(theta)
        lr, ls = float(logsumexp(-t*self.r/a)), float(logsumexp(-t*self.s/b))
        constant = t*self.gap + a*lr + b*ls
        if constant <= 1e-13:
            raise ArithmeticError('The common-type root is numerically degenerate.')

        def equation(rate):
            return constant - rate*(logsumexp(-t*self.ex/rate) + logsumexp(-t*self.ey/rate))

        high = 4.
        while equation(high) > 0:
            high *= 2
            if high > 1e6:
                raise ArithmeticError('The fixed-profile root was not bracketed.')
        rate = float(brentq(equation, 1e-10, high, xtol=2e-12))
        if not gradient:
            return rate
        pr, ps = softmax(-t*self.r/a), softmax(-t*self.s/b)
        px, py = softmax(-t*self.ex/rate), softmax(-t*self.ey/rate)
        entropy = -sum(float(p)*math.log(float(p)) for p in (*px, *py) if p > 0)
        if entropy < 1e-10:
            raise ArithmeticError('The physical type entropy is degenerate.')
        dt = t*(self.gap-pr@self.r-ps@self.s+px@self.ex+py@self.ey)
        hr, hs = lr + t*(pr@self.r)/a, ls + t*(ps@self.s)/b
        return rate, np.asarray((dt, a*b*(-b*hr+a*hs)))/entropy
