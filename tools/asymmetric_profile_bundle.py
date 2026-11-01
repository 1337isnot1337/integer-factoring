#!/usr/bin/env python3
"""Existing fixed-type bundle ingredients for complete asymmetric crossings.

Raw coordinates are the 32 actual heights and the positive intrinsic scale.
At fixed supporting angle, complete profile pressure is convex in these
coordinates. The common-type root is an infimum over physical types; its
gradient must have the physical term removed before making a global plane.
This module does not treat a root gradient as a convex supporting plane.
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from asymmetric_crossing_pressure import AsymmetricPressure


def raw_from_scaled(point):
    point = np.asarray(point, dtype=float)
    t = np.exp(point[-2])
    return np.r_[t*point[:-2], t]


def scaled_from_raw(point, logtan):
    point = np.asarray(point, dtype=float)
    if point[-1] <= 0:
        raise ValueError('A positive intrinsic scale is required.')
    return np.r_[point[:-1]/point[-1], np.log(point[-1]), logtan]


class AsymmetricBundleOracle:
    def __init__(self, seed, logtan):
        self.base = AsymmetricPressure(seed)
        self.logtan = float(logtan)
        self.angle = float(np.arctan(np.exp(logtan)))

    def query(self, raw):
        raw = np.asarray(raw, dtype=float)
        rate, derivative, details = self.base.evaluate(np.r_[raw, self.angle])
        a, b = details['weights']
        r, s = (np.asarray(v) for v in details['scaled_orders'])
        value = float(a*logsumexp(-r/a)+b*logsumexp(-s/b))
        px, py = np.asarray(details['type_x']), np.asarray(details['type_y'])
        physical = np.r_[[v for _ in range(2) for leg in range(4)
                          for v in (px if leg in (0, 3) else py)[1:]], 0.]
        gradient = details['entropy_sum']*derivative[:-1]-physical
        return rate, value, gradient, details, physical
