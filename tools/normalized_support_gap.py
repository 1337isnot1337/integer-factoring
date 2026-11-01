#!/usr/bin/env python3
"""A search objective with the same threshold as the actual support root.

For y=log(tan theta), 1/(a*b)=2*cosh(y). This normalizes the
distance to support one, without changing its sign or certifying it.
"""
from __future__ import annotations

import math
import numpy as np


def normalized_gap(rate, point, derivative=None):
    y = float(point[-1])
    multiplier = 2*math.cosh(y)
    value = (float(rate)-1)*multiplier
    if derivative is None:
        return value
    result = np.asarray(derivative, float)*multiplier
    result[-1] += (float(rate)-1)*2*math.sinh(y)
    return value, result
