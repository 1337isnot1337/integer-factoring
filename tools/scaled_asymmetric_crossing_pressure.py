#!/usr/bin/env python3
"""Stable scale/angle coordinates for the existing complete support oracle."""
from __future__ import annotations

import math
import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure


class ScaledAsymmetricPressure:
    def __init__(self, seed, cut=1):
        self.base = AsymmetricPressure(seed, cut)
        self.parameters = self.base.parameters

    def evaluate(self, variables, gradient=True):
        variables = np.asarray(variables, dtype=float)
        if variables.shape != (self.parameters,) or not np.all(np.isfinite(variables)):
            raise ValueError('Supply normalized heights, log scale and log tangent angle.')
        t, tangent = math.exp(variables[-2]), math.exp(variables[-1])
        angle = math.atan(tangent)
        native = np.r_[t*variables[:-2], t, angle]
        if gradient:
            rate, derivative, details = self.base.evaluate(native)
            result = np.r_[t*derivative[:-2], derivative[:-1] @ native[:-1],
                           derivative[-1]*math.sin(angle)*math.cos(angle)]
        else:
            rate, details = self.base.evaluate(native, gradient=False)
        details = dict(details, native_variables=native.tolist(),
                       optimizer_variables=variables.tolist())
        return (rate, result, details) if gradient else (rate, details)
