#!/usr/bin/env python3
"""Gibbs-optimized balanced pressure with dual mirrored compound gauges."""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp, softmax


def evaluate(oracle, heights, scale):
    """Use one gauge and its y/u-swapped partner at the opposite place.

    After dividing numerator and denominator by two, the objective is
    (log Z + <pX,cX> + <pY,cY>)/(HX+HY), where
    cX=gx+gv+s*EX/2 and cY=gy+gu+s*EY/2.
    """
    heights = np.asarray(heights, dtype=float)
    raw, gradient, scale_gradient, orders = oracle.log_profile_with_gradient(heights, scale)
    offsets = oracle.offsets
    cx = heights[offsets[0]:offsets[1]] + heights[offsets[3]:offsets[4]] + scale * oracle.energies_x / 2
    cy = heights[offsets[1]:offsets[2]] + heights[offsets[2]:offsets[3]] + scale * oracle.energies_y / 2
    minx, miny = float(min(cx)), float(min(cy))
    constant = raw + minx + miny

    def equation(rate):
        return constant - rate * (float(logsumexp(-(cx - minx) / rate))
                                  + float(logsumexp(-(cy - miny) / rate)))

    if equation(1e-10) <= 0:
        raise ArithmeticError('The mirrored pressure has no positive common-type root.')
    high = 4.0
    while equation(high) > 0:
        high *= 2
        if high > 1e6:
            raise ArithmeticError('The mirrored common-type root was not bracketed.')
    rho = float(brentq(equation, 1e-10, high, xtol=2e-13))
    px, py = softmax(-cx / rho), softmax(-cy / rho)
    entropy = -sum(float(value) * np.log(float(value)) for value in np.r_[px, py] if value > 0)
    if entropy < 1e-10:
        raise ArithmeticError('The mirrored common type has degenerate entropy.')
    gradient = (gradient + np.r_[px, py, py, px]) / entropy
    scale_gradient = (scale_gradient + .5 * (px @ oracle.energies_x + py @ oracle.energies_y)) / entropy
    details = {'rho': rho, 'scale': float(scale), 'heights': heights.tolist(),
               'type_x': px.tolist(), 'type_y': py.tolist(), 'entropy_sum': float(entropy),
               'orders': orders.tolist(), 'local_log_sum': raw}
    return rho, gradient, float(scale_gradient), details
