#!/usr/bin/env python3
"""Paired Gibbs quotient with specified lower bounds on physical type mass.

All calculations in this module are numerical. The active-set formula
solves the entropy minimization for fixed real energies and a positive
trial rate; it does not certify a polynomial profile or global gauge
optimum.
"""
from __future__ import annotations

import math


def bounded_distribution(energy, rate, lower):
    import numpy as np
    energy, lower = np.asarray(energy, dtype=float), np.asarray(lower, dtype=float)
    if rate <= 0 or len(energy) != len(lower) or any(lower < 0) or sum(lower) >= 1:
        raise ValueError('Use a positive rate and nonnegative lower masses with sum below one.')
    fixed = np.zeros(len(energy), dtype=bool)
    result = np.zeros(len(energy))
    for _ in range(len(energy)+1):
        free = ~fixed
        remaining = 1.-float(sum(lower[fixed]))
        if not any(free):
            raise ArithmeticError('All type coordinates became fixed despite strictly feasible bounds.')
        logarithms = -energy[free]/rate
        weights = np.exp(logarithms-max(logarithms))
        result[fixed] = lower[fixed]
        result[free] = remaining*weights/sum(weights)
        violated = free & (result < lower-1e-15)
        if not any(violated):
            assert abs(float(sum(result))-1) < 1e-12
            return result
        fixed |= violated
    raise ArithmeticError('The bounded Gibbs active set did not terminate.')


def paired_value_gradient(oracle, parameters, lower_x, lower_y, cut=1):
    import numpy as np
    from scipy.optimize import brentq
    from scipy.special import logsumexp
    params = np.asarray(parameters, dtype=float)
    _, gradient, orders = oracle.evaluate(params, [1/oracle.A]*oracle.A,
                                          [1/oracle.B]*oracle.B, cut, True)
    gx, gy = np.zeros(oracle.A), np.zeros(oracle.B)
    for side in range(2):
        for leg in range(4):
            energy = gx if leg in (0, 3) else gy
            start = side*oracle.stride+oracle.offsets[leg]
            energy[1:] += params[start:start+len(energy)-1]
    constant = float(logsumexp(-np.asarray(orders[:oracle.n]))+
                     logsumexp(-np.asarray(orders[oracle.n:])))

    def candidate(rate):
        px = bounded_distribution(gx, rate, lower_x)
        py = bounded_distribution(gy, rate, lower_y)
        entropy = -sum(float(p)*math.log(float(p)) for p in list(px)+list(py) if p > 0)
        value = constant+float(px@gx)+float(py@gy)-rate*entropy
        return value, entropy, px, py

    low, high = 1e-10, 8.
    if candidate(low)[0] <= 0:
        raise ArithmeticError('Nonpositive bounded Gibbs root needs exact profile investigation.')
    while candidate(high)[0] > 0:
        high *= 2
        if high > 1e8:
            raise ArithmeticError('The bounded paired Gibbs root was not bracketed.')
    rate = brentq(lambda trial: candidate(trial)[0], low, high, xtol=2e-13, rtol=2e-13)
    _, entropy, px, py = candidate(rate)
    if entropy < 1e-10:
        raise ArithmeticError('The bounded type entropy is too small for stable optimization.')
    derivative = np.asarray(gradient, dtype=float)
    for side in range(2):
        for leg in range(4):
            probabilities = px if leg in (0, 3) else py
            start = side*oracle.stride+oracle.offsets[leg]
            derivative[start:start+len(probabilities)-1] += probabilities[1:]-1/len(probabilities)
    derivative /= 2*entropy
    return rate/2, derivative, px.tolist(), py.tolist(), list(orders), gx.tolist(), gy.tolist()
