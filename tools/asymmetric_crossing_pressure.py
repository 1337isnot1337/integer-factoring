#!/usr/bin/env python3
"""Exploratory supporting-line pressure for a complete two-cut crossing.

For positive a,b of Euclidean norm one, this computes the Gibbs-optimized
bound on a*r+b*c, with r,c normalized by twice the external label width.
A value below one is a candidate for a stronger rectangular recursion,
not a factoring theorem or an exact pressure certificate by itself.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp, softmax

from native_rectangular_smith_f4 import RectangularOracle


class AsymmetricPressure:
    def __init__(self, seed, cut=1):
        self.seed, self.cut = seed, cut
        self.oracle = RectangularOracle(seed['values'], seed['q_x'], seed['q_u'])
        self.qx, self.qy = seed['q_x'], seed['q_u']
        self.n = self.qx * self.qy
        self.stride = 2 * (self.qx - 1) + 2 * (self.qy - 1)
        self.alphabets = (self.qx, self.qy, self.qy, self.qx)
        self.offsets = np.cumsum([0] + [q - 1 for q in self.alphabets[:-1]])
        self.parameters = 2 * self.stride + 2

    def evaluate(self, variables, gradient=True):
        variables = np.asarray(variables, dtype=float)
        if len(variables) != self.parameters:
            raise ValueError('incorrect number of height, scale, and angle variables')
        params, angle = variables[:-1], float(variables[-1])
        a, b = math.cos(angle), math.sin(angle)
        if min(a, b) <= 0:
            raise ValueError('both supporting-line coefficients must be positive')
        prefixes = self.oracle.prefixes(params, self.cut)
        orders = np.asarray(prefixes['orders'], float).reshape(2, self.n)
        features = np.zeros((2, self.n, len(params)))
        for side in range(2):
            chosen_cut = self.cut if side == 0 else 1 - self.cut
            previous = 0
            for i in range(self.n):
                position = side * self.n + i
                row, col = prefixes['rows'][position], prefixes['cols'][position]
                x, second = divmod(row, self.qy)
                third, v = divmod(col, self.qx)
                letters = (x, second, third, v) if chosen_cut == 0 else (x, third, second, v)
                for leg, symbol in enumerate(letters):
                    if symbol:
                        features[side, i, side * self.stride + self.offsets[leg] + symbol - 1] += 1
                cumulative = prefixes['base_prefix'][position]
                features[side, i, -1] = cumulative - previous
                previous = cumulative
        gx, gy = np.zeros(self.qx), np.zeros(self.qy)
        for side in range(2):
            for leg, alphabet in enumerate(self.alphabets):
                values = gx if leg in (0, 3) else gy
                begin = side * self.stride + self.offsets[leg]
                values[1:] += params[begin:begin + alphabet - 1]
        scales = (a, b)
        logs = [float(logsumexp(-orders[side] / scales[side])) for side in range(2)]
        channel_probabilities = [softmax(-orders[side] / scales[side]) for side in range(2)]
        raw = a * logs[0] + b * logs[1]
        minx, miny = float(min(gx)), float(min(gy))
        constant = raw + minx + miny

        def equation(rate):
            return constant - rate * (logsumexp(-(gx - minx) / rate)
                                      + logsumexp(-(gy - miny) / rate))

        if equation(1e-10) <= 0:
            raise ArithmeticError('no positive common-type root')
        high = 4.0
        while equation(high) > 0:
            high *= 2
            if high > 1e6:
                raise ArithmeticError('common-type root was not bracketed')
        rate = float(brentq(equation, 1e-10, high, xtol=2e-12))
        px, py = softmax(-gx / rate), softmax(-gy / rate)
        entropy = -sum(float(p) * math.log(float(p)) for p in list(px) + list(py) if p > 0)
        if entropy < 1e-10:
            raise ArithmeticError('degenerate type entropy')
        data = {'support_value': rate, 'weights': [a, b], 'angle': angle,
                'type_x': px.tolist(), 'type_y': py.tolist(), 'entropy_sum': entropy,
                'scaled_orders': orders.tolist(), 'scaled_energy_x': gx.tolist(),
                'scaled_energy_y': gy.tolist()}
        if not gradient:
            return rate, data
        derivative = np.zeros(len(variables))
        for side in range(2):
            derivative[:-1] -= channel_probabilities[side] @ features[side]
            for leg, alphabet in enumerate(self.alphabets):
                probability = px if leg in (0, 3) else py
                begin = side * self.stride + self.offsets[leg]
                derivative[begin:begin + alphabet - 1] += probability[1:]
        entropies = [logs[side] + float(channel_probabilities[side] @ orders[side]) / scales[side]
                     for side in range(2)]
        derivative[-1] = -b * entropies[0] + a * entropies[1]
        return rate, derivative / entropy, data
