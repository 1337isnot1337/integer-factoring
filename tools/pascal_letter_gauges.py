#!/usr/bin/env python3
"""Complete original-Fock Pascal profiles with arbitrary symbol potentials.

These floating weighted searches allow correlations between occupied modes.
Every proposed bound still requires an independent exact pressure certificate.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp, softmax

from native_matrix_smith_f4 import MatrixOracle
from tensor_pascal_box_crossings import exterior_tensor
from weighted_pascal_flags import natural_heights


class PascalLetterPressure:
    def __init__(self, seed, precision=512, alphabet=None):
        self.seed = seed
        self.n = seed['n']
        if self.n > 4:
            raise ValueError('The literal original-Fock search is limited to four modes.')
        tensor = exterior_tensor(seed)
        self.alphabet = list(tensor['alphabet']) if alphabet is None else list(alphabet)
        if (len(self.alphabet) < 2 or len(set(self.alphabet)) != len(self.alphabet)
                or any(type(value) is not int or not 0 <= value < (1 << self.n)
                       for value in self.alphabet)):
            raise ValueError('Supply at least two distinct ambient occupation words.')
        index = {symbol: i for i, symbol in enumerate(self.alphabet)}
        self.q = len(self.alphabet)
        self.size = 4 * self.q
        self.queries, self.last = 0, None
        self.cuts = []
        order = self.q ** 2
        for place, legs in enumerate(((0, 2, 1, 3), (0, 1, 2, 3))):
            matrix = [[(0, 0) for _ in range(order)] for _ in range(order)]
            for labels, value in tensor['values'].items():
                if any(symbol not in index for symbol in labels):
                    continue
                row = index[labels[legs[0]]] + self.q * index[labels[legs[1]]]
                column = index[labels[legs[2]]] + self.q * index[labels[legs[3]]]
                matrix[row][column] = (value, 0)
            self.cuts.append((MatrixOracle(matrix, order, precision=precision), legs))

    def lift_mode_parameters(self, parameters):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (8 * self.n + 1,):
            raise ValueError('Supply eight mode-height vectors and their common scale.')
        modes = parameters[:-1].reshape(2, 4, self.n)
        heights = np.asarray([[[sum(modes[c, leg, bit] for bit in range(self.n)
                                    if symbol >> bit & 1) for symbol in self.alphabet]
                                for leg in range(4)] for c in range(2)])
        return np.r_[heights.ravel(), parameters[-1]]

    def natural(self, scale):
        mode_parameters = np.r_[[scale * value for place in natural_heights(self.seed)
                                for value in place], scale]
        return self.lift_mode_parameters(mode_parameters)

    def evaluate(self, parameters, gradient=False):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (8 * self.q + 1,) or not np.all(np.isfinite(parameters)):
            raise ValueError('Supply eight complete symbol-height vectors and their common scale.')
        scale = parameters[-1]
        if scale <= 0:
            raise ValueError('The profile scale must be positive.')
        heights = parameters[:-1].reshape(2, 4, self.q)
        derivative = np.zeros_like(parameters)
        log_b, records = 0., []
        for place, (oracle, legs) in enumerate(self.cuts):
            row = [heights[place, legs[0], a] + heights[place, legs[1], b]
                   for b in range(self.q) for a in range(self.q)]
            column = [heights[place, legs[2], a] + heights[place, legs[3], b]
                      for b in range(self.q) for a in range(self.q)]
            result = oracle.profile(row, column, scale, infinity=bool(place))
            orders = np.asarray(result['orders'])
            log_b += float(logsumexp(-orders))
            probabilities = softmax(-orders)
            for weight, r, c in zip(probabilities, result['rows'], result['cols']):
                for leg, symbol in ((legs[0], r % self.q), (legs[1], r // self.q),
                                    (legs[2], c % self.q), (legs[3], c // self.q)):
                    derivative[place * self.size + leg * self.q + symbol] -= weight
            derivative[-1] -= float(probabilities @ np.diff([0] + result['base_prefix']))
            records.append(result)
        physical = np.asarray([heights[:, 0, :].sum(axis=0) + heights[:, 3, :].sum(axis=0),
                               heights[:, 1, :].sum(axis=0) + heights[:, 2, :].sum(axis=0)])

        def equation(rate):
            return log_b - rate * float(logsumexp(-physical / rate, axis=1).sum())

        if equation(1e-8) <= 0:
            raise ArithmeticError('The positive Gibbs root is unresolved near zero entropy.')
        high = 4.
        while equation(high) > 0:
            high *= 2
            if high > 1e8:
                raise ArithmeticError('The paired pressure root was not bracketed.')
        rate = brentq(equation, 1e-8, high, xtol=1e-13)
        probabilities = softmax(-physical / rate, axis=1)
        entropy = float((logsumexp(-physical / rate, axis=1)
                         + (physical / rate * probabilities).sum(axis=1)).sum())
        if entropy < 1e-8:
            raise ArithmeticError('The physical entropy is too small for a stable rate.')
        for place in range(2):
            for leg in range(4):
                offset = place * self.size + leg * self.q
                derivative[offset:offset + self.q] += probabilities[0 if leg in (0, 3) else 1]
        self.queries += 1
        self.last = {'rho': float(rate / 2), 'entropy': entropy,
                     'critical_pressure': equation(math.sqrt(2)),
                     'parameters': parameters.tolist(),
                     'pair_probabilities': probabilities.tolist(),
                     'complete_original_profiles': records}
        return (rate / 2, derivative / (2 * entropy)) if gradient else rate / 2
