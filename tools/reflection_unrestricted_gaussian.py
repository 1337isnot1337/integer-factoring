#!/usr/bin/env python3
"""Complete Gaussian Smith profiles without an internal regularity cone.

The Clifford Smith theorem reduces each 2**M Fock matrix to a 2*M
Bogoliubov matrix. Every physical mode may have its own real height at
each place. Field arithmetic in the native oracle is exact; its weights
and pressure optimization are floating and require exact replay before
any strict bound claim. See research/clifford-gaussian-smith-profiles.md.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import laurent_kernel, logical_heights, validate


def bogoliubov_data(interval, side):
    """Exact F2 Laurent entries, in creation-then-annihilation order.

    The variable remains w on both sides; the native oracle performs
    w -> z**-1 when side=1. Each coefficient has at most two monomials.
    """
    validate(interval)
    if side not in (0, 1):
        raise ValueError('Choose side zero or one.')
    m = interval // 2
    a = [[5 << (4 * (i & ~j)) for j in range(m)] for i in range(m)]
    c = [[5 << (4 * (i & j)) for j in range(m)] for i in range(m)]
    d = [[5 << (4 * (m - 1 - (i | j))) for j in range(m)] for i in range(m)]
    at = [list(row) for row in zip(*a)]
    identity = [[int(i == j) for j in range(m)] for i in range(m)]
    zero = [[0] * m for _ in range(m)]
    if side == 0:
        blocks = [[a, identity, c, zero], [identity, a, zero, c],
                  [d, zero, at, identity], [zero, d, identity, at]]
        block_shifts = [[-2*interval, -2*interval, -2*interval, 0],
                        [-2*interval, -2*interval, 0, -2*interval],
                        [0, 0, -2*interval, -2*interval],
                        [0, 0, -2*interval, -2*interval]]
        legs = (0, 1)
    else:
        blocks = [[c, identity, a, zero], [identity, d, zero, at],
                  [at, zero, d, identity], [zero, a, identity, c]]
        block_shifts = [[0, 2*interval, 0, 0], [2*interval, 2*interval, 0, 0],
                        [0, 0, 2*interval, 2*interval], [0, 0, 2*interval, 0]]
        legs = (0, 2)
    matrix = [[(blocks[bi][bj][i][j], 0)
               for bj in range(4) for j in range(m)]
              for bi in range(4) for i in range(m)]
    shifts = [[block_shifts[bi][bj] for bj in range(4) for j in range(m)]
              for bi in range(4) for i in range(m)]
    rows = [leg*m + i for leg in legs for i in range(m)]
    columns = [i for i in range(2*interval) if i not in rows]
    return matrix, shifts, rows, columns


def bogoliubov_oracle(interval, side, precision=None):
    matrix, shifts, rows, columns = bogoliubov_data(interval, side)
    if precision is None:
        precision = max(512, 1 << (16*interval - 1).bit_length())
    return MatrixOracle(matrix, 2*interval, shifts, precision), rows, columns


def fock_spectrum(interval, side, heights, scale=1., precision=None):
    """Return the scalar and positive factors of the full Fock spectrum.

    All full orders are minimum + subset sums of positive_factors.
    The mean formula includes the original common scalar at this place.
    """
    oracle, rows, columns = bogoliubov_oracle(interval, side, precision)
    heights = np.asarray(heights, dtype=float)
    if heights.shape != (2*interval,):
        raise ValueError('One height is required for every physical mode.')
    row = heights[rows]
    column = heights[columns]
    profile = oracle.profile(np.r_[row, -row], np.r_[column, -column],
                             scale, infinity=bool(side))
    orders = np.asarray(profile['orders'])
    if np.max(np.abs(orders + orders[::-1])) > 2e-7:
        raise ArithmeticError('The full Bogoliubov pairing was not resolved.')
    positive = orders[interval:]
    mean = (interval*interval*scale if side == 0 else 0.) + heights.sum()/2
    return {'minimum': float(mean - positive.sum()/2),
            'mean': float(mean), 'positive_factors': positive.tolist(),
            'bogoliubov_profile': profile}


class ReflectionGaussianPressure:
    """Both complete cut profiles, all physical type costs, unrestricted heights."""
    def __init__(self, interval, precision=None):
        validate(interval)
        self.interval = interval
        self.modes = interval // 2
        self.size = 2*interval
        self.cuts = []
        self.queries = 0
        self.last = None
        self.kernel = laurent_kernel(interval)
        constraints = []
        for side in (0, 1):
            oracle, rows, columns = bogoliubov_oracle(interval, side, precision)
            self.cuts.append((oracle, np.asarray(rows), np.asarray(columns)))
            for group in (rows, columns):
                for position, i in enumerate(group):
                    for j in group[position+1:]:
                        entry = self.kernel[i][j]
                        if entry:
                            order = min(entry) if side == 0 else -max(entry)
                            constraints.append((side*self.size+i, side*self.size+j, order))
        self.margin_first = np.asarray([first for first, second, order in constraints])
        self.margin_second = np.asarray([second for first, second, order in constraints])
        self.margin_orders = np.asarray([order for first, second, order in constraints])

    def natural(self, scale):
        return np.asarray([scale*h for place in ('zero', 'infinity')
                           for h in logical_heights(self.interval, place)] + [scale])

    def internal_margin(self, parameters):
        """Diagnostic only: negative values are allowed by this model."""
        parameters = np.asarray(parameters)
        return float(np.min(parameters[self.margin_first]+parameters[self.margin_second]
                            +self.margin_orders*parameters[-1]))

    def evaluate(self, parameters, gradient=False):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (2*self.size+1,) or not np.all(np.isfinite(parameters)):
            raise ValueError('Supply finite heights at both places and one scale.')
        scale = parameters[-1]
        if scale <= 0:
            raise ValueError('A positive profile scale is required.')
        # The two original common scalar orders cancel. The two cross
        # determinants together contribute M**2 to the full mean.
        log_b = -self.interval*self.interval*scale - parameters[:-1].sum()/2
        derivatives = np.r_[np.full(2*self.size, -.5), -float(self.interval*self.interval)]
        records = []
        for side, (oracle, rows, columns) in enumerate(self.cuts):
            heights = parameters[side*self.size:(side+1)*self.size]
            row, column = heights[rows], heights[columns]
            result = oracle.profile(np.r_[row, -row], np.r_[column, -column],
                                    scale, infinity=bool(side))
            orders = np.asarray(result['orders'])
            if np.max(np.abs(orders + orders[::-1])) > 2e-7:
                raise ArithmeticError('The full Bogoliubov pairing was not resolved.')
            log_b += .5*np.logaddexp(orders/2, -orders/2).sum()
            coefficients = .25*np.tanh(orders/2)
            base = np.diff([0] + result['base_prefix'])
            for coefficient, r, c, b in zip(coefficients, result['rows'], result['cols'], base):
                derivatives[side*self.size + rows[r % self.interval]] += coefficient*(1 if r < self.interval else -1)
                derivatives[side*self.size + columns[c % self.interval]] += coefficient*(1 if c < self.interval else -1)
                derivatives[-1] += coefficient*b
            records.append(result)
        heights = parameters[:-1].reshape(2, 4, self.modes)
        physical = np.asarray([heights[:, 0, :].sum(axis=0) + heights[:, 3, :].sum(axis=0),
                               heights[:, 1, :].sum(axis=0) + heights[:, 2, :].sum(axis=0)])

        def equation(rate):
            return float(log_b - rate*np.logaddexp(0, -physical/rate).sum())

        if equation(1e-8) <= 0:
            raise ArithmeticError('The positive Gibbs root is unresolved near zero entropy.')
        high = 4.
        while equation(high) > 0:
            high *= 2
        rate = brentq(equation, 1e-8, high, xtol=1e-13)
        probabilities = expit(-physical/rate)
        entropy = float(np.logaddexp(0, -physical/rate).sum()
                        + (physical/rate*probabilities).sum())
        if entropy < 1e-8:
            raise ArithmeticError('Physical entropy is too small for a stable rate.')
        for side in range(2):
            for leg in range(4):
                start = side*self.size + leg*self.modes
                derivatives[start:start+self.modes] += probabilities[0 if leg in (0, 3) else 1]
        self.queries += 1
        self.last = {'rho': float(rate/2), 'entropy': entropy,
                     'critical_pressure': equation(math.sqrt(2)),
                     'parameters': parameters.tolist(),
                     'internal_regularity_margin': self.internal_margin(parameters),
                     'bogoliubov_profiles': records}
        return (rate/2, derivatives/(2*entropy)) if gradient else rate/2
