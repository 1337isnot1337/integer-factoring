#!/usr/bin/env python3
"""Complete weighted local profiles at arbitrary closed places over F4.

For an irreducible f of degree d, expand each matrix entry at alpha+t in
F4[alpha]/f, then restrict coefficients to F4. Every true local order
occurs d times. The logarithmic channel sum divides that multiplicity out;
Hermite allocation must still charge d conjugate places separately.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

import numpy as np

from affine_place_profile_oracle import polynomial
from check_fused_closed_profiles import add, gcd, polynomial_modulo, power_modulo
from f16_weighted_local_oracle import F16WeightedLocalOracle
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import F4_PRODUCT, flatten


@lru_cache(maxsize=None)
def checked_factor(point):
    """Return the monic irreducible coefficient tuple, or None at infinity."""
    if point is None:
        return None
    point = tuple(point)
    if len(point) < 2 or point[-1] != 1 or any(c not in range(4) for c in point):
        raise ValueError('Supply a monic positive-degree polynomial over F4.')
    factor = {i: value for i, value in enumerate(point) if value}
    degree = len(point) - 1
    x = polynomial_modulo({1: 1}, factor)
    current = x
    for j in range(1, degree + 1):
        current = power_modulo(current, 4, factor)
        if j <= degree // 2 and max(gcd(factor, add(current, x))) != 0:
            raise ValueError('The supplied closed-place polynomial is reducible.')
    if current != x:
        raise ValueError('The supplied closed-place polynomial is reducible.')
    return point


def alpha_powers(point, bound):
    """Powers of the residue class alpha in the basis 1,alpha,...,alpha^(d-1)."""
    point = checked_factor(tuple(point))
    d = len(point) - 1
    powers = [[1] + [0] * (d - 1)]
    for _ in range(bound):
        old = powers[-1]
        top = old[-1]
        value = [0] + old[:-1]
        if top:
            value = [a ^ F4_PRODUCT[top][b] for a, b in zip(value, point[:-1])]
        powers.append(value)
    return powers


def expanded_local_matrix(matrix, degree, point):
    """F4 polynomial multiplication matrix for the root-local expansion.

    The output ordering is original coordinate first and residue-basis
    coordinate second. This ordering is also used for repeated heights.
    """
    point = checked_factor(None if point is None else tuple(point))
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError('Supply a nonempty square polynomial matrix.')
    if any(k < 0 or k > degree or c not in range(4)
           for row in matrix for entry in row for k, c in entry.items()):
        raise ValueError('The entry coefficients or homogeneous degree are invalid.')
    if point is None:
        return [[(sum((c & 1) << (degree - k) for k, c in entry.items()),
                  sum((c >> 1) << (degree - k) for k, c in entry.items()))
                 for entry in row] for row in matrix], 1
    d = len(point) - 1
    powers = alpha_powers(point, degree + d - 1)
    result = [[(0, 0) for _ in range(n * d)] for _ in range(n * d)]
    for i, row in enumerate(matrix):
        for j, entry in enumerate(row):
            block = [[[0, 0] for _ in range(d)] for _ in range(d)]
            for k, coefficient in entry.items():
                for order in range(k + 1):
                    if order & ~k:
                        continue  # Lucas parity of binomial(k,order).
                    for column in range(d):
                        vector = powers[k - order + column]
                        for output, scalar in enumerate(vector):
                            value = F4_PRODUCT[coefficient][scalar]
                            block[output][column][0] ^= (value & 1) << order
                            block[output][column][1] ^= (value >> 1) << order
            for output, column in product(range(d), repeat=2):
                result[i * d + output][j * d + column] = tuple(block[output][column])
    return result, d


class F4ClosedPlaceOracle:
    def __init__(self, seed, precision=512):
        self.source = seed
        self.qx, self.qy = seed['q_x'], seed['q_u']
        self.q_x, self.q_y = self.qx, self.qy
        self.alphabets = (self.qx, self.qy, self.qy, self.qx)
        self.dimensions = self.alphabets
        self.offsets = np.cumsum([0] + list(self.alphabets[:-1]))
        self.alphabet = self.qx if self.qx == self.qy else None
        self.degree, self.n = seed['degree'], self.qx * self.qy
        if list(seed['dimensions']) != list(self.alphabets):
            raise ValueError('The tensor dimensions do not match its paired alphabets.')
        if len(seed['values']) != self.n * self.n:
            raise ValueError('Supply every tensor coefficient.')
        for index, (low, high) in enumerate(seed['values']):
            x = index % self.qx
            y = index // self.qx % self.qy
            u = index // (self.qx * self.qy) % self.qy
            v = index // (self.qx * self.qy * self.qy)
            if (low.bit_count() & 1, high.bit_count() & 1) != (int(x == v and y == u), 0):
                raise ArithmeticError('The tensor fails its complete strict endpoint.')
        self.cuts = [[[polynomial(value) for value in row] for row in flatten(seed, cut)[0]]
                     for cut in (0, 1)]
        self.precision = precision
        self.native, self.exact = {}, {}

    def weights(self, cut, heights):
        if len(heights) != sum(self.alphabets) or cut not in (0, 1):
            raise ValueError('Supply all physical heights and a valid cut.')
        row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
        rows = [heights[self.offsets[row_legs[0]] + a]
                + heights[self.offsets[row_legs[1]] + b]
                for a in range(self.qx) for b in range(self.qy)]
        columns = [heights[self.offsets[column_legs[0]] + a]
                   + heights[self.offsets[column_legs[1]] + b]
                   for a in range(self.qy) for b in range(self.qx)]
        return rows, columns, row_legs, column_legs

    def local_matrix(self, cut, point):
        return expanded_local_matrix(self.cuts[cut], self.degree,
                                     None if point is None else tuple(point))

    def profile(self, cut, point, heights, scale=1):
        point = checked_factor(None if point is None else tuple(point))
        if cut not in self.exact:
            self.exact[cut] = F16WeightedLocalOracle(self.cuts[cut], self.degree,
                                                    max_cached_states=96)
        rows, columns, _, _ = self.weights(cut, heights)
        return self.exact[cut].profile(rows, columns, point, scale)['orders']

    def log_profile_with_gradient(self, cut, point, heights, scale=1.0):
        point = checked_factor(None if point is None else tuple(point))
        key = cut, point
        if key not in self.native:
            matrix, multiplicity = self.local_matrix(cut, point)
            self.native[key] = (MatrixOracle(matrix, self.n * multiplicity,
                                             precision=self.precision), multiplicity)
        oracle, multiplicity = self.native[key]
        rows, columns, row_legs, column_legs = self.weights(cut, heights)
        profile = oracle.profile(np.repeat(rows, multiplicity),
                                 np.repeat(columns, multiplicity), scale)
        orders = np.asarray(profile['orders'])
        minimum = float(np.min(orders))
        exponentials = np.exp(-(orders - minimum))
        probabilities = exponentials / np.sum(exponentials)
        value = -minimum + float(np.log(np.sum(exponentials))) - np.log(float(multiplicity))
        derivative = np.zeros(sum(self.alphabets))
        chosen_rows = np.asarray(profile['rows']) // multiplicity
        chosen_columns = np.asarray(profile['cols']) // multiplicity
        for leg, labels in ((row_legs[0], chosen_rows // self.qy),
                            (row_legs[1], chosen_rows % self.qy),
                            (column_legs[0], chosen_columns // self.qx),
                            (column_legs[1], chosen_columns % self.qx)):
            start = self.offsets[leg]
            derivative[start:start + self.alphabets[leg]] -= np.bincount(
                labels, weights=probabilities, minlength=self.alphabets[leg])
        increments = np.diff(np.r_[0.0, profile['base_prefix']])
        scale_derivative = -float(probabilities @ increments)
        ordered = np.sort(orders).reshape(self.n, multiplicity)
        if not np.allclose(ordered, ordered[:, :1], atol=2e-8, rtol=2e-10):
            raise ArithmeticError('Restriction of scalars did not repeat every local channel.')
        return value, derivative, scale_derivative, ordered[:, 0]

    def type_shift(self, heights, type_x, type_y):
        if len(type_x) != self.qx or len(type_y) != self.qy:
            raise ValueError('The two paired types have incorrect dimensions.')
        probabilities = (type_x, type_y, type_y, type_x)
        return sum(probabilities[leg][a] * heights[self.offsets[leg] + a]
                   for leg in range(4) for a in range(self.alphabets[leg]))

    def linear_heights(self, slopes):
        if len(slopes) != 4:
            raise ValueError('Supply one slope per physical leg.')
        return tuple(slopes[leg] * a for leg in range(4)
                     for a in range(self.alphabets[leg]))
