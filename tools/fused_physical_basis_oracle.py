#!/usr/bin/env python3
"""Complete physical basis changes and weighted profiles for fused seeds.

Native profiles use exact field arithmetic and floating valuation comparisons.
The exact path uses the independently controlled F16 bordered-minor oracle.
Only rational F16 places and infinity are supported in this search adapter.
"""
from __future__ import annotations

from itertools import product

import numpy as np

from compressed_fused_six_vertex import construct_compressed
from f16_identity_operator_profiles import MUL, add, evaluate_one, inverse_matrix
from fused_six_vertex_crossings import expanded, power
from check_fused_multipoint_profiles import local_polynomial
from f16_weighted_local_oracle import F16WeightedLocalOracle
from native_matrix_smith_f4 import MatrixOracle


def identity(size):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def vandermonde(values):
    return tuple(tuple(power(value, j) for j in range(len(values))) for value in values)


def quantum_pascal(size, parameter):
    rows = [[0] * size for _ in range(size)]
    for n in range(size):
        rows[n][0] = rows[n][n] = 1
        for k in range(1, n):
            rows[n][k] = rows[n - 1][k] ^ MUL[power(parameter, n - k)][rows[n - 1][k - 1]]
    return tuple(tuple(row) for row in rows)


def basis_change(matrix, first, second):
    """Apply A_x B_y B^-T_u A^-T_v, or conjugate F by B tensor A."""
    q = len(first)
    if len(second) != q:
        raise ValueError('Both physical bases must have the seed alphabet size.')
    first_inverse, second_inverse = inverse_matrix(first), inverse_matrix(second)
    size = q * q
    left = [[MUL[second[y][v]][first[x][u]] for v, u in product(range(q), repeat=2)]
            for y, x in product(range(q), repeat=2)]
    right = [[MUL[second_inverse[y][v]][first_inverse[x][u]] for v, u in product(range(q), repeat=2)]
             for y, x in product(range(q), repeat=2)]
    temporary = [[{} for _ in range(size)] for _ in range(size)]
    result = [[{} for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for k in range(size):
            if left[i][k]:
                for j in range(size):
                    if matrix[k][j]:
                        temporary[i][j] = add(temporary[i][j], matrix[k][j], left[i][k])
    for i in range(size):
        for k in range(size):
            if temporary[i][k]:
                for j in range(size):
                    if right[k][j]:
                        result[i][j] = add(result[i][j], temporary[i][k], right[k][j])
    if not all(evaluate_one(value) == int(i == j)
               for i, row in enumerate(result) for j, value in enumerate(row)):
        raise ArithmeticError('The four physical bases did not preserve the strict endpoint.')
    return result


def cut_matrix(matrix, q, cut):
    result = [[{} for _ in range(q * q)] for _ in range(q * q)]
    for x, y, u, v in product(range(q), repeat=4):
        row, column = (q * x + y, q * u + v) if cut == 0 else (q * x + u, q * y + v)
        result[row][column] = matrix[q * y + x][q * u + v]
    return result


class FusedPhysicalBasisOracle:
    def __init__(self, m, theta, first, second, precision=512):
        self.source = construct_compressed(m, theta, normalize_joint=True)
        self.alphabet, self.degree = m + 1, self.source['degree']
        if len(first) != self.alphabet or len(second) != self.alphabet:
            raise ValueError('Physical bases must match the actual seed alphabet.')
        self.first, self.second = first, second
        self.matrix = basis_change(self.source['matrix'], first, second)
        self.cuts = [cut_matrix(self.matrix, self.alphabet, cut) for cut in (0, 1)]
        self.native, self.exact = {}, {}
        self.precision = precision

    def local_matrix(self, cut, point):
        if point is not None:
            if len(point) != 2 or point[1] != 1:
                raise ValueError('This adapter supports only rational F16 places and infinity.')
            point = point[0]
        return [[local_polynomial(value, point, self.degree) for value in row]
                for row in self.cuts[cut]]

    def weights(self, cut, heights):
        q = self.alphabet
        row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
        row_weights = [heights[row_legs[0] * q + a] + heights[row_legs[1] * q + b]
                       for a, b in product(range(q), repeat=2)]
        column_weights = [heights[column_legs[0] * q + a] + heights[column_legs[1] * q + b]
                          for a, b in product(range(q), repeat=2)]
        return row_weights, column_weights, row_legs, column_legs

    def log_profile_with_gradient(self, cut, point, heights, scale=1.0):
        q = self.alphabet
        key = cut, point
        if key not in self.native:
            # Each cut is invertible at the strict endpoint, proving its
            # generic rank. Restriction F16/F4 doubles that full rank.
            self.native[key] = MatrixOracle(expanded(self.local_matrix(cut, point)),
                                            2 * q * q, precision=self.precision)
        row_weights, column_weights, row_legs, column_legs = self.weights(cut, heights)
        profile = self.native[key].profile(np.repeat(row_weights, 2), np.repeat(column_weights, 2), scale)
        orders = np.asarray(profile['orders'])
        minimum = float(np.min(orders))
        exponentials = np.exp(-(orders - minimum))
        probabilities = exponentials / np.sum(exponentials)
        value = -minimum + float(np.log(np.sum(exponentials))) - np.log(2.0)
        gradient = np.zeros(4 * q)
        rows, columns = np.asarray(profile['rows']) // 2, np.asarray(profile['cols']) // 2
        for leg, labels in ((row_legs[0], rows // q), (row_legs[1], rows % q),
                            (column_legs[0], columns // q), (column_legs[1], columns % q)):
            gradient[leg * q:(leg + 1) * q] -= np.bincount(labels, weights=probabilities, minlength=q)
        increments = np.diff(np.r_[0.0, profile['base_prefix']])
        gradient_scale = -float(probabilities @ increments)
        orders = np.sort(orders)
        if not np.allclose(orders[::2], orders[1::2], atol=2e-8, rtol=2e-10):
            raise ArithmeticError('Restriction of scalars did not duplicate the complete Smith profile.')
        return value, gradient, gradient_scale, orders[::2]

    def profile(self, cut, point, heights, scale=1):
        key = cut, point
        if key not in self.exact:
            self.exact[key] = F16WeightedLocalOracle(self.local_matrix(cut, point), self.degree,
                                                    max_cached_states=96)
        row_weights, column_weights, _, _ = self.weights(cut, heights)
        return self.exact[key].profile(row_weights, column_weights, (0, 1), scale)['orders']

    def linear_heights(self, slopes):
        return tuple(slope * a for slope in slopes for a in range(self.alphabet))

    def type_shift(self, heights, type_x, type_y):
        q = self.alphabet
        return (sum(type_x[a] * (heights[a] + heights[3 * q + a]) for a in range(q))
                + sum(type_y[a] * (heights[q + a] + heights[2 * q + a]) for a in range(q)))
