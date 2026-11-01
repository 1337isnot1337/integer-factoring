#!/usr/bin/env python3
"""Complete weighted local profiles of an explicit square F4 crossing.

The input is the full four-letter polynomial tensor.  The strict endpoint
proves the generic rank of both cuts.  F4 places use a direct native matrix;
other rational F16 places use the checked restriction-of-scalars adapter.
An independent rational-height path remains available through ``profile``.
"""
from __future__ import annotations

import numpy as np

from check_reflection_laurent_frames import cut_matrix
from f16_identity_operator_profiles import evaluate_one
from fused_physical_basis_oracle import FusedPhysicalBasisOracle
from fused_six_vertex_crossings import gcd, pdiv, pscale
from native_matrix_smith_f4 import MatrixOracle


def polynomial(value):
    a, b = value
    return {k: ((a >> k) & 1) + 2 * ((b >> k) & 1)
            for k in range(max(a.bit_length(), b.bit_length()))
            if ((a | b) >> k) & 1}


def bitplanes(value):
    if any(coefficient not in range(4) for coefficient in value.values()):
        raise ValueError('The direct native matrix requires coefficients in F4.')
    return (sum((coefficient & 1) << k for k, coefficient in value.items()),
            sum((coefficient >> 1) << k for k, coefficient in value.items()))


def primitive_seed(seed):
    """Divide every common scalar factor and restore the strict endpoint."""
    entries = [polynomial(value) for value in seed['values']]
    common = {}
    for entry in entries:
        if entry:
            common = gcd(common, entry)
            if common == {0: 1}:
                break
    endpoint = evaluate_one(common)
    if not endpoint:
        raise ArithmeticError('A common factor vanishes at the strict endpoint.')
    values = []
    for entry in entries:
        quotient, remainder = pdiv(entry, common)
        if remainder:
            raise ArithmeticError('The scalar-content division was not exact.')
        values.append(bitplanes(pscale(quotient, endpoint)))
    degree = max(max(value[0].bit_length(), value[1].bit_length()) - 1 for value in values)
    result = dict(seed, values=values, degree=degree, monomial_content=0)
    return result, common


class AffinePlaceOracle(FusedPhysicalBasisOracle):
    def __init__(self, seed, precision=512):
        self.source = seed
        self.alphabet, self.degree = seed['q_x'], seed['degree']
        if seed['q_u'] != self.alphabet:
            raise ValueError('This adapter expects equal physical alphabets.')
        q = self.alphabet
        for index, value in enumerate(seed['values']):
            x, y, u, v = ((index // q ** leg) % q for leg in range(4))
            if evaluate_one(polynomial(value)) != int(x == v and y == u):
                raise ArithmeticError('The explicit seed fails its strict crossing endpoint.')
        self.cuts = [cut_matrix(seed, cut) for cut in (0, 1)]
        self.native, self.exact = {}, {}
        self.direct = {}
        self.precision = precision

    def log_profile_with_gradient(self, cut, point, heights, scale=1.0):
        if point is not None and (len(point) != 2 or point[1] != 1 or point[0] not in range(4)):
            return super().log_profile_with_gradient(cut, point, heights, scale)
        q = self.alphabet
        key = cut, point
        if key not in self.direct:
            local = self.local_matrix(cut, point)
            matrix = [[bitplanes(value) for value in row] for row in local]
            self.direct[key] = MatrixOracle(matrix, q * q, precision=self.precision)
        rw, cw, row_legs, column_legs = self.weights(cut, heights)
        profile = self.direct[key].profile(rw, cw, scale)
        orders = np.asarray(profile['orders'])
        minimum = float(np.min(orders))
        exponentials = np.exp(-(orders - minimum))
        probabilities = exponentials / np.sum(exponentials)
        value = -minimum + float(np.log(np.sum(exponentials)))
        derivative = np.zeros(4 * q)
        rows, columns = np.asarray(profile['rows']), np.asarray(profile['cols'])
        for leg, labels in ((row_legs[0], rows // q), (row_legs[1], rows % q),
                            (column_legs[0], columns // q), (column_legs[1], columns % q)):
            derivative[leg * q:(leg + 1) * q] -= np.bincount(labels, weights=probabilities, minlength=q)
        increments = np.diff(np.r_[0.0, profile['base_prefix']])
        scale_derivative = -float(probabilities @ increments)
        return value, derivative, scale_derivative, np.sort(orders)
