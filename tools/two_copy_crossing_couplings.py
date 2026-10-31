#!/usr/bin/env python3
"""Literal endpoint-zero couplings between two fixed scalar crossing seeds.

The physical alphabet is the ordered pair of source letters. Both cut
matrices are formed after contraction, with the complete polynomial and
physical weights retained. Integer profile queries use the native exact
coefficient engine; a favorable numerical pressure still needs replay.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from closed_auxiliary_blocks import multiply
from native_matrix_smith_f4 import MatrixOracle


@lru_cache(None)
def binary_product(a, b):
    return multiply(a, b)


def endpoint_power(power):
    """The binary polynomial (1+z)^power, including all Lucas coefficients."""
    if power < 0:
        raise ValueError('A nonnegative power is required.')
    result = 1
    place = 1
    while power:
        if power & 1:
            result ^= result << place
        power >>= 1
        place <<= 1
    return result


def reverse(polynomial, degree):
    if polynomial.bit_length() > degree+1:
        raise ValueError('Polynomial exceeds the declared degree.')
    return int(format(polynomial, '0'+str(degree+1)+'b')[::-1], 2)


def constant_slab_tensors(generators, q):
    """Put each constant triple slab at each possible first physical label."""
    answer = []
    for basis, slab in enumerate(generators):
        if len(slab) != q**3 or any(value not in (0, 1) for value in slab):
            raise ValueError('Supply constant binary triple slabs.')
        for label in range(q):
            support = []
            for index, value in enumerate(slab):
                if value:
                    support.append((label, index % q, index//q % q, index//q**2))
            answer.append({'basis': basis, 'first_label': label, 'support': support})
    return answer


class CoupledCopies:
    def __init__(self, q, values, degree, heights, precision=512):
        if q < 2 or q**4 > 1024 or len(values) != q**4:
            raise ValueError('The source alphabet must fit the complete native matrices.')
        if len(heights) != 8 or any(len(row) != q for row in heights):
            raise ValueError('Eight complete physical height lists are required.')
        if any(not isinstance(v, int) for row in heights for v in row):
            raise ValueError('This exact profile interface requires integer heights.')
        if any(high or low < 0 or low.bit_length() > degree+1 for low, high in values):
            raise ValueError('The current source constructor requires binary polynomials.')
        self.q, self.n, self.degree = q, q*q, degree
        self.values, self.heights, self.precision = tuple(values), tuple(map(tuple, heights)), precision
        self.coordinates = tuple(tuple((index//q**j) % q for j in range(4)) for index in range(q**4))
        for (low, _), (x, y, u, v) in zip(values, self.coordinates):
            if low.bit_count() % 2 != int(x == v and y == u):
                raise ValueError('The source must have the literal crossing endpoint.')
        self.scalars, self.local, self.base, self.single_profiles = [], [], [], []
        for side in (0, 1):
            h = self.heights[4*side:4*side+4]
            scalar = -sum(min(row) for row in h)+(degree if side == 0 else 0)
            self.scalars.append(scalar)
            matrix = [[0]*self.n for _ in range(self.n)]
            for bits, (polynomial, _) in zip(self.coordinates, values):
                i, j = self.location(bits, side)
                offset = scalar+sum(row[a] for row, a in zip(h, bits))
                if side == 0:
                    polynomial = reverse(polynomial, degree)
                    offset -= degree
                assert offset >= 0
                matrix[i][j] = polynomial << offset
            self.local.append(matrix)
            native = MatrixOracle([[(p, 0) for p in row] for row in matrix], self.n,
                                  precision=precision).profile()['orders']
            self.single_profiles.append(self.integer_orders(native, scalar))
            tensor_matrix = [[(binary_product(matrix[i1][j1], matrix[i2][j2]), 0)
                              for j1 in range(self.n) for j2 in range(self.n)]
                             for i1 in range(self.n) for i2 in range(self.n)]
            self.base.append(tensor_matrix)

    def location(self, bits, side):
        x, y, u, v = bits
        return (self.q*x+u, self.q*y+v) if side == 0 else (self.q*x+y, self.q*u+v)

    @staticmethod
    def integer_orders(raw, scalar):
        if any(abs(value-round(value)) > 1e-8 for value in raw):
            raise ArithmeticError('Integer input produced a nonintegral reported order.')
        return tuple(round(value)-scalar for value in raw)

    def weighted_support(self, support, side):
        h = self.heights[4*side:4*side+4]
        result = []
        for bits in support:
            if len(bits) != 4 or any(not 0 <= a < self.q for a in bits):
                raise ValueError('Invalid physical support coordinate.')
            i, j = self.location(bits, side)
            result.append((i, j, sum(row[a] for row, a in zip(h, bits))))
        if len(set(support)) != len(support):
            raise ValueError('A support must not repeat a coordinate.')
        return result

    def matrices(self, first=(), second=(), polynomial=0, coefficient=1):
        if coefficient not in (0, 1, 2, 3):
            raise ValueError('An F4 coefficient is required.')
        if polynomial < 0 or polynomial.bit_length() > 2*self.degree+1:
            raise ValueError('Coupling polynomial exceeds the complete degree bound.')
        if polynomial.bit_count() % 2:
            raise ValueError('The coupling polynomial must vanish at one.')
        result = []
        for side in (0, 1):
            matrix = [row[:] for row in self.base[side]]
            if polynomial and coefficient:
                local_polynomial = reverse(polynomial, 2*self.degree) if side == 0 else polynomial
                for i1, j1, w1 in self.weighted_support(first, side):
                    for i2, j2, w2 in self.weighted_support(second, side):
                        shift = 2*self.scalars[side]+w1+w2-(2*self.degree if side == 0 else 0)
                        assert shift >= 0
                        term = local_polynomial << shift
                        i, j = self.n*i1+i2, self.n*j1+j2
                        low, high = matrix[i][j]
                        matrix[i][j] = (low ^ (term if coefficient & 1 else 0),
                                        high ^ (term if coefficient & 2 else 0))
            result.append(matrix)
        return result

    def profiles(self, first=(), second=(), polynomial=0, coefficient=1):
        profiles = []
        for side, matrix in enumerate(self.matrices(first, second, polynomial, coefficient)):
            raw = MatrixOracle(matrix, self.n*self.n, precision=self.precision).profile()['orders']
            profiles.append(self.integer_orders(raw, 2*self.scalars[side]))
        return tuple(profiles)

    def baseline_control(self):
        actual = self.profiles()
        for side in (0, 1):
            expected = tuple(sorted(a+b for a in self.single_profiles[side]
                                    for b in self.single_profiles[side]))
            assert tuple(sorted(actual[side])) == expected
        return actual

