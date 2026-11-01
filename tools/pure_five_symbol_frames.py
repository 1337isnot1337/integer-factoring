#!/usr/bin/env python3
"""Fast framed profiles for a pure F4 tensor correction of the q5 reference.

The correction is f(z) a_X b_Y c_U d_V, with f(1)=0. It has rank one
at both physical cuts. The structured engine is a search accelerator;
retained results must be replayed on the complete original matrices.
"""
from __future__ import annotations

from functools import lru_cache

from f4_recentered_lattices import pack, regular_frame, scale_row, unpack
from native_rank_one_smith_f4 import rank_one_profile


MULTIPLY = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))


def reverse(value, degree):
    if value < 0 or value.bit_length() > degree+1:
        raise ValueError('The polynomial exceeds the declared degree.')
    return int(format(value, '0'+str(degree+1)+'b')[::-1], 2)


def validate_vectors(vectors):
    vectors = tuple(map(tuple, vectors))
    if len(vectors) != 4 or any(len(row) != 5 or not any(row)
            or any(value not in (0, 1, 2, 3) for value in row) for row in vectors):
        raise ValueError('Supply four nonzero F4 vectors of length five.')
    return vectors


def scalar_polynomial(polynomial):
    polynomial = (polynomial, 0) if isinstance(polynomial, int) else tuple(polynomial)
    if (len(polynomial) != 2 or not any(polynomial)
            or any(not isinstance(value, int) or value < 0 or value.bit_length() > 19
                   or value.bit_count() % 2 for value in polynomial)):
        raise ValueError('Use a nonzero endpoint-zero F4 polynomial of degree at most18.')
    return polynomial


def literal_values(base, vectors, polynomial, coefficient=1):
    vectors = validate_vectors(vectors)
    polynomial = scalar_polynomial(polynomial)
    if coefficient not in (1, 2, 3):
        raise ValueError('Use a nonzero F4 scalar.')
    result = []
    for index, (low, high) in enumerate(base):
        value = coefficient
        for leg in range(4):
            value = MULTIPLY[value][vectors[leg][(index//5**leg) % 5]]
        first, second = scale_row(*polynomial, value)
        result.append((low ^ first, high ^ second))
    return tuple(result)


class PureFiveSymbolFrames:
    def __init__(self, base, energies, jet_precision=128, work_precision=512):
        if len(base) != 625 or len(energies) != 5 or any(value % 2 or value < 0 for value in energies):
            raise ValueError('Use the complete q5 tensor and nonnegative even paired energies.')
        if jet_precision >= work_precision:
            raise ValueError('Work precision must exceed the frame precision.')
        self.base = tuple(map(tuple, base))
        self.energies = tuple(energies)
        self.half = tuple(value//2 for value in energies)
        self.jet_precision, self.work_precision = jet_precision, work_precision
        self.mask = (1 << jet_precision)-1
        self.diagonals, self.frames, self.reference_profiles = [], [], []
        for side in (0, 1):
            matrix = [[0]*25 for _ in range(25)]
            for index, (low, high) in enumerate(base):
                x, y, u, v = ((index//5**leg) % 5 for leg in range(4))
                assert low.bit_count() % 2 == int(x == v and y == u) and high.bit_count() % 2 == 0
                if side == 0:
                    i, j, weight = 5*x+u, 5*y+v, self.half[u]+self.half[v]
                    low, high = reverse(low, 18), reverse(high, 18)
                else:
                    i, j, weight = 5*x+y, 5*u+v, self.half[x]+self.half[y]
                matrix[i][j] = pack((low << weight, high << weight))
            orders, left, right, diagonal = regular_frame(matrix, jet_precision)
            self.diagonals.append(tuple(unpack(diagonal[i][i]) for i in range(25)))
            self.frames.append((tuple(tuple(unpack(value) for value in row) for row in left),
                                tuple(tuple(unpack(value) for value in row) for row in right)))
            scalar = 18 if side == 0 else 0
            self.reference_profiles.append(tuple(order-scalar for order in orders))
        self.reference_profiles = tuple(self.reference_profiles)

    @lru_cache(maxsize=4096)
    def transformed_pair(self, side, role, first, second):
        if side not in (0, 1) or role not in ('row', 'column'):
            raise ValueError('Invalid local side or vector role.')
        left, right = self.frames[side]
        entries = []
        for i, a in enumerate(first):
            for j, b in enumerate(second):
                coefficient = MULTIPLY[a][b]
                if not coefficient:
                    continue
                if side == 0:
                    weight = self.half[j]
                elif role == 'row':
                    weight = self.half[i]+self.half[j]
                else:
                    weight = 0
                entries.append((5*i+j, coefficient, weight))
        result = []
        for out in range(25):
            low = high = 0
            for index, coefficient, weight in entries:
                a, b = left[out][index] if role == 'row' else right[index][out]
                a, b = scale_row(a, b, coefficient)
                low ^= a << weight
                high ^= b << weight
            result.append((low & self.mask, high & self.mask))
        return tuple(result)

    def structured_data(self, vectors, polynomial, coefficient, side):
        vectors = validate_vectors(vectors)
        polynomial = scalar_polynomial(polynomial)
        if coefficient not in (1, 2, 3):
            raise ValueError('Use a nonzero F4 scalar.')
        a, b, c, d = vectors
        first, second = (a, c) if side == 0 else (a, b)
        third, fourth = (b, d) if side == 0 else (c, d)
        left = self.transformed_pair(side, 'row', first, second)
        right = self.transformed_pair(side, 'column', third, fourth)
        polynomial = tuple(reverse(value, 18) for value in polynomial) if side == 0 else polynomial
        terms = [(power, ((polynomial[0] >> power) & 1) | (((polynomial[1] >> power) & 1) << 1))
                 for power in range(19) if ((polynomial[0] | polynomial[1]) >> power) & 1]
        scaled = []
        for low, high in left:
            a = b = 0
            for power, value in terms:
                first, second = scale_row(low, high, value)
                a ^= first << power
                b ^= second << power
            a, b = scale_row(a & self.mask, b & self.mask, coefficient)
            scaled.append((a, b))
        return self.diagonals[side], tuple(scaled), right

    def profiles(self, vectors, polynomial, coefficient=1):
        result, counts = [], []
        for side in (0, 1):
            diagonal, left, right = self.structured_data(vectors, polynomial, coefficient, side)
            answer = rank_one_profile(diagonal, left, right, cutoff=self.jet_precision,
                                      work_precision=self.work_precision)
            scalar = 18 if side == 0 else 0
            result.append(tuple(order-scalar for order in answer['orders']))
            counts.append(answer['pivot_counts'])
        return tuple(result), counts

