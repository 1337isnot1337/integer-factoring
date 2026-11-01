#!/usr/bin/env python3
"""Pure physical updates whose individual letters can correlate two copies."""
from __future__ import annotations

from functools import lru_cache

from f4_recentered_lattices import pack, unpack
from native_rank_one_smith_f4 import rank_one_profile
from native_tensor_frame_f4 import sandwich
from pure_two_copy_frames import PureTwoCopyFrames
from two_copy_crossing_couplings import reverse
from weighted_smith_oracle_f4 import scalar_multiply


FIELD = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))


class EntangledTwoCopyFrames(PureTwoCopyFrames):
    def __init__(self, compiler, jet_precision=128, work_precision=512):
        super().__init__(compiler, jet_precision, work_precision)
        self.binary_frames = []
        for left, right in self.frames:
            left = [[unpack(a) for a in row] for row in left]
            right = [[unpack(a) for a in row] for row in right]
            self.binary_frames.append((left, [list(row) for row in zip(*left)],
                                       right, [list(row) for row in zip(*right)]))

    def validate_vectors(self, vectors):
        if len(vectors) != 4 or any(len(v) != self.compiler.n or not any(v)
                                   or any(c not in (0, 1, 2, 3) for c in v) for v in vectors):
            raise ValueError('Use four nonzero F4 vectors on the paired physical alphabet.')

    @lru_cache(maxsize=256)
    def weighted_pair(self, first, second, leg1, leg2, side):
        q = self.compiler.q; n = q*q
        h = self.compiler.heights[4*side:4*side+4]
        a = [x-min(h[leg1]) for x in h[leg1]]
        b = [x-min(h[leg2]) for x in h[leg2]]
        matrix = [[(0, 0)]*n for _ in range(n)]
        for i in range(n):
            x1, y1 = divmod(i, q)
            for j in range(n):
                x2, y2 = divmod(j, q)
                coefficient = FIELD[first[q*x1+x2]][second[q*y1+y2]]
                power = a[x1]+a[x2]+b[y1]+b[y2]
                matrix[i][j] = ((1 << power) if coefficient & 1 else 0,
                                (1 << power) if coefficient & 2 else 0)
        return matrix

    def transformed_vectors(self, vectors, side):
        vectors = tuple(map(tuple, vectors))
        self.validate_vectors(vectors)
        row_legs, column_legs = ((0, 2), (1, 3)) if side == 0 else ((0, 1), (2, 3))
        row = self.weighted_pair(*(vectors[k] for k in row_legs), *row_legs, side)
        column = self.weighted_pair(*(vectors[k] for k in column_legs), *column_legs, side)
        left, left_t, right, right_t = self.binary_frames[side]
        transformed_row = sandwich(left, row, left_t, self.jet_precision)
        transformed_column = sandwich(right_t, column, right, self.jet_precision)
        return ([a for row in transformed_row for a in row],
                [a for row in transformed_column for a in row])

    def profiles(self, vectors, polynomial, coefficient=1):
        if coefficient not in (1, 2, 3) or polynomial < 0 or polynomial.bit_count() % 2:
            raise ValueError('Use a nonzero F4 scalar and an endpoint-zero polynomial.')
        if polynomial.bit_length() > 2*self.compiler.degree+1:
            raise ValueError('The complete polynomial degree bound is exceeded.')
        result, counts = [], []
        for side in (0, 1):
            row, column = self.transformed_vectors(vectors, side)
            local = reverse(polynomial, 2*self.compiler.degree) if side == 0 else polynomial
            local = scalar_multiply(pack((local, 0)), coefficient)
            row = [unpack(self.mul(local, pack(a))) for a in row]
            answer = rank_one_profile(self.diagonals[side], row, column,
                                      cutoff=self.jet_precision, work_precision=self.work_precision)
            result.append(tuple(order-2*self.compiler.scalars[side] for order in answer['orders']))
            counts.append(answer['pivot_counts'])
        return tuple(result), counts

    def original_matrices(self, vectors, polynomial, coefficient=1):
        """Independently insert the literal physical four-vector correction."""
        vectors = tuple(map(tuple, vectors)); self.validate_vectors(vectors)
        if coefficient not in (1, 2, 3) or polynomial < 0 or polynomial.bit_count() % 2:
            raise ValueError('Invalid endpoint-zero correction.')
        if polynomial.bit_length() > 2*self.compiler.degree+1:
            raise ValueError('Polynomial degree exceeds the complete bound.')
        result = []
        for side in (0, 1):
            matrix = [row[:] for row in self.compiler.base[side]]
            row_legs, column_legs = ((0, 2), (1, 3)) if side == 0 else ((0, 1), (2, 3))
            row = self.weighted_pair(*(vectors[k] for k in row_legs), *row_legs, side)
            column = self.weighted_pair(*(vectors[k] for k in column_legs), *column_legs, side)
            def entries(matrix):
                values = []
                for i, pair in enumerate(a for row in matrix for a in row):
                    low, high = pair
                    active = low | high
                    if active:
                        assert active & (active-1) == 0
                        values.append((i, int(bool(low)) | (int(bool(high)) << 1), active.bit_length()-1))
                return values
            local = reverse(polynomial, 2*self.compiler.degree) if side == 0 else polynomial
            for i, a, exponent_a in entries(row):
                for j, b, exponent_b in entries(column):
                    c = FIELD[coefficient][FIELD[a][b]]
                    term = local << (exponent_a+exponent_b)
                    low, high = matrix[i][j]
                    matrix[i][j] = (low ^ (term if c & 1 else 0), high ^ (term if c & 2 else 0))
            result.append(matrix)
        return result
