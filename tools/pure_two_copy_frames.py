#!/usr/bin/env python3
"""Regular Smith frames for pure two-copy corrections.

The finite-jet frame identities are checked exactly. The structured profile
routine remains a search accelerator; retained candidates require a replay
against the complete original polynomial matrices.
"""
from __future__ import annotations

from functools import lru_cache

from f4_recentered_lattices import pack, unpack, regular_frame
from native_rank_one_smith_f4 import rank_one_profile
from two_copy_crossing_couplings import reverse
from weighted_smith_oracle_f4 import polynomial_multiply, scalar_multiply


class PureTwoCopyFrames:
    def __init__(self, compiler, jet_precision=128, work_precision=512):
        if jet_precision >= work_precision:
            raise ValueError('Work precision must exceed the complete frame jet.')
        self.compiler = compiler
        self.jet_precision = jet_precision
        self.work_precision = work_precision
        self.mask = (1 << (2*jet_precision))-1
        self.frames = []
        self.diagonals = []
        self.controls = []
        for side, matrix in enumerate(compiler.local):
            source = [[pack((value, 0)) for value in row] for row in matrix]
            orders, left, right, diagonal = regular_frame(source, jet_precision)
            actual = tuple(value+compiler.scalars[side] for value in compiler.single_profiles[side])
            assert orders == actual
            self.frames.append((left, right))
            entries = [diagonal[i][i] for i in range(compiler.n)]
            self.diagonals.append([unpack(self.mul(a, b)) for a in entries for b in entries])
            self.controls.append({'side': side, 'one_copy_orders': list(orders),
                                  'complete_frame_identity_entries': compiler.n**2,
                                  'regular_frame_constant_ranks': [compiler.n, compiler.n]})

    def mul(self, a, b):
        return polynomial_multiply(a, b) & self.mask

    @lru_cache(None)
    def transformed_support(self, support, side):
        """Factor the literal binary support and transform its two vectors."""
        compiler = self.compiler
        q, n = compiler.q, compiler.n
        entries = {compiler.location(bits, side) for bits in support}
        if not entries or len(entries) != len(support):
            raise ValueError('Use a nonempty, unrepeated support.')
        rows, columns = sorted({i for i, _ in entries}), sorted({j for _, j in entries})
        if entries != {(i, j) for i in rows for j in columns}:
            raise ValueError('The correction is not rank one at both physical cuts.')
        h = compiler.heights[4*side:4*side+4]
        h = [[value-min(row) for value in row] for row in h]
        if side == 0:
            row_legs, column_legs = (0, 2), (1, 3)
        else:
            row_legs, column_legs = (0, 1), (2, 3)
        row = [0]*n; column = [0]*n
        for i in rows:
            a, b = divmod(i, q)
            row[i] = 1 << (2*(h[row_legs[0]][a]+h[row_legs[1]][b]))
        for j in columns:
            a, b = divmod(j, q)
            column[j] = 1 << (2*(h[column_legs[0]][a]+h[column_legs[1]][b]))
        left, right = self.frames[side]
        transformed_row, transformed_column = [0]*n, [0]*n
        for i in range(n):
            for j in rows:
                transformed_row[i] ^= self.mul(left[i][j], row[j])
        for j in range(n):
            for i in columns:
                transformed_column[j] ^= self.mul(column[i], right[i][j])
        return tuple(transformed_row), tuple(transformed_column)

    def structured_data(self, first, second, polynomial, coefficient, side):
        compiler = self.compiler
        if coefficient not in (0, 1, 2, 3):
            raise ValueError('Supply an F4 coefficient.')
        if polynomial < 0 or polynomial.bit_length() > 2*compiler.degree+1 or polynomial.bit_count() % 2:
            raise ValueError('Use an endpoint-zero polynomial inside the complete degree bound.')
        first = tuple(map(tuple, first)); second = tuple(map(tuple, second))
        u1, v1 = self.transformed_support(first, side)
        u2, v2 = self.transformed_support(second, side)
        local = reverse(polynomial, 2*compiler.degree) if side == 0 else polynomial
        scalar = scalar_multiply(pack((local, 0)), coefficient)
        left = [unpack(self.mul(scalar, self.mul(a, b))) for a in u1 for b in u2]
        right = [unpack(self.mul(a, b)) for a in v1 for b in v2]
        return self.diagonals[side], left, right

    def profiles(self, first, second, polynomial, coefficient=1):
        result, counts = [], []
        for side in (0, 1):
            d, u, v = self.structured_data(first, second, polynomial, coefficient, side)
            answer = rank_one_profile(d, u, v, cutoff=self.jet_precision,
                                      work_precision=self.work_precision)
            result.append(tuple(order-2*self.compiler.scalars[side] for order in answer['orders']))
            counts.append(answer['pivot_counts'])
        return tuple(result), counts
