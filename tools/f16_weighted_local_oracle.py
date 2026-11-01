#!/usr/bin/env python3
"""Cached fraction-free weighted local Smith profiles over F16.

The polynomial matrix is kept literal. Cached bordered minors are independent
of the local place and of physical weights. This avoids enumerating every
minor when a fused support block has seven or more rows.
"""
from __future__ import annotations

from functools import lru_cache

from f16_identity_operator_profiles import multiply, add
from fused_six_vertex_crossings import pdiv
from check_fused_closed_profiles import factor_valuation, from_canonical, canonical


class F16WeightedLocalOracle:
    def __init__(self, matrix, homogeneous_degree, max_cached_states=2500):
        self.size = len(matrix)
        if not self.size or any(len(row) != self.size for row in matrix):
            raise ValueError('A nonempty square polynomial matrix is required.')
        self.matrix = tuple(tuple(dict(value) for value in row) for row in matrix)
        self.degree = homogeneous_degree
        if any(k < 0 or k > homogeneous_degree for row in matrix for value in row for k in value):
            raise ValueError('The homogeneous degree must contain every polynomial entry.')
        self.cache = {(0, 0): (self.matrix, {0: 1})}
        self.max_cached_states = max_cached_states
        self.cache_resets = 0

    @lru_cache(maxsize=65536)
    def value(self, point, polynomial):
        if point is None:
            return 1 - len(polynomial)
        if point == (0, 1):
            return next(i for i, coefficient in enumerate(polynomial) if coefficient)
        return factor_valuation(from_canonical(polynomial), from_canonical(point))

    def profile(self, row_weights, column_weights, point, scale=1):
        if len(row_weights) != self.size or len(column_weights) != self.size or scale < 0:
            raise ValueError('Supply complete row/column weights and a nonnegative scale.')
        if len(self.cache) > self.max_cached_states:
            self.cache = {(0, 0): (self.matrix, {0: 1})}
            self.cache_resets += 1
        row_mask = column_mask = 0
        orders, row_pivots, column_pivots, base_prefix = [], [], [], []
        for step in range(self.size):
            matrix, previous = self.cache[row_mask, column_mask]
            previous_order = self.value(point, canonical(previous))
            candidates = []
            for i in range(self.size):
                if row_mask >> i & 1:
                    continue
                for j in range(self.size):
                    if column_mask >> j & 1 or not matrix[i][j]:
                        continue
                    base = self.value(point, canonical(matrix[i][j])) - previous_order
                    if point is None:
                        base += self.degree
                    cost = scale * base + row_weights[i] + column_weights[j]
                    candidates.append((cost, i, j))
            if not candidates:
                raise ArithmeticError('The polynomial matrix is not generically invertible.')
            cost, i, j = min(candidates)
            pivot = matrix[i][j]
            orders.append(cost)
            row_pivots.append(i)
            column_pivots.append(j)
            prefix = self.value(point, canonical(pivot))
            if point is None:
                prefix += self.degree * (step + 1)
            base_prefix.append(prefix)
            next_rows, next_columns = row_mask | (1 << i), column_mask | (1 << j)
            if step + 1 < self.size and (next_rows, next_columns) not in self.cache:
                updated = [[{} for _ in range(self.size)] for _ in range(self.size)]
                for a in range(self.size):
                    if next_rows >> a & 1:
                        continue
                    for b in range(self.size):
                        if next_columns >> b & 1:
                            continue
                        numerator = add(multiply(pivot, matrix[a][b]),
                                        multiply(matrix[a][j], matrix[i][b]))
                        quotient, remainder = pdiv(numerator, previous)
                        if remainder:
                            raise ArithmeticError('The bordered-minor division was not exact.')
                        updated[a][b] = quotient
                self.cache[next_rows, next_columns] = (tuple(tuple(row) for row in updated), pivot)
            row_mask, column_mask = next_rows, next_columns
        if all(isinstance(value, (int,)) for value in orders):
            assert orders == sorted(orders)
        else:
            assert all(float(orders[i]) <= float(orders[i + 1]) + 1e-9
                       for i in range(len(orders) - 1))
        return {'orders': tuple(orders), 'rows': tuple(row_pivots),
                'columns': tuple(column_pivots), 'base_prefix': tuple(base_prefix),
                'determinant': pivot}
