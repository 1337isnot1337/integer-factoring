#!/usr/bin/env python3
"""Exact weighted Smith orders by fraction-free minimum-valuation pivots.

Polynomial entries are bit masks in F2[z]. Row/column weights and scale may
be integers or Fractions for exact output. Floats are accepted for search,
but their comparisons are not certificates. The returned prefix minors
attain every weighted determinantal minimum for nonnegative scale.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


def polynomial_multiply(a: int, b: int) -> int:
    if a.bit_count() < b.bit_count():
        a, b = b, a
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def polynomial_exact_divide(numerator: int, denominator: int) -> int:
    if not denominator:
        raise ZeroDivisionError('zero polynomial denominator')
    if denominator == 1:
        return numerator
    answer = 0
    degree = denominator.bit_length()
    while numerator.bit_length() >= degree:
        shift = numerator.bit_length() - degree
        answer ^= 1 << shift
        numerator ^= denominator << shift
    if numerator:
        raise ArithmeticError('fraction-free division was not exact')
    return answer


@dataclass(frozen=True)
class PrefixMinor:
    size: int
    row_mask: int
    column_mask: int
    polynomial: int
    valuation: int


class WeightedSmithOracle:
    """Repeated gauge queries on one generically invertible square matrix.

    Cached states depend only on selected row/column sets. A state stores
    the next bordered minors and the preceding pivot minor; coefficients
    therefore remain polynomials throughout the elimination.
    """

    def __init__(self, matrix: Sequence[Sequence[int]], infinity: bool = False):
        self.size = len(matrix)
        if not self.size or any(len(row) != self.size for row in matrix):
            raise ValueError('a nonempty square matrix is required')
        if any(value < 0 for row in matrix for value in row):
            raise ValueError('polynomial bit masks must be nonnegative')
        self.infinity = infinity
        self.matrix = tuple(tuple(row) for row in matrix)
        self.cache = {(0, 0): (self.matrix, 1)}

    def valuation(self, value: int) -> int:
        if not value:
            raise ValueError('zero polynomial has infinite valuation')
        if self.infinity:
            return 1 - value.bit_length()
        return (value & -value).bit_length() - 1

    def clear_cache(self) -> None:
        self.cache = {(0, 0): (self.matrix, 1)}

    def profile(self, row_weights: Sequence[int | Fraction | float],
                column_weights: Sequence[int | Fraction | float],
                scale: int | Fraction | float = 1):
        n = self.size
        if len(row_weights) != n or len(column_weights) != n:
            raise ValueError('one weight per row and column is required')
        if scale < 0:
            raise ValueError('valuation scale must be nonnegative')
        row_mask = column_mask = 0
        orders = []
        prefixes = []
        for step in range(n):
            matrix, previous = self.cache[row_mask, column_mask]
            previous_order = self.valuation(previous)
            candidates = []
            for i in range(n):
                if row_mask >> i & 1:
                    continue
                for j in range(n):
                    if column_mask >> j & 1 or not matrix[i][j]:
                        continue
                    order = (scale * (self.valuation(matrix[i][j]) - previous_order)
                             + row_weights[i] + column_weights[j])
                    candidates.append((order, i, j))
            if not candidates:
                raise ValueError('matrix is not generically invertible')
            order, i, j = min(candidates)
            orders.append(order)
            pivot = matrix[i][j]
            next_rows = row_mask | (1 << i)
            next_columns = column_mask | (1 << j)
            prefixes.append(PrefixMinor(step + 1, next_rows, next_columns,
                                        pivot, self.valuation(pivot)))
            if step + 1 < n and (next_rows, next_columns) not in self.cache:
                updated = [[0] * n for _ in range(n)]
                for a in range(n):
                    if next_rows >> a & 1:
                        continue
                    for b in range(n):
                        if next_columns >> b & 1:
                            continue
                        numerator = (polynomial_multiply(pivot, matrix[a][b])
                                     ^ polynomial_multiply(matrix[a][j], matrix[i][b]))
                        updated[a][b] = polynomial_exact_divide(numerator, previous)
                self.cache[next_rows, next_columns] = (tuple(tuple(row) for row in updated), pivot)
            row_mask, column_mask = next_rows, next_columns
        return tuple(orders), tuple(prefixes)
