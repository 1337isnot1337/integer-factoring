#!/usr/bin/env python3
"""Deterministic all-subset Gram scalings over characteristic-two towers.

The mathematical construction works over every finite characteristic-two
base. This implementation represents bases in the explicit quadratic tower
over F2. Field elements are bit vectors; no conjugation is used in Gram
matrices. The construction does not enumerate column subsets.
"""
from __future__ import annotations


class TowerField:
    """F_(2^(2^level)) using z_j^2 + z_j + product_(i<j) z_i.

    The product of all generators has absolute trace one at every level.
    Thus each successive quadratic is irreducible. Low coordinates retain
    the preceding field, giving canonical compatible embeddings.
    """

    def __init__(self, level: int):
        if type(level) is not int or level < 0:
            raise ValueError('Tower level must be a nonnegative integer.')
        self.level = level
        self.bits = 1 << level
        self.order = 1 << self.bits

    @staticmethod
    def _mul(a: int, b: int, level: int) -> int:
        if not a or not b:
            return 0
        if a == 1:
            return b
        if b == 1:
            return a
        if level == 0:
            return a & b
        half = 1 << (level - 1)
        mask = (1 << half) - 1
        a0, a1 = a & mask, a >> half
        b0, b1 = b & mask, b >> half
        low = TowerField._mul(a0, b0, level - 1)
        high = TowerField._mul(a1, b1, level - 1)
        cross = TowerField._mul(a0 ^ a1, b0 ^ b1, level - 1) ^ low
        alpha = 1 << (half - 1)
        return (low ^ TowerField._mul(high, alpha, level - 1)) | (cross << half)

    def mul(self, a: int, b: int) -> int:
        return self._mul(a, b, self.level)

    def pow(self, a: int, exponent: int) -> int:
        if exponent < 0:
            return self.pow(self.inv(a), -exponent)
        answer = 1
        while exponent:
            if exponent & 1:
                answer = self.mul(answer, a)
            exponent >>= 1
            if exponent:
                a = self.mul(a, a)
        return answer

    def inv(self, a: int) -> int:
        if not a:
            raise ZeroDivisionError('Zero has no field inverse.')
        return self.pow(a, self.order - 2)

    def absolute_trace(self, a: int) -> int:
        total, term = 0, a
        for _ in range(self.bits):
            total ^= term
            term = self.mul(term, term)
        if total not in (0, 1):
            raise ArithmeticError('Absolute trace is outside F2.')
        return total


def rank(matrix: list[list[int]], field: TowerField) -> int:
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError('Matrix is not rectangular.')
    A = [row[:] for row in matrix]
    cursor = 0
    for column in range(width):
        pivot = next((i for i in range(cursor, len(A)) if A[i][column]), None)
        if pivot is None:
            continue
        A[cursor], A[pivot] = A[pivot], A[cursor]
        inverse = field.inv(A[cursor][column])
        A[cursor] = [field.mul(x, inverse) for x in A[cursor]]
        for i in range(cursor + 1, len(A)):
            scale = A[i][column]
            if scale:
                A[i] = [x ^ field.mul(scale, y) for x, y in zip(A[i], A[cursor])]
        cursor += 1
        if cursor == len(A):
            break
    return cursor


def gram(matrix: list[list[int]], field: TowerField) -> list[list[int]]:
    result = [[0] * len(matrix) for _ in matrix]
    for i, row in enumerate(matrix):
        for j, other in enumerate(matrix):
            for a, b in zip(row, other):
                result[i][j] ^= field.mul(a, b)
    return result


def orthogonal_scaling(matrix: list[list[int]], base_level: int = 0) -> dict:
    """Return scales making every column-subset Gram have maximal rank.

    No subset checks or searches occur here. One column basis has weight
    zero; the k remaining columns have weights 1,...,k. Every restriction
    has a unique minimum-weight basis. The polynomial degree is at most
    a(2k-a+1)/2, where a=min(rank(matrix),k).
    """
    if not matrix or not matrix[0]:
        raise ValueError('Expected a nonempty coefficient matrix.')
    base = TowerField(base_level)
    m = len(matrix[0])
    if any(len(row) != m or any(type(a) is not int or not 0 <= a < base.order for a in row)
           for row in matrix):
        raise ValueError('Invalid coefficient matrix for the base field.')
    R = rank(matrix, base)
    chosen = []
    for j in range(m):
        columns = chosen + [j]
        if rank([[row[i] for i in columns] for row in matrix], base) > len(chosen):
            chosen.append(j)
    assert len(chosen) == R
    chosen_set = set(chosen)
    weights, next_weight = [], 1
    for j in range(m):
        if j in chosen_set:
            weights.append(0)
        else:
            weights.append(next_weight)
            next_weight += 1
    corank = m - R
    a = min(R, corank)
    degree_bound = a * (2 * corank - a + 1) // 2
    old_index_degree_bound = max(r * (m - r) for r in range(R + 1))
    assert degree_bound <= old_index_degree_bound
    relative_degree = 1
    extra_levels = 0
    while relative_degree <= degree_bound:
        relative_degree *= 2
        extra_levels += 1
    extension = TowerField(base_level + extra_levels)
    theta = 1 if extra_levels == 0 else 1 << (extension.bits // 2)
    powers = []
    current = 1
    for _ in range(corank + 1):
        powers.append(current)
        current = extension.mul(current, theta)
    scales = [powers[weight] for weight in weights]
    scaled = [[extension.mul(a, scales[i]) for i, a in enumerate(row)] for row in matrix]
    return {'base_level': base_level, 'base_bits': base.bits,
            'extension_level': extension.level, 'extension_bits': extension.bits,
            'relative_degree': relative_degree, 'degree_bound': degree_bound,
            'old_index_degree_bound': old_index_degree_bound,
            'coefficient_rank': R, 'corank': corank, 'column_basis': chosen,
            'weights': weights, 'theta': theta, 'scales': scales,
            'scaled_matrix': scaled}
