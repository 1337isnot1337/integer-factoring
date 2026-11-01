#!/usr/bin/env python3
"""Streamed characteristic-two squarefree product queries via one Pfaffian.

The implementation uses explicitly represented finite binary fields. The
identity itself holds over every field of characteristic two. It evaluates
one scalar entry and makes no coefficient-matrix rank assertion.
"""

from dataclasses import dataclass


class BinaryField:
    """GF(2^degree), represented modulo a supplied irreducible polynomial."""

    def __init__(self, modulus):
        self.modulus = modulus
        self.degree = modulus.bit_length() - 1
        self.size = 1 << self.degree

    def mul(self, first, second):
        result = 0
        while second:
            if second & 1:
                result ^= first
            second >>= 1
            first <<= 1
            if first & self.size:
                first ^= self.modulus
        return result

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            value = self.mul(value, value)
            exponent >>= 1
        return result

    def inverse(self, value):
        if not value:
            raise ZeroDivisionError('A field inverse requires a nonzero element.')
        return self.power(value, self.size - 2)


def pfaffian(matrix, field):
    """Alternating-matrix Pfaffian by cubic-time paired elimination."""
    size = len(matrix)
    if size % 2:
        raise ValueError('The Pfaffian matrix must have even dimension.')
    values = [list(row) for row in matrix]
    if any(len(row) != size for row in values):
        raise ValueError('The Pfaffian matrix must be square.')
    if any(values[i][i] or any(values[i][j] != values[j][i] for j in range(i))
           for i in range(size)):
        raise ValueError('The matrix must be alternating in characteristic two.')
    result = 1
    for start in range(0, size, 2):
        pivot = next((j for j in range(start + 1, size) if values[start][j]), None)
        if pivot is None:
            return 0
        second = start + 1
        if pivot != second:
            values[second], values[pivot] = values[pivot], values[second]
            for row in values:
                row[second], row[pivot] = row[pivot], row[second]
        scalar = values[start][second]
        result = field.mul(result, scalar)
        reciprocal = field.inverse(scalar)
        for i in range(start + 2, size):
            for j in range(i + 1, size):
                correction = field.mul(values[start][i], values[second][j])
                correction ^= field.mul(values[second][i], values[start][j])
                values[i][j] ^= field.mul(reciprocal, correction)
                values[j][i] = values[i][j]
    return result


@dataclass(frozen=True)
class QueryResult:
    value: int
    forced_sites: int
    residual_forms: int
    zero_reason: str | None = None


def squarefree_coefficient_stream(forms, sites, field):
    """Return [z_1...z_k] product_i(m_i+b_i sum_a A_ai z_a).

    Each streamed site is (column_of_A, m_i, b_i). No complete input matrix
    is retained. The result includes small diagnostic metadata for controls.
    Field work is O(P*k^2+k^3), and field-element workspace is O(k^2),
    for k>=1. Addition is XOR in the supplied binary field representation.
    """
    if forms < 0:
        raise ValueError('The number of forms cannot be negative.')
    first = [0] * forms
    second = [[0] * forms for _ in range(forms)]
    forced = []
    scalar = 1
    for column, ordinary, inserted in sites:
        if len(column) != forms:
            raise ValueError('Each site column must contain one weight per form.')
        if not ordinary:
            if not inserted:
                return QueryResult(0, len(forced) + 1, max(0, forms - len(forced) - 1),
                                   'zero local factor')
            if len(forced) == forms:
                return QueryResult(0, len(forced) + 1, 0, 'too many forced sites')
            forced.append(list(column))
            scalar = field.mul(scalar, inserted)
            continue
        scalar = field.mul(scalar, ordinary)
        beta = field.mul(inserted, field.inverse(ordinary))
        weighted = [field.mul(beta, value) for value in column]
        for a, value in enumerate(weighted):
            first[a] ^= value
            for b in range(a, forms):
                second[a][b] ^= field.mul(value, weighted[b])
    for a in range(forms):
        for b in range(a):
            second[a][b] = second[b][a]

    count = len(forced)
    reduced = [[forced[j][i] for j in range(count)] for i in range(forms)]
    change = [[int(i == j) for j in range(forms)] for i in range(forms)]
    inverse_determinant = 1
    for col in range(count):
        pivot = next((i for i in range(col, forms) if reduced[i][col]), None)
        if pivot is None:
            return QueryResult(0, count, forms - count, 'dependent forced columns')
        reduced[col], reduced[pivot] = reduced[pivot], reduced[col]
        change[col], change[pivot] = change[pivot], change[col]
        value = reduced[col][col]
        inverse_determinant = field.mul(inverse_determinant, value)
        reciprocal = field.inverse(value)
        reduced[col] = [field.mul(reciprocal, x) for x in reduced[col]]
        change[col] = [field.mul(reciprocal, x) for x in change[col]]
        for i in range(forms):
            if i == col or not reduced[i][col]:
                continue
            multiple = reduced[i][col]
            reduced[i] = [x ^ field.mul(multiple, y)
                          for x, y in zip(reduced[i], reduced[col])]
            change[i] = [x ^ field.mul(multiple, y)
                         for x, y in zip(change[i], change[col])]

    projection = change[count:]
    residual = forms - count
    projected_first = []
    left_second = []
    for row in projection:
        value = 0
        product_row = [0] * forms
        for i, weight in enumerate(row):
            value ^= field.mul(weight, first[i])
            for j in range(forms):
                product_row[j] ^= field.mul(weight, second[i][j])
        projected_first.append(value)
        left_second.append(product_row)
    alternating = [[0] * residual for _ in range(residual)]
    for a in range(residual):
        for b in range(a + 1, residual):
            value = field.mul(projected_first[a], projected_first[b])
            for j in range(forms):
                value ^= field.mul(left_second[a][j], projection[b][j])
            alternating[a][b] = alternating[b][a] = value
    if residual % 2:
        for row, value in zip(alternating, projected_first):
            row.append(value)
        alternating.append(projected_first + [0])
    answer = field.mul(field.mul(scalar, inverse_determinant), pfaffian(alternating, field))
    return QueryResult(answer, count, residual)
