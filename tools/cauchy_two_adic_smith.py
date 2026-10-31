#!/usr/bin/env python3
"""Explicit two-adic channels for an odd-denominator Cauchy matrix."""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial, prod


def valuation2(value: int) -> int:
    """Two-adic valuation of a nonzero integer."""
    if value == 0:
        raise ValueError("zero has no finite valuation")
    value = abs(value)
    return (value & -value).bit_length() - 1


def smith_exponent(index: int) -> int:
    if type(index) is not int or index < 0:
        raise ValueError("index must be a nonnegative integer")
    return 4 * index - 2 * index.bit_count()


def validate(rows: int, columns: int, odd_shift: int) -> None:
    if type(rows) is not int or rows < 1:
        raise ValueError("rows must be a positive integer")
    if type(columns) is not int or columns < 1:
        raise ValueError("columns must be a positive integer")
    if type(odd_shift) is not int or odd_shift % 2 == 0:
        raise ValueError("shift must be an odd integer")


def cauchy_matrix(rows: int, columns: int, odd_shift: int):
    validate(rows, columns, odd_shift)
    return [[Fraction(1, odd_shift + 2 * (a - b))
             for b in range(columns)] for a in range(rows)]


def rational_ldu(rows: int, columns: int, odd_shift: int):
    """Return H=L diag(d) U over Z_(2), using closed formulas.

    L is rows by r and U is r by columns, where r=min(rows,columns).
    Their leading square blocks are unit triangular. Every denominator
    is odd, and the pivot valuation is smith_exponent(k).
    """
    validate(rows, columns, odd_shift)
    c = odd_shift
    rank = min(rows, columns)
    left = [[Fraction(0) for _ in range(rank)] for _ in range(rows)]
    right = [[Fraction(0) for _ in range(columns)] for _ in range(rank)]
    diagonal = []
    for k in range(rank):
        positive = prod(c + 2 * j for j in range(k + 1))
        negative = prod(c - 2 * j for j in range(k + 1))
        for a in range(k, rows):
            denominator = prod(c + 2 * (a - j) for j in range(k + 1))
            left[a][k] = Fraction(comb(a, k) * positive, denominator)
        for b in range(k, columns):
            denominator = prod(c + 2 * (j - b) for j in range(k + 1))
            right[k][b] = Fraction(comb(b, k) * negative, denominator)
        denominator = c * prod(c * c - 4 * j * j for j in range(1, k + 1))
        diagonal.append(Fraction((-4) ** k * factorial(k) ** 2, denominator))
    return left, diagonal, right


def residue(value: Fraction, precision: int) -> int:
    if type(precision) is not int or precision < 1:
        raise ValueError("precision must be a positive integer")
    if value.denominator % 2 == 0:
        raise ValueError("only odd denominators have residues here")
    modulus = 1 << precision
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def weighted_monomial_count(dimension: int, copies: int, shift: int, budget: int):
    """Exact count of all mode monomials of total weight at most budget."""
    if any(type(x) is not int for x in (dimension, copies, shift, budget)):
        raise ValueError("parameters must be integers")
    if dimension < 1 or copies < 1 or shift < 1 or budget < 0:
        raise ValueError("invalid dimension, copies, shift, or budget")
    coefficients = [0] * (budget + 1)
    coefficients[0] = 1
    for k in range(dimension):
        weight = shift + smith_exponent(k)
        for _ in range(copies):
            for total in range(weight, budget + 1):
                coefficients[total] += coefficients[total - weight]
    return sum(coefficients)


def count_nonzero_channels(dimension: int, precision: int, shift: int = 0) -> int:
    if type(dimension) is not int or dimension < 1:
        raise ValueError("dimension must be a positive integer")
    if type(precision) is not int or precision < 1:
        raise ValueError("precision must be a positive integer")
    if type(shift) is not int or shift < 0:
        raise ValueError("shift must be a nonnegative integer")
    return sum(shift + smith_exponent(k) < precision for k in range(dimension))
