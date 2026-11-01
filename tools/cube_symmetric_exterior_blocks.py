#!/usr/bin/env python3
"""Exact blocks for symmetric multiplication in a characteristic-two exterior algebra.

The algebra has commuting generators e_i with e_i**2=0. A polynomial
coefficient is encoded as a binary integer; its bit j is [x**j].
These one-variable profiles are specialization screens for the multiaxis
Cayley family, not certificates of favorable generic multivariable ranks.
"""
from __future__ import annotations

from collections import Counter
from math import comb, lcm


def carryless_product(first: int, second: int, mask: int) -> int:
    if first.bit_count() < second.bit_count():
        first, second = second, first
    out = 0
    while second:
        bit = second & -second
        out ^= first << (bit.bit_length() - 1)
        second ^= bit
    return out & mask


def unit_inverse(value: int, precision: int) -> int:
    if not value & 1:
        raise ValueError("Only a power-series unit has a unit inverse.")
    inverse, width = 1, 1
    while width < precision:
        width = min(2 * width, precision)
        square, pending = 0, inverse
        while pending:
            bit = pending & -pending
            square ^= 1 << (2 * (bit.bit_length() - 1))
            pending ^= bit
        inverse = carryless_product(value, square, (1 << width) - 1)
    return inverse


def series_smith_orders(matrix: list[list[int]], precision: int) -> list[int]:
    """Complete orders when the caller proves every order is below precision."""
    n = len(matrix)
    if not n or precision < 1 or any(len(row) != n for row in matrix):
        raise ValueError("Supply a nonempty square matrix and positive precision.")
    mask = (1 << precision) - 1
    work = [[value & mask for value in row] for row in matrix]
    orders = []
    for k in range(n):
        exponent, pivot_i, pivot_j = precision, -1, -1
        for i in range(k, n):
            for j in range(k, n):
                value = work[i][j]
                if value:
                    current = (value & -value).bit_length() - 1
                    if current < exponent:
                        exponent, pivot_i, pivot_j = current, i, j
        if pivot_i < 0:
            raise ArithmeticError("Singular matrix or insufficient absolute precision.")
        orders.append(exponent)
        work[k], work[pivot_i] = work[pivot_i], work[k]
        for row in work:
            row[k], row[pivot_j] = row[pivot_j], row[k]
        inverse = unit_inverse(work[k][k] >> exponent, precision - exponent)
        nonzero_columns = []
        for j in range(k + 1, n):
            if work[k][j]:
                work[k][j] = carryless_product(work[k][j], inverse, mask)
                nonzero_columns.append(j)
        for i in range(k + 1, n):
            factor = work[i][k] >> exponent
            if factor:
                for j in nonzero_columns:
                    work[i][j] ^= carryless_product(factor, work[k][j], mask)
    if orders != sorted(orders):
        raise ArithmeticError("The elimination did not return sorted orders.")
    return orders


def multiplication_matrix(coefficients: list[int]) -> list[list[int]]:
    """Coefficients are indexed by subsets of the square-zero generators."""
    n = len(coefficients)
    if n < 1 or n & (n - 1):
        raise ValueError("The coefficient array must have power-of-two length.")
    return [[coefficients[i ^ j] if i & j == j else 0 for j in range(n)]
            for i in range(n)]


def symmetric_blocks(coefficients: list[int]) -> list[tuple[int, list[int]]]:
    """Return (multiplicity, subset-coefficient array) for a_0+...+a_r e_r."""
    r = len(coefficients) - 1
    if r < 1:
        raise ValueError("Use at least one exterior generator.")
    blocks = []
    if r & 1:
        s = r // 2
        for k in range(s + 1):
            local = [coefficients[2 * (mask >> 1).bit_count() + (mask & 1)]
                     for mask in range(1 << (k + 1))]
            blocks.append((comb(s, k) * (1 << (s - k)), local))
    else:
        s = r // 2 - 1
        for k in range(s + 1):
            local = []
            for mask in range(1 << (k + 2)):
                low = mask & 3
                local.append(0 if low == 2 else
                             coefficients[2 * (mask >> 2).bit_count() + low.bit_count()])
            blocks.append((comb(s, k) * (1 << (s - k)), local))
    assert sum(mult * len(local) for mult, local in blocks) == 1 << r
    return blocks


def symmetric_profile(coefficients: list[int], precision: int) -> Counter:
    result = Counter()
    for multiplicity, local in symmetric_blocks(coefficients):
        for order in series_smith_orders(multiplication_matrix(local), precision):
            result[order] += multiplicity
    return result


def natural_cayley_profile(d: int) -> tuple[int, Counter]:
    if d < 2:
        raise ValueError("Use at least two cube axes.")
    degree = lcm(*range(1, d + 1))
    coefficients = [1 << (degree // (j + 1)) for j in range(d)]
    profile = symmetric_profile(coefficients, 2 * degree + 1)
    order = 1 << (d - 1)
    assert sum(profile.values()) == order
    assert sum(value * multiplicity for value, multiplicity in profile.items()) == order * degree
    assert profile == Counter({2 * degree - value: multiplicity
                               for value, multiplicity in profile.items()})
    return degree, profile


def exterior_product(first: int, second: int) -> int:
    """Exterior polynomials use bit position S for the squarefree monomial S."""
    result = 0
    while first:
        bit = first & -first
        left = bit.bit_length() - 1
        pending = second
        while pending:
            other = pending & -pending
            right = other.bit_length() - 1
            if not left & right:
                result ^= 1 << (left | right)
            pending ^= other
        first ^= bit
    return result


def linear_form(vector: int) -> int:
    result = 0
    while vector:
        bit = vector & -vector
        result ^= 1 << bit
        vector ^= bit
    return result


def _symplectic_pairs(even_dimension: int, first_all_ones: bool) -> list[tuple[int, int]]:
    if even_dimension & 1:
        raise ValueError("The symplectic dimension must be even.")
    if not even_dimension:
        return []
    if first_all_ones:
        pool = [(1 << even_dimension) - 1] + [1 << j for j in range(even_dimension - 1)]
    else:
        pool = [1 << j for j in range(even_dimension)]

    def inner(first, second):
        return ((first & second).bit_count() + first.bit_count() * second.bit_count()) & 1

    pairs = []
    while pool:
        first = pool.pop(0)
        index = next(j for j, value in enumerate(pool) if inner(first, value))
        second = pool.pop(index)
        pairs.append((first, second))
        pool = [value ^ (first if inner(value, second) else 0)
                      ^ (second if inner(value, first) else 0) for value in pool]
    return pairs


def canonical_vectors(r: int) -> tuple[int, list[tuple[int, int]]]:
    """u=sum e_i and pairs with sum v_i*w_i equal to the elementary 2-form."""
    if r < 1:
        raise ValueError("Use at least one generator.")
    u = (1 << r) - 1
    if r & 1:
        def into_even_hyperplane(value):
            return value | ((value.bit_count() & 1) << (r - 1))
        pairs = [(into_even_hyperplane(v), into_even_hyperplane(w))
                 for v, w in _symplectic_pairs(r - 1, False)]
    else:
        pairs = _symplectic_pairs(r, True)
        assert pairs[0][0] == u
    return u, pairs
