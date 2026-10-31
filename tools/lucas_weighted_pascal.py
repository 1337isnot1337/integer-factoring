#!/usr/bin/env python3
"""Gap-skipping Pascal complements and weighted profiles in prime fields.

The field characteristic is a fixed prime. Exact profile calls should use
integer or Fraction heights. Float heights are only suitable for search.
"""

from __future__ import annotations

from math import comb, isqrt


def _check_prime(prime):
    if prime < 2 or any(prime % divisor == 0 for divisor in range(2, isqrt(prime) + 1)):
        raise ValueError('the characteristic must be prime')


def lucas_binomial(top, bottom, prime):
    """binom(top, bottom) modulo prime, without constructing the integer."""
    if prime < 2:
        raise ValueError('the characteristic must be at least two')
    if top < 0 or bottom < 0:
        raise ValueError('Pascal indices must be nonnegative')
    result = 1
    while top or bottom:
        a, b = top % prime, bottom % prime
        if b > a:
            return 0
        result = result * comb(a, b) % prime
        top //= prime
        bottom //= prime
    return result


class FieldBasis:
    def __init__(self, prime):
        self.prime = prime
        self.rows = {}

    def independent(self, vector, retain=False):
        value = list(vector)
        for pivot in range(len(value)):
            coefficient = value[pivot]
            if not coefficient:
                continue
            if pivot in self.rows:
                old = self.rows[pivot]
                value = [(a - coefficient * b) % self.prime
                         for a, b in zip(value, old)]
            else:
                if retain:
                    inverse = pow(coefficient, -1, self.prime)
                    self.rows[pivot] = [a * inverse % self.prime for a in value]
                return True
        return False


def greedy_complement(indices, prime, *, statistics=None):
    """Greedy outside Pascal columns, skipping complete radix-prime blocks."""
    _check_prime(prime)
    modes = tuple(sorted(indices))
    if not modes or modes[0] < 0 or len(set(modes)) != len(modes):
        raise ValueError('modes must be a nonempty set of nonnegative integers')
    excluded = set(modes)
    basis = FieldBasis(prime)
    result = []
    count = {'blocks_visited': 0, 'blocks_pruned': 0, 'leaves_visited': 0}

    def visit(start, power):
        if len(result) == len(modes):
            return
        count['blocks_visited'] += 1
        groups = {}
        for row, mode in enumerate(modes):
            coefficient = lucas_binomial(start // power, mode // power, prime)
            if not coefficient:
                continue
            residue = mode % power
            if residue not in groups:
                groups[residue] = [0] * len(modes)
            groups[residue][row] = coefficient
        if all(not basis.independent(vector) for vector in groups.values()):
            count['blocks_pruned'] += 1
            return
        if power == 1:
            count['leaves_visited'] += 1
            if start not in excluded:
                column = [lucas_binomial(start, mode, prime) for mode in modes]
                if basis.independent(column, retain=True):
                    result.append(start)
            return
        step = power // prime
        for digit in range(prime):
            visit(start + digit * step, step)

    power = 1
    while power <= modes[-1]:
        power *= prime
    # This range contains [0, 2 max(S)+1], which already spans all rows.
    visit(0, power * prime)
    if len(result) != len(modes):
        raise ArithmeticError('the complement did not attain full row rank')
    if statistics is not None:
        statistics.update(count)
    return tuple(result)


def weighted_constant_profile(matrix, row_weights, column_weights, prime):
    """Smith orders of diag(t^rows) matrix diag(t^columns).

    Complete minimum-weight pivoting preserves this monomial form.
    The returned matching identifies the row and column contributing each
    order, useful for exploratory piecewise derivatives.
    """
    _check_prime(prime)
    size = len(matrix)
    if size != len(row_weights) or size != len(column_weights) or any(len(row) != size for row in matrix):
        raise ValueError('a square matrix and matching weight lists are required')
    work = [[value % prime for value in row] for row in matrix]
    rows, columns = list(range(size)), list(range(size))
    orders, matching = [], []
    while rows:
        candidates = [(row_weights[i] + column_weights[j], i, j)
                      for i in rows for j in columns if work[i][j]]
        if not candidates:
            raise ArithmeticError('the constant matrix is singular')
        order, i, j = min(candidates)
        inverse = pow(work[i][j], -1, prime)
        orders.append(order)
        matching.append((i, j))
        rows.remove(i)
        columns.remove(j)
        for k in rows:
            multiplier = work[k][j] * inverse % prime
            if multiplier:
                for ell in columns:
                    work[k][ell] = (work[k][ell] - multiplier * work[i][ell]) % prime
    if orders != sorted(orders):
        raise ArithmeticError('minimum-weight elimination returned unsorted orders')
    return tuple(orders), tuple(matching)


def weighted_pascal_profile(indices, heights, prime):
    """Normalized infinity orders with equal added heights on both row blocks.

    Regular elimination requires h_a >= h_b whenever binom(b,a) is nonzero.
    Lucas antichains satisfy this requirement for every height vector.
    """
    _check_prime(prime)
    modes = tuple(indices)
    if modes != tuple(sorted(set(modes))) or len(modes) != len(heights):
        raise ValueError('modes must be increasing and have one height each')
    for i, a in enumerate(modes):
        for j, b in enumerate(modes):
            if lucas_binomial(b, a, prime) and heights[i] < heights[j]:
                raise ValueError('the height vector violates regular Pascal elimination')
    complement = greedy_complement(modes, prime)
    matrix = [[lucas_binomial(c, a, prime) for c in complement] for a in modes]
    orders, matching = weighted_constant_profile(matrix, heights, complement, prime)
    combined = sorted([a + h for a, h in zip(modes, heights)] + list(orders))
    return tuple(combined), complement, matching
