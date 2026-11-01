#!/usr/bin/env python3
"""Recover an integer factor polynomial from supplied regularized moments.

This module does not produce the displaced Bezout jet or factor moments.
Its input contract is the unit-weight moment representation proved in
research/regularized-bezout-factor-moments.md.
"""

from __future__ import annotations


def extraction_precision(factor_bound: int, max_factor: int,
                         separation_bits: int | None = None) -> dict:
    if any(not isinstance(value, int) or value < 1
           for value in (factor_bound, max_factor)):
        raise ValueError("Require a positive factor bound and positive node bound")
    if separation_bits is None:
        separation_bits = max(1, max_factor.bit_length())
    if not isinstance(separation_bits, int) or separation_bits < 1:
        raise ValueError("Require a positive separation valuation bound")
    coefficient_bits = factor_bound * (max_factor + 1).bit_length() + 2
    determinant_guard = factor_bound * (factor_bound - 1) * separation_bits
    return {"factor_bound": factor_bound, "max_factor": max_factor,
            "separation_bits": separation_bits,
            "coefficient_bits": coefficient_bits,
            "determinant_guard_bits": determinant_guard,
            "precision": coefficient_bits + determinant_guard + 1,
            "moment_count": 2 * factor_bound,
            "decision_precision": (factor_bound - 1) * separation_bits + 1,
            "decision_moment_count": factor_bound}


def determinant_mod(matrix: list[list[int]], modulus: int) -> int:
    """Fraction-free determinant of integer lifts, reduced at the end."""
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("The determinant input must be square")
    if size == 0:
        return 1 % modulus
    values = [[entry % modulus for entry in row] for row in matrix]
    previous, sign = 1, 1
    for column in range(size - 1):
        pivot_row = next((row for row in range(column, size)
                          if values[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            values[column], values[pivot_row] = values[pivot_row], values[column]
            sign = -sign
        pivot = values[column][column]
        for row in range(column + 1, size):
            for col in range(column + 1, size):
                numerator = (pivot * values[row][col]
                             - values[row][column] * values[column][col])
                if numerator % previous:
                    raise ArithmeticError("Fraction-free determinant division failed")
                values[row][col] = numerator // previous
            values[row][column] = 0
        previous = pivot
    return sign * values[-1][-1] % modulus


def decode_factor_moments(moments: list[int] | tuple[int, ...], *,
                          factor_bound: int, max_factor: int,
                          precision: int,
                          separation_bits: int | None = None) -> dict:
    """Decode sum_a unit_weight[a]*a^j modulo 2^precision.

    The distinct positive integer nodes must be <=max_factor, with at most
    factor_bound nodes. Pairwise difference valuations must be bounded by
    separation_bits. The weights have valuation zero. The routine returns
    the exact monic integer polynomial with those nodes as roots.
    """
    budget = extraction_precision(factor_bound, max_factor, separation_bits)
    if not isinstance(precision, int) or precision < budget["precision"]:
        raise ValueError("Insufficient precision for determinant and coefficient recovery")
    if len(moments) < budget["moment_count"]:
        raise ValueError("Require twice the declared factor bound in supplied moments")
    if any(not isinstance(value, int) for value in moments):
        raise ValueError("Moment residues must be integers")
    modulus = 1 << precision
    data = [value % modulus for value in moments[:budget["moment_count"]]]
    determinants = []
    degree = 0
    for size in range(1, factor_bound + 1):
        matrix = [[data[i + j] for j in range(size)] for i in range(size)]
        determinant = determinant_mod(matrix, modulus)
        determinants.append(determinant)
        if determinant:
            degree = size
    if degree == 0:
        if any(data):
            raise ValueError("Supplied moments violate the declared unit-weight model")
        return {**budget, "supplied_precision": precision, "degree": 0,
                "factor_polynomial": [1], "hankel_determinant_valuation": None,
                "nonzero_leading_hankel_sizes": []}
    determinant = determinants[degree - 1]
    determinant_order = (determinant & -determinant).bit_length() - 1
    if determinant_order > budget["determinant_guard_bits"]:
        raise ValueError("The moment determinant exceeds the declared conditioning bound")
    bits = budget["coefficient_bits"]
    coefficient_modulus = 1 << bits
    unit_inverse = pow((determinant >> determinant_order) % coefficient_modulus,
                       -1, coefficient_modulus)
    matrix = [[data[i + j] for j in range(degree)] for i in range(degree)]
    right = [-data[degree + i] for i in range(degree)]
    polynomial = []
    for column in range(degree):
        replaced = [row[:] for row in matrix]
        for row in range(degree):
            replaced[row][column] = right[row]
        numerator = determinant_mod(replaced, modulus)
        if numerator % (1 << determinant_order):
            raise ValueError("The supplied recurrence has nonintegral coefficients")
        residue = ((numerator >> determinant_order) * unit_inverse) % coefficient_modulus
        if residue >= coefficient_modulus // 2:
            residue -= coefficient_modulus
        polynomial.append(residue)
    polynomial.append(1)
    for offset in range(len(data) - degree):
        if sum(coefficient * data[offset + i]
               for i, coefficient in enumerate(polynomial)) % modulus:
            raise ValueError("Recovered polynomial does not satisfy the supplied moments")
    if any(abs(value) > (max_factor + 1) ** factor_bound for value in polynomial):
        raise ValueError("Recovered coefficients violate the declared integer height bound")
    return {**budget, "supplied_precision": precision, "degree": degree,
            "factor_polynomial": polynomial,
            "hankel_determinant_valuation": determinant_order,
            "nonzero_leading_hankel_sizes": [i + 1 for i, value in enumerate(determinants) if value]}
