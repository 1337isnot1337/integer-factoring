#!/usr/bin/env python3
"""Exact polynomial kernels, physical restrictions, and profiles for Pascal boxes.

Binary integers encode F2[u] polynomials. The full-Fock helper is deliberately
limited to four modes; the one-particle constructor also handles larger boxes.
The theorem for arbitrary fixed boxes is in tensor-pascal-box-crossings.md.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb

from cube_symmetric_exterior_blocks import series_smith_orders


def polynomial_product(first: int, second: int) -> int:
    if first.bit_count() < second.bit_count():
        first, second = second, first
    out = 0
    while second:
        bit = second & -second
        out ^= first << (bit.bit_length() - 1)
        second ^= bit
    return out


def polynomial_power(value: int, exponent: int) -> int:
    if not isinstance(exponent, int) or exponent < 0:
        raise ValueError("A polynomial exponent must be a nonnegative integer.")
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, value)
        value = polynomial_product(value, value)
        exponent >>= 1
    return result


def reverse_polynomial(value: int) -> int:
    return int(bin(value)[2:][::-1], 2) if value else 0


def matrix_product(first, second):
    columns = list(zip(*second))
    result = []
    for row in first:
        out = []
        for column in columns:
            value = 0
            for a, b in zip(row, column):
                value ^= polynomial_product(a, b)
            out.append(value)
        result.append(out)
    return result


def box_kernel(coordinates: int, side: int):
    if any(not isinstance(v, int) or v < 1 or v & (v - 1)
           for v in (coordinates, side)):
        raise ValueError("The coordinate count and side must be positive powers of two.")
    n = side ** coordinates
    if n > 128:
        raise ValueError("This finite-control constructor is limited to 128 modes.")
    modes = list(product(range(side), repeat=coordinates))
    outer = coordinates * side
    energies = [2 * sum(a) + coordinates for a in modes]
    pascal = [[int(all(y & x == y for x, y in zip(a, b))) for b in modes]
              for a in modes]
    matrix = [[0] * (2 * n) for _ in range(2 * n)]
    for a, left in enumerate(modes):
        for b, right in enumerate(modes):
            matrix[a][b] = polynomial_power(3, energies[a]) if pascal[b][a] else 0
            matrix[a][b + n] = (1 << (2 * outer)) if pascal[a][b] else 0
            matrix[a + n][b] = int(a == b)
            shift = 2 * (coordinates * (side - 1) - sum(x | y for x, y in zip(left, right)))
            matrix[a + n][b + n] = polynomial_power(3, coordinates) << shift
    return {"coordinates": coordinates, "side": side, "modes": modes,
            "n": n, "outer": outer, "energies": energies,
            "pascal": pascal, "matrix": matrix}


def mode_heights(seed, center=None):
    outer = seed["outer"]
    center = outer if center is None else center
    return [((0, energy - 2 * outer, 0, energy),
             (energy - center, 2 * outer - center, energy - center, -center))
            for energy in seed["energies"]]


def physical_heights(seed, symbol, center=None):
    heights = mode_heights(seed, center)
    return tuple(tuple(sum(heights[j][place][leg] for j in range(seed["n"])
                           if symbol >> j & 1) for leg in range(4))
                 for place in range(2))


def weighted_infinity_matrix(seed, center=None):
    """Return a polynomial shift of the actual gauged one-particle matrix."""
    n = seed["n"]
    heights = mode_heights(seed, center)
    entries, minimum = [], 0
    for i, row in enumerate(seed["matrix"]):
        for j, value in enumerate(row):
            if not value:
                continue
            row_mode, row_leg = (i, 2) if i < n else (i - n, 3)
            column_mode, column_leg = (j, 0) if j < n else (j - n, 1)
            shift = (heights[row_mode][1][row_leg]
                     + heights[column_mode][1][column_leg] - (value.bit_length() - 1))
            reversed_value = reverse_polynomial(value)
            valuation = (reversed_value & -reversed_value).bit_length() - 1
            minimum = min(minimum, valuation + shift)
            entries.append((i, j, reversed_value, shift))
    matrix = [[0] * (2 * n) for _ in range(2 * n)]
    for i, j, value, shift in entries:
        adjusted = shift - minimum
        if adjusted >= 0:
            matrix[i][j] = value << adjusted
        else:
            assert not value & ((1 << -adjusted) - 1)
            matrix[i][j] = value >> -adjusted
    return matrix, minimum


def exterior_tensor(seed):
    if seed["n"] > 4:
        raise ValueError("The full physical tensor helper is limited to four modes.")
    n, matrix = seed["n"], seed["matrix"]

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        pending, value = columns, 0
        while pending:
            other = pending & -pending
            j = other.bit_length() - 1
            if matrix[i][j]:
                value ^= polynomial_product(matrix[i][j], minor(rows ^ bit, columns ^ other))
            pending ^= other
        return value

    q = 1 << n
    values = {(b, x, v, y): minor(v | (y << n), b | (x << n))
              for b, x, v, y in product(range(q), repeat=4)}
    return {"values": values, "alphabet": list(range(q))}


def retract(tensor, left, right, alphabet):
    """Maps on (B,X,V,Y) are (left,left,right,right), stored as rows."""
    dimension = len(left)
    if len(right) != dimension or len(alphabet) != dimension:
        raise ValueError("The physical maps and height labels need matching dimensions.")
    ambient = len(tensor["alphabet"])
    if any(len(row) != ambient or any(value not in (0, 1) for value in row)
           for row in left + right):
        raise ValueError("Use binary maps with one column for each ambient physical state.")
    if any(sum(a * b for a, b in zip(l, r)) % 2 != int(i == j)
           for i, l in enumerate(left) for j, r in enumerate(right)):
        raise ValueError("The binary physical maps must be a split projection.")
    supports_l = [[i for i, value in enumerate(row) if value] for row in left]
    supports_r = [[i for i, value in enumerate(row) if value] for row in right]
    values = {}
    for key in product(range(dimension), repeat=4):
        value = 0
        for inputs in product(supports_l[key[0]], supports_l[key[1]],
                              supports_r[key[2]], supports_r[key[3]]):
            value ^= tensor["values"][inputs]
        values[key] = value
    return {"values": values, "alphabet": list(alphabet)}


def physical_profiles(seed, tensor, center=None):
    """Exact complete profiles at the declared integer inherited heights.

    The typed endpoint proof supplies full rank for the permitted homogeneous
    split maps. A complete list of pivots below the absolute window certifies
    each returned valuation. Missing pivots increase the absolute precision.
    """
    alphabet = tensor["alphabet"]
    q = len(alphabet)
    heights = [physical_heights(seed, symbol, center) for symbol in alphabet]
    result, precisions = [], []
    for place in range(2):
        minimum, maximum, entries = 0, 0, []
        for key, value in tensor["values"].items():
            if not value:
                continue
            shift = sum(heights[key[leg]][place][leg] for leg in range(4))
            if place:
                shift -= value.bit_length() - 1
                value = reverse_polynomial(value)
            valuation = (value & -value).bit_length() - 1
            minimum = min(minimum, valuation + shift)
            maximum = max(maximum, value.bit_length() - 1 + shift)
            entries.append((key, value, shift))
        matrix = [[0] * (q * q) for _ in range(q * q)]
        for (b, x, v, y), value, shift in entries:
            row, column = (b + q * v, x + q * y) if not place else (b + q * x, v + q * y)
            adjusted = shift - minimum
            if adjusted >= 0:
                value <<= adjusted
            else:
                assert not value & ((1 << -adjusted) - 1)
                value >>= -adjusted
            matrix[row][column] = value
        bound = q * q * (maximum - minimum)
        precision = min(bound + 1, max(64, 2 * (maximum - minimum) + 1))
        while True:
            try:
                orders = series_smith_orders(matrix, precision)
                break
            except ArithmeticError:
                if precision >= bound + 1:
                    raise ArithmeticError("The supplied physical flattening is singular.")
                precision = min(2 * precision, bound + 1)
        result.append([order + minimum for order in orders])
        precisions.append(precision)
    total_heights = [sum(map(sum, item)) for item in heights]
    return {"profiles": result, "total_heights": total_heights,
            "absolute_precisions": precisions}


ENERGY_PAIRS = ((2, 4), (3, 5), (6, 9), (10, 12), (11, 13))
SPLIT_CHOICES = tuple((left, right)
                     for left in ((1, 0), (0, 1), (1, 1))
                     for right in ((1, 0), (0, 1), (1, 1))
                     if sum(a * b for a, b in zip(left, right)) % 2 == 1)


def small_box_code(tensor, choices):
    omitted = {ENERGY_PAIRS[index][1] for index in choices}
    alphabet = [symbol for symbol in range(16) if symbol not in omitted]
    left, right = [], []
    for symbol in alphabet:
        first, second = [0] * 16, [0] * 16
        first[symbol] = second[symbol] = 1
        for index, choice in choices.items():
            a, b = ENERGY_PAIRS[index]
            if symbol == a:
                (first[a], first[b]), (second[a], second[b]) = SPLIT_CHOICES[choice]
        left.append(first)
        right.append(second)
    return retract(tensor, left, right, alphabet)


def three_quarter_certificate(profiles, heights):
    """P(1+v) has nonnegative integer coefficients, proving the all-scale floor."""
    def multiply(first, second):
        result = Counter()
        for a, x in first.items():
            for b, y in second.items():
                result[a + b] += x * y
        return result

    first = Counter(-3 * value for value in profiles[0])
    second = Counter(-3 * value for value in profiles[1])
    physical = Counter(-value for value in heights)
    numerator = multiply(first, second)
    denominator = multiply(multiply(physical, physical), physical)
    for exponent, coefficient in denominator.items():
        numerator[exponent] -= coefficient
    minimum = min(exponent for exponent, coefficient in numerator.items() if coefficient)
    polynomial = {exponent - minimum: coefficient for exponent, coefficient in numerator.items()
                  if coefficient}
    translated = [sum(coefficient * comb(exponent, index)
                      for exponent, coefficient in polynomial.items() if exponent >= index)
                  for index in range(max(polynomial) + 1)]
    if min(translated) < 0 or not any(translated):
        raise ArithmeticError("This polynomial did not prove the three-quarter floor.")
    return {"laurent_shift": -minimum, "degree": max(polynomial),
            "translated_coefficients": translated}
