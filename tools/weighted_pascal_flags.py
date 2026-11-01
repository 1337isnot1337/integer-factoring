#!/usr/bin/env python3
"""Exact weighted Pascal kernels and their two small Clifford matrices.

This finite-control implementation limits the dense mode count to128 and
the outer scale to1024. The algebraic constructor applies to every fixed
positive weight list whose sum is less than its power-of-two outer scale.
"""
from __future__ import annotations

from tensor_pascal_box_crossings import matrix_product, polynomial_power, polynomial_product


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def weighted_kernel(weights, outer):
    weights = list(weights)
    if (not weights or any(not isinstance(value, int) or value < 1 for value in weights)
            or not isinstance(outer, int) or outer < 1 or outer & (outer - 1)
            or sum(weights) >= outer):
        raise ValueError("Use positive integer weights with sum below a power-of-two outer scale.")
    n = 1 << len(weights)
    if n > 128 or outer > 1024:
        raise ValueError("The finite-control dense implementation limits modes to128 and scale to1024.")
    low_bits = [value & -value for value in weights]
    offset = outer - sum(weights)
    zero_offset = outer - sum(low_bits)
    energies = [offset + 2 * sum(weight for j, weight in enumerate(weights) if a >> j & 1)
                for a in range(n)]
    zero_energies = [zero_offset + 2 * sum(weight for j, weight in enumerate(low_bits) if a >> j & 1)
                     for a in range(n)]
    pascal = [[int(a & b == b) for b in range(n)] for a in range(n)]
    identity = [[int(a == b) for b in range(n)] for a in range(n)]
    first = [[polynomial_power(3, energies[a]) if pascal[b][a] else 0 for b in range(n)]
             for a in range(n)]
    second = []
    for a in range(n):
        row = []
        for b in range(n):
            value = polynomial_power(3, offset)
            for j, weight in enumerate(weights):
                if not (a | b) >> j & 1:
                    value = polynomial_product(value, polynomial_power(3, 2 * weight) ^ 1)
            row.append(value)
        second.append(row)
    upper_right = [[value << (2 * outer) for value in row] for row in pascal]
    matrix = [first[i] + upper_right[i] for i in range(n)]
    matrix += [identity[i] + second[i] for i in range(n)]
    return {"weights": weights, "low_bits": low_bits, "n": n, "outer": outer,
            "energies": energies, "zero_energies": zero_energies,
            "pascal": pascal, "matrix": matrix, "flag_changed": False}


def natural_heights(seed, center=None):
    outer = seed["outer"]
    center = outer if center is None else center
    zero, infinity = [], []
    for energy_zero, energy_infinity in zip(seed["zero_energies"], seed["energies"]):
        zero.append((0, energy_zero - 2 * outer, 0, energy_zero))
        infinity.append((energy_infinity - center, 2 * outer - center,
                         energy_infinity - center, -center))
    return [[heights[a][leg] for leg in range(4) for a in range(seed["n"])]
            for heights in (zero, infinity)]


def inverse_kernel(seed):
    n, outer = seed["n"], seed["outer"]
    matrix, pascal = seed["matrix"], seed["pascal"]
    first = [row[:n] for row in matrix[:n]]
    second = [row[n:] for row in matrix[n:]]
    pa, bp = matrix_product(pascal, first), matrix_product(second, pascal)
    scaled_identity = [[int(i == j) << (2 * outer) for j in range(n)] for i in range(n)]
    return [bp[i] + scaled_identity[i] for i in range(n)] + [pascal[i] + pa[i] for i in range(n)]


def clifford_data(seed, place):
    """Binary polynomial entries, Laurent shifts, and physical index maps.

    The variable is u at both places; reverse it for the infinity valuation.
    Physical modes use boundary order (B,X,V,Y). The matrix order is4n.
    """
    n, outer = seed["n"], seed["outer"]
    matrix, pascal = seed["matrix"], seed["pascal"]
    if place == 0:
        first = [row[:n] for row in matrix[:n]]
        second = [row[n:] for row in matrix[n:]]
        identity = [[int(i == j) for j in range(n)] for i in range(n)]
        zero = [[0] * n for _ in range(n)]
        pa, bp = matrix_product(pascal, first), matrix_product(second, pascal)
        blocks = [[zero, identity, transpose(pa), zero],
                  [pascal, zero, zero, first],
                  [second, zero, zero, identity],
                  [zero, transpose(bp), transpose(pascal), zero]]
        shifts = [[0, -2 * outer, -2 * outer, 0], [0, 0, 0, 0],
                  [0, 0, 0, 0], [0, -2 * outer, -2 * outer, 0]]
        values = [[blocks[bi][bj][i][j] for bj in range(4) for j in range(n)]
                  for bi in range(4) for i in range(n)]
        shifts = [[shifts[bi][bj] for bj in range(4) for _ in range(n)]
                  for bi in range(4) for _ in range(n)]
        rows = list(range(n)) + list(range(2 * n, 3 * n))
        columns = list(range(n, 2 * n)) + list(range(3 * n, 4 * n))
    elif place == 1:
        transposed, inverse = transpose(matrix), inverse_kernel(seed)
        zero = [[0] * (2 * n) for _ in range(2 * n)]
        values = [transposed[i] + zero[i] for i in range(2 * n)]
        values += [zero[i] + inverse[i] for i in range(2 * n)]
        shifts = [[0] * (4 * n) for _ in range(4 * n)]
        rows, columns = list(range(2 * n)), list(range(2 * n, 4 * n))
    else:
        raise ValueError("Choose place zero or one.")
    return values, shifts, rows, columns


def binary_inverse(matrix):
    n = len(matrix)
    if not n or any(len(row) != n or any(value not in (0, 1) for value in row) for row in matrix):
        raise ValueError("Use a nonempty square binary matrix.")
    rows = [sum(value << j for j, value in enumerate(row)) | (1 << (n + i))
            for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i] >> j & 1), None)
        if pivot is None:
            raise ValueError("The constant flag must be invertible.")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(n):
            if i != j and rows[i] >> j & 1:
                rows[i] ^= rows[j]
    return [[rows[i] >> (n + j) & 1 for j in range(n)] for i in range(n)]


def change_flag(seed, constant):
    """A -> A C^-1, B -> C B; retain the endpoint and Clifford formulas.

    The displayed natural spectra need not survive this operation. Always
    recompute both profiles for a changed flag and any proposed gauges.
    """
    n, matrix = seed["n"], seed["matrix"]
    if len(constant) != n:
        raise ValueError("Supply one flag coordinate per mode.")
    inverse = binary_inverse(constant)
    first = matrix_product([row[:n] for row in matrix[:n]], inverse)
    second = matrix_product(constant, [row[n:] for row in matrix[n:]])
    changed = dict(seed)
    changed["matrix"] = [first[i] + matrix[i][n:] for i in range(n)]
    changed["matrix"] += [matrix[n + i][:n] + second[i] for i in range(n)]
    changed["flag_changed"] = True
    changed["constant_flag"] = [list(row) for row in constant]
    return changed
