#!/usr/bin/env python3
"""Exact zero energies for arbitrary invertible binary linear exponent maps."""
from __future__ import annotations

from tensor_pascal_box_crossings import matrix_product, polynomial_power
from weighted_pascal_flags import weighted_kernel


def binary_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            leading = value.bit_length() - 1
            if leading in pivots:
                value ^= pivots[leading]
            else:
                pivots[leading] = value
                break
    return len(pivots)


def validate_linear_map(rows):
    rows = list(rows)
    dimension = len(rows)
    if (dimension < 1 or any(not isinstance(row, int) or not 0 <= row < 1 << dimension for row in rows)
            or binary_rank(rows) != dimension):
        raise ValueError("Use the rows of an invertible binary square matrix as integer masks.")
    return rows


def projected_basis_weight(rows, mask):
    """Minimum sum of 2^j over a row basis after projection onto mask's columns."""
    rows = validate_linear_map(rows)
    if not isinstance(mask, int) or not 0 <= mask < 1 << len(rows):
        raise ValueError("Use a column mask within the specified binary dimension.")
    pivots, weight = {}, 0
    for index, row in enumerate(rows):
        value = row & mask
        while value:
            leading = value.bit_length() - 1
            if leading in pivots:
                value ^= pivots[leading]
            else:
                pivots[leading] = value
                weight += 1 << index
                break
    if len(pivots) != mask.bit_count():
        raise ArithmeticError("An invertible matrix must have independent projected columns.")
    return weight


def linear_images(rows):
    rows = validate_linear_map(rows)
    return [sum(((word & row).bit_count() & 1) << j for j, row in enumerate(rows))
            for word in range(1 << len(rows))]


def linear_pascal_kernel(rows):
    rows = validate_linear_map(rows)
    dimension, length = len(rows), 1 << len(rows)
    if length > 128:
        raise ValueError("The finite-control dense implementation limits modes to128.")
    seed = weighted_kernel([1 << j for j in range(dimension)], length)
    images = linear_images(rows)
    exponents = [2 * value + 1 for value in images]
    basis_weights = [projected_basis_weight(rows, mask) for mask in range(length)]
    zero = [2 * length - 1 - 2 * basis_weights[length - 1 - i] for i in range(length)]
    pascal = seed["pascal"]
    transposed = [list(row) for row in zip(*pascal)]
    first = [[polynomial_power(3, exponents[i]) if transposed[i][j] else 0
              for j in range(length)] for i in range(length)]
    diagonal_pascal = [[polynomial_power(3, 2 * length - exponents[i]) if pascal[i][j] else 0
                       for j in range(length)] for i in range(length)]
    second = matrix_product(transposed, diagonal_pascal)
    original = seed["matrix"]
    seed["matrix"] = [first[i] + original[i][length:] for i in range(length)]
    seed["matrix"] += [original[length + i][:length] + second[i] for i in range(length)]
    seed.update({"linear_rows": rows, "images": images, "projected_basis_weights": basis_weights,
                 "exponents": exponents, "energies": exponents, "zero_energies": zero,
                 "infinity_spectrum_is_explicit": False})
    return seed
