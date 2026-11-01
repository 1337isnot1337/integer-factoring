#!/usr/bin/env python3
"""Pascal crossings from all Boolean maps preserving the binary power flag.

The exact zero-cut profile is retained. Infinity profiles must be computed
from the actual kernel; the exponent list is not an infinity Smith list.
"""
from __future__ import annotations

from tensor_pascal_box_crossings import matrix_product, polynomial_power


def triangular_map(diagonal, tables):
    """Return phi_j(x)=diagonal[j]*x_j + tables[j](x_0,...,x_{j-1})."""
    diagonal, tables = list(diagonal), [list(row) for row in tables]
    if (not diagonal or len(tables) != len(diagonal)
            or any(not isinstance(value, int) or value not in (0, 1) for value in diagonal)
            or any(len(row) != 1 << j or any(not isinstance(value, int) or value not in (0, 1) for value in row)
                   for j, row in enumerate(tables))):
        raise ValueError("Use one binary diagonal bit and one lower-bit truth table per coordinate.")
    return [sum(((diagonal[j] * ((i >> j) & 1)) ^ tables[j][i & ((1 << j) - 1)]) << j
                for j in range(len(diagonal))) for i in range(1 << len(diagonal))]


def decode_triangular_map(images):
    images = list(images)
    length = len(images)
    if (length < 2 or length & (length - 1)
            or any(not isinstance(value, int) or not 0 <= value < length for value in images)):
        raise ValueError("Use a map of a binary cube of size at least two to itself.")
    dimension = length.bit_length() - 1
    diagonal = [((images[1 << j] ^ images[0]) >> j) & 1 for j in range(dimension)]
    tables = [[(images[i] >> j) & 1 for i in range(1 << j)] for j in range(dimension)]
    if triangular_map(diagonal, tables) != images:
        raise ValueError("The map does not preserve the complete binary power filtration.")
    return diagonal, tables


def upset_orders(images):
    """Literal orders of sum_{i containing q}(1+u)^images[i]; None means zero."""
    images = list(images)
    powers = [polynomial_power(3, value) for value in images]
    result = []
    for q in range(len(images)):
        value = 0
        for i, power in enumerate(powers):
            if i & q == q:
                value ^= power
        result.append((value & -value).bit_length() - 1 if value else None)
    return result


def filtered_pascal_kernel(images):
    """Construct an actual polynomial kernel from a filtration-preserving map.

    Dense control matrices are limited to128 modes. The algebraic result
    applies to every power-of-two number of modes.
    """
    images = list(images)
    diagonal, tables = decode_triangular_map(images)
    length = len(images)
    if length > 128:
        raise ValueError("The finite-control implementation limits the dense mode count to128.")
    exponents = [2 * value + 1 for value in images]
    pascal = [[int(i & j == j) for j in range(length)] for i in range(length)]
    transposed = [list(row) for row in zip(*pascal)]
    first = [[polynomial_power(3, exponents[i]) if transposed[i][j] else 0
              for j in range(length)] for i in range(length)]
    diagonal_pascal = [[polynomial_power(3, 2 * length - exponents[i]) if pascal[i][j] else 0
                       for j in range(length)] for i in range(length)]
    second = matrix_product(transposed, diagonal_pascal)
    identity = [[int(i == j) for j in range(length)] for i in range(length)]
    upper_right = [[value << (2 * length) for value in row] for row in pascal]
    matrix = [first[i] + upper_right[i] for i in range(length)]
    matrix += [identity[i] + second[i] for i in range(length)]
    return {"n": length, "outer": length, "images": images,
            "diagonal": diagonal, "lower_bit_tables": tables,
            "exponents": exponents, "energies": exponents,
            "zero_energies": list(range(1, 2 * length, 2)),
            "pascal": pascal, "matrix": matrix, "flag_changed": False,
            "infinity_spectrum_is_explicit": False}
