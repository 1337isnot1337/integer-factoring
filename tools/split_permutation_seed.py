#!/usr/bin/env python3
"""Exact split-variable permutation crossing and its finite curve seeds."""

from itertools import permutations, product

from fermionic_piece_polynomial import polynomial_product


def split_polynomial(a, b, u, y):
    """Multilinear polynomial in P_i=X_i^2 and Q_j=Y_j^2.

    Coefficient positions are subset masks on the 2*n variables.
    The arguments are permutations, ordered as operator A, operator B,
    output U, input Y.
    """
    size = len(a)
    if not len(b) == len(u) == len(y) == size:
        raise ValueError("All four labels must have equal size.")
    edges = [frozenset(enumerate(label)) for label in (a, b, u, y)]
    if edges[0] ^ edges[1] ^ edges[2] ^ edges[3]:
        return 0
    active = {i for i in range(size) if a[i] != y[i]}
    result = 1 << sum(1 << i for i in active)
    inverse_y = {value: i for i, value in enumerate(y)}
    remaining = set(active)
    while remaining:
        start = next(iter(remaining))
        cycle = []
        current = start
        while current not in cycle:
            cycle.append(current)
            remaining.remove(current)
            current = inverse_y[a[current]]
        assert current == start
        if all(b[i] == y[i] for i in cycle):
            for i in cycle:
                q = 1 << (size + a[i])
                result = polynomial_product(result, 1 ^ (1 << q))
        else:
            assert all(b[i] == a[i] for i in cycle)
    for i in range(size):
        if i in active:
            continue
        q = 1 << (size + a[i])
        factor = 1 ^ (1 << q)
        if a[i] == b[i]:
            factor ^= 1 << (q + (1 << i))
        result = polynomial_product(result, factor)
    return result


def coefficient_tensors(size):
    """Return (labels, bit-packed row flattenings for every subset degree)."""
    labels = tuple(permutations(range(size)))
    count = len(labels)
    matrices = [0] * (1 << (2 * size))
    for ai, bi, ui, yi in product(range(count), repeat=4):
        polynomial = split_polynomial(labels[ai], labels[bi], labels[ui], labels[yi])
        entry = (ai * count + yi) * count**2 + bi * count + ui
        while polynomial:
            bit = polynomial & -polynomial
            polynomial ^= bit
            matrices[bit.bit_length() - 1] ^= 1 << entry
    return labels, tuple(matrices)


def unpack_matrix(packed, dimension):
    mask = (1 << dimension) - 1
    return tuple((packed >> (dimension * row)) & mask for row in range(dimension))


def column_matrix(rows, labels):
    result = [0] * (labels * labels)
    for row, value in enumerate(rows):
        a, y = divmod(row, labels)
        while value:
            bit = value & -value
            value ^= bit
            b, u = divmod(bit.bit_length() - 1, labels)
            result[a * labels + u] ^= 1 << (b * labels + y)
    return tuple(result)


def curve_matrices(tensors, variable_polynomials, labels):
    """Substitute F2[z] bitsets for each P/Q variable.

    Every variable polynomial must evaluate to1 at z=1.
    Common zero powers and unused trailing coefficients are removed.
    """
    if len(tensors) != 1 << len(variable_polynomials):
        raise ValueError("The number of variables does not match the tensor family.")
    if any(polynomial.bit_count() % 2 != 1 for polynomial in variable_polynomials):
        raise ValueError("Every curve coordinate must equal1 at z=1.")
    coefficients = []
    for subset, tensor in enumerate(tensors):
        if not tensor:
            continue
        polynomial = 1
        for variable, value in enumerate(variable_polynomials):
            if subset >> variable & 1:
                polynomial = polynomial_product(polynomial, value)
        while polynomial:
            bit = polynomial & -polynomial
            polynomial ^= bit
            degree = bit.bit_length() - 1
            if len(coefficients) <= degree:
                coefficients.extend([0] * (degree + 1 - len(coefficients)))
            coefficients[degree] ^= tensor
    while coefficients and not coefficients[-1]:
        coefficients.pop()
    while coefficients and not coefficients[0]:
        coefficients.pop(0)
    if not coefficients:
        raise AssertionError("A crossing curve cannot vanish identically.")
    dimension = labels * labels
    rows = [unpack_matrix(value, dimension) for value in coefficients]
    columns = [column_matrix(matrix, labels) for matrix in rows]
    return rows, columns
