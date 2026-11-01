#!/usr/bin/env python3
"""Actual Pascal crossings for every exponent permutation, with explicit cuts."""
from __future__ import annotations

from tensor_pascal_box_crossings import matrix_product, polynomial_power
from weighted_pascal_flags import weighted_kernel


def validate_permutation(images):
    images = list(images)
    length = len(images)
    if (length < 2 or length & (length - 1)
            or any(type(value) is not int for value in images)
            or sorted(images) != list(range(length))):
        raise ValueError("Supply a permutation of a power-of-two interval of size at least two.")
    return images


def permutation_upset_orders(images):
    """Compute all orders by the Boolean zeta transform of actual polynomials."""
    images = validate_permutation(images)
    values = [polynomial_power(3, value) for value in images]
    for bit in range(len(images).bit_length() - 1):
        for mask in range(len(images)):
            if not mask >> bit & 1:
                values[mask] ^= values[mask | (1 << bit)]
    if any(value == 0 for value in values):
        raise ArithmeticError("A nonempty subset of distinct monomials cannot vanish.")
    return [(value & -value).bit_length() - 1 for value in values]


def permuted_pascal_kernel(images):
    images = validate_permutation(images)
    length = len(images)
    if length > 128:
        raise ValueError("The finite-control dense implementation limits modes to 128.")
    dimension = length.bit_length() - 1
    seed = weighted_kernel([1 << j for j in range(dimension)], length)
    orders = permutation_upset_orders(images)
    exponents = [2 * value + 1 for value in images]
    pascal = seed['pascal']
    transposed = [list(row) for row in zip(*pascal)]
    first = [[polynomial_power(3, exponents[i]) if transposed[i][j] else 0
              for j in range(length)] for i in range(length)]
    diagonal_pascal = [[polynomial_power(3, 2 * length - exponents[i]) if pascal[i][j] else 0
                       for j in range(length)] for i in range(length)]
    second = matrix_product(transposed, diagonal_pascal)
    original = seed['matrix']
    seed['matrix'] = [first[i] + original[i][length:] for i in range(length)]
    seed['matrix'] += [original[length + i][:length] + second[i] for i in range(length)]
    seed.update({'images': images, 'upset_orders': orders, 'exponents': exponents,
                 'energies': exponents, 'zero_energies': [2 * length - 1 - 2 * value for value in orders],
                 'infinity_spectrum_is_explicit': False})
    return seed


def explicit_infinity_gauge(seed, beta):
    """Boolean-isotonic beta gives the proved two-by-two infinity reduction."""
    beta = list(beta)
    length = seed['n']
    if (len(beta) != length or any(type(value) is not int for value in beta)
            or any(beta[i] > beta[j] for i in range(length) for j in range(length) if i & j == i)):
        raise ValueError("Use integer beta values increasing under Boolean containment.")
    exponents = seed['exponents']
    heights = ([value - length for value in beta] + [length] * length
               + [value - length for value in exponents] + [-length] * length)
    orders = sorted([min(a, b) - 2 * length for a, b in zip(beta, exponents)]
                    + [max(a, b) for a, b in zip(beta, exponents)])
    pairs = [[e + b - 2 * length for e, b in zip(seed['zero_energies'], beta)],
             [e + b - 2 * length for e, b in zip(seed['zero_energies'], exponents)]]
    return {'beta': beta, 'heights': heights, 'one_particle_orders': orders,
            'physical_pair_energies': pairs}
