#!/usr/bin/env python3
"""Fixed-field Baxterized fusion and an explicit three-by-one scalar slide."""

from __future__ import annotations

import numpy as np

from bounded_spin_crossing import paths
from paired_path_crossing import MULTIPLY, power
from quantum_braid_path import local_action

FIELD = np.array(MULTIPLY, dtype=np.uint8)
Q = 8


def matrix_product(a, b):
    return np.bitwise_xor.reduce(FIELD[a[:, :, None], b[None, :, :]], axis=1)


def physical_generator(strands, site, q=Q):
    result = np.zeros((1 << strands, 1 << strands), dtype=np.uint8)
    for word in range(1 << strands):
        bits = word >> site & 3
        if bits in (1, 2):
            result[word, word] = power(q, -1) if bits == 1 else q
            result[word ^ (3 << site), word] = 1
    return result


def symmetrizer(strands, q=Q):
    """Ordinary JW projectors through four strands; all denominators are units."""
    assert 1 <= strands <= 4
    delta = q ^ power(q, -1)
    inverse_delta = power(delta, -1)

    def integer(h):
        return int(FIELD[power(q, h) ^ power(q, -h), inverse_delta])

    projector = np.eye(2, dtype=np.uint8)
    for h in range(1, strands):
        up = np.kron(np.eye(2, dtype=np.uint8), projector)
        generator = physical_generator(h + 1, h - 1, q)
        scalar = FIELD[integer(h), power(integer(h + 1), -1)]
        projector = (
            up ^ FIELD[scalar, matrix_product(up, matrix_product(generator, up))]
        )
    return projector


def apply_linear_gate(coefficients, constant, linear):
    result = np.zeros((len(coefficients) + 1,) + coefficients.shape[1:], dtype=np.uint8)
    for j, coefficient in enumerate(coefficients):
        result[j] ^= matrix_product(constant, coefficient)
        result[j + 1] ^= matrix_product(linear, coefficient)
    return result


def slide_scalars(q=Q):
    """p=q^2(1+z)(1+q^3*z)(1+q^4*z); r=q^2(1+q^3)(1+z)(1+q^3*z)."""
    common = [power(q, 2), power(q, 2) ^ 1, 1]
    p = common + [0]
    for j, c in enumerate(common):
        p[j + 1] ^= int(FIELD[power(q, 4), c])
    r = [int(FIELD[1 ^ power(q, 3), c]) for c in common] + [0]
    return p, r


def three_one_coefficients(q=Q):
    """W(z)F_L with application sites2,1,0 and spectral shifts0,-2,-4."""
    identity = np.eye(16, dtype=np.uint8)
    left = np.kron(np.eye(2, dtype=np.uint8), symmetrizer(3, q))
    coefficients, constant_braid = left[None, ...], left.copy()
    for site, exponent in zip((2, 1, 0), (0, -2, -4)):
        generator = physical_generator(4, site, q)
        braid = FIELD[q, identity] ^ generator
        inverse_braid = FIELD[power(q, -1), identity] ^ generator
        coefficients = apply_linear_gate(
            coefficients, inverse_braid, FIELD[power(q, exponent), braid]
        )
        constant_braid = matrix_product(braid, constant_braid)
    return coefficients, constant_braid, symmetrizer(4, q)


def normalized_three_one_coefficients(q=Q):
    """Polynomial after exact cancellation of q^2(1+z)(1+q^3*z)."""
    _, braid, f4 = three_one_coefficients(q)
    return np.array(
        [braid ^ FIELD[1 ^ power(q, 3), f4], FIELD[power(q, 4), braid]],
        dtype=np.uint8,
    )


def scalar_product(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] ^= int(FIELD[x, y])
    return result


def string_scalar(free_strands, q=Q):
    result = [1]
    p, _ = slide_scalars(q)
    for j in range(free_strands):
        result = scalar_product(
            result, [int(FIELD[c, power(q, j * k)]) for k, c in enumerate(p)]
        )
    return result


def normalized_string_scalar(free_strands, q=Q):
    result = [1]
    for j in range(free_strands):
        result = scalar_product(result, [1, power(q, j + 4)])
    return result


def fused_string_coefficients(free_strands, q=Q):
    """Vacuum quotient controls for a three-strand block crossing odd m strands.

    The external spectral parameters are q^j*z. Return the polynomial and
    the corresponding constant block braid, both already on the input f3 image.
    """
    assert free_strands >= 1 and free_strands % 2 == 1
    ps = paths(free_strands + 3)
    indices = {p: i for i, p in enumerate(ps)}
    delta = q ^ power(q, -1)
    embed = (0, 1, delta, delta ^ 1)
    generators = []
    for site in range(free_strands + 2):
        generator = np.zeros((len(ps), len(ps)), dtype=np.uint8)
        for j, p in enumerate(ps):
            for target, c in local_action(p, site):
                generator[indices[target], j] = embed[c]
        generators.append(generator)
    identity = np.eye(len(ps), dtype=np.uint8)
    e0, e1 = generators[:2]
    projector = (
        identity
        ^ e0
        ^ e1
        ^ FIELD[power(delta, -1), matrix_product(e0, e1) ^ matrix_product(e1, e0)]
    )
    coefficients, constant_braid = projector[None, ...], projector.copy()
    for j in range(free_strands):
        for site, exponent in zip((j + 2, j + 1, j), (j, j - 2, j - 4)):
            generator = generators[site]
            braid = FIELD[q, identity] ^ generator
            inverse = FIELD[power(q, -1), identity] ^ generator
            coefficients = apply_linear_gate(
                coefficients, inverse, FIELD[power(q, exponent), braid]
            )
            constant_braid = matrix_product(braid, constant_braid)
    return coefficients, constant_braid
