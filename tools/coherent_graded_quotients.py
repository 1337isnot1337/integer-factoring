#!/usr/bin/env python3
"""Coherent physical codes dual to homogeneous Boolean algebra quotients.

The first k physical generators map identically to k virtual generators.
Every later physical generator can map to an arbitrary odd homogeneous
virtual polynomial of the same height. Dualizing this surjection gives
an embedding and a coordinate retraction, with deletion closure built in.
"""
from __future__ import annotations

from coherent_hot_compounds import height, support
from coherent_hot_mixed import validate_maps


def multiply(first, second):
    answer = 0
    for a in support(first):
        for b in support(second):
            if a & b == 0:
                answer ^= 1 << (a | b)
    return answer


def monomial_parameters(virtual_modes, physical_modes):
    if not 1 <= virtual_modes <= physical_modes <= 16:
        raise ValueError('Use one through sixteen modes with virtual modes first.')
    return [(a, word) for a in range(virtual_modes, physical_modes)
            for word in range(1, 1 << virtual_modes) if height(word) == 2 * a + 1]


def parameter_images(virtual_modes, physical_modes, mask):
    parameters = monomial_parameters(virtual_modes, physical_modes)
    if not isinstance(mask, int) or not 0 <= mask < 1 << len(parameters):
        raise ValueError('The parameter mask does not fit the homogeneous monomial list.')
    images = [1 << (1 << a) if a < virtual_modes else 0 for a in range(physical_modes)]
    for index, (a, word) in enumerate(parameters):
        if mask >> index & 1:
            images[a] ^= 1 << word
    return images


def quotient_maps(images, virtual_modes):
    images = tuple(images)
    modes = len(images)
    if not 1 <= virtual_modes <= modes <= 16:
        raise ValueError('Use one through sixteen modes with virtual modes first.')
    q = 1 << virtual_modes
    if images[:virtual_modes] != tuple(1 << (1 << a) for a in range(virtual_modes)):
        raise ValueError('The first virtual generators must map identically.')
    for a, polynomial in enumerate(images):
        if not isinstance(polynomial, int) or polynomial < 0 or polynomial.bit_length() > q:
            raise ValueError('Generator images must be virtual Boolean polynomials.')
        if any(height(word) != 2 * a + 1 for word in support(polynomial)):
            raise ValueError('Each generator image must preserve height.')
        if multiply(polynomial, polynomial) != 0:
            raise ArithmeticError('A generator image did not square to zero.')
    monomial_images = [0] * (1 << modes)
    monomial_images[0] = 1
    embedding = [0] * q
    for word in range(1 << modes):
        if word:
            bit = word & -word
            a = bit.bit_length() - 1
            monomial_images[word] = multiply(monomial_images[word ^ bit], images[a])
        for virtual_word in support(monomial_images[word]):
            embedding[virtual_word] ^= 1 << word
    retraction = [1 << word for word in range(q)]
    _, _, heights, checks = validate_maps(embedding, retraction, modes)
    assert list(heights) == [height(word) for word in range(q)]
    return {'modes': modes, 'virtual_modes': virtual_modes, 'q': q,
            'generator_images': list(images), 'embedding': embedding, 'retraction': retraction,
            'heights': list(heights), 'deletion_checks': checks,
            'nonzero_monomial_images': sum(value != 0 for value in monomial_images)}
