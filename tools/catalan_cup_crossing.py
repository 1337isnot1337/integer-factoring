#!/usr/bin/env python3
"""Full cup-diagram primal frames, coordinate duals, and ordered crossings."""

from __future__ import annotations

from functools import lru_cache
from math import comb

import numpy as np

from paired_path_crossing import power
from quantum_block_braid import braid_word
from quantum_braid_spectral_fusion import FIELD, matrix_product


def label_count(length):
    assert length >= 0 and length % 2 == 0
    half = length // 2
    return comb(length, half) // (half + 1)


def binary_length(width):
    assert width >= 0
    length = 2
    while label_count(length) < 1 << width:
        length += 2
    return length


@lru_cache(maxsize=8)
def cup_basis(length):
    """Noncrossing mates and opener words, in increasing Dyck-area order."""
    assert length > 0 and length % 2 == 0
    result = []

    def visit(position, opens, stack, mates, word, area):
        if position == length:
            result.append((area, word, tuple(mates)))
            return
        if opens < length // 2:
            visit(
                position + 1,
                opens + 1,
                stack + (position,),
                mates,
                word | (1 << position),
                area + len(stack) + 1,
            )
        if stack:
            following = mates.copy()
            left = stack[-1]
            following[left], following[position] = position, left
            visit(
                position + 1, opens, stack[:-1], following, word, area + len(stack) - 1
            )

    visit(0, 0, (), [-1] * length, 0, 0)
    result.sort()
    return tuple(x[2] for x in result), tuple(x[1] for x in result)


def cup_frame(mates, words, q=1):
    """Rows are physical words, columns are cup tensors; no Gram inversion."""
    assert q
    words = np.array(words, dtype=object)
    frame = np.zeros((len(words), len(mates)), dtype=np.uint8)
    powers = np.array(
        [power(q, k) for k in range(len(mates[0]) // 2 + 1)], dtype=np.uint8
    )
    for column, mate in enumerate(mates):
        valid = np.ones(len(words), dtype=bool)
        exponent = np.zeros(len(words), dtype=int)
        for left, right in enumerate(mate):
            if left < right:
                a, b = (words >> left) & 1, (words >> right) & 1
                valid &= a != b
                exponent += a == 0
        frame[:, column] = np.where(valid, powers[exponent], 0)
    return frame


def triangular_dual(minor):
    dimension = len(minor)
    assert np.all(np.diag(minor) == 1) and not np.any(np.tril(minor, -1))
    inverse = np.eye(dimension, dtype=np.uint8)
    for i in range(dimension - 1, -1, -1):
        if i + 1 < dimension:
            inverse[i] ^= np.bitwise_xor.reduce(
                FIELD[minor[i, i + 1 :, None], inverse[i + 1 :]], axis=0
            )
    return inverse


def field_kronecker(a, b):
    return FIELD[a[:, None, :, None], b[None, :, None, :]].reshape(
        a.shape[0] * b.shape[0], a.shape[1] * b.shape[1]
    )


def multiply(a, b, q):
    # Modulo-256 overflow preserves parity for binary matrices.
    if q == 1:
        return (a @ b) & 1
    return matrix_product(a, b)


def cup_action(mate, site, q=1):
    if mate[site] == site + 1:
        return mate, q ^ power(q, -1)
    a, b = mate[site], mate[site + 1]
    target = list(mate)
    target[site], target[site + 1] = site + 1, site
    target[a], target[b] = b, a
    return tuple(target), 1


@lru_cache(maxsize=1)
def crossing_layout(length, q=1):
    mates, words = cup_basis(length)
    global_mates, _ = cup_basis(2 * length)
    dimension = len(mates)
    inverse = triangular_dual(cup_frame(mates, words, q))
    index = {p: i for i, p in enumerate(global_mates)}
    boundary = np.array(
        [index[p + tuple(i + length for i in r)] for p in mates for r in mates]
    )
    selected = [p | (r << length) for p in words for r in words]
    output = multiply(
        field_kronecker(inverse, inverse), cup_frame(global_mates, selected, q), q
    )
    routes = []
    for site in range(2 * length - 1):
        targets, scalars = zip(*(cup_action(mate, site, q) for mate in global_mates))
        targets = np.array([index[p] for p in targets])
        order = np.argsort(targets, kind="stable")
        unique, starts = np.unique(targets[order], return_index=True)
        routes.append((order, unique, starts, np.array(scalars, dtype=np.uint8)))
    return dimension, len(global_mates), boundary, output, routes


def crossing_tensor(length, q=1):
    """Axes degree,P,Q,R,S; strict identity/crossing endpoints, over F2 at q=1."""
    dimension, states_count, boundary, output, routes = crossing_layout(length, q)
    states = np.zeros((states_count, dimension * dimension, 1), dtype=np.uint8)
    states[boundary, np.arange(dimension * dimension), 0] = 1
    for site in braid_word(length, length):
        order, targets, starts, scalars = routes[site]
        following = np.zeros(
            states.shape[:-1] + (states.shape[-1] + 1,), dtype=np.uint8
        )
        following[..., :-1] = FIELD[q, states]
        action = FIELD[scalars[:, None, None], states]
        following[targets, :, 1:] ^= np.bitwise_xor.reduceat(
            action[order], starts, axis=0
        )
        states = following
    tensor = np.array(
        [multiply(output, states[..., j], q) for j in range(states.shape[-1])],
        dtype=np.uint8,
    )
    gamma = power(q, -length * length // 2)
    tensor = FIELD[gamma, tensor]
    strict = np.zeros((len(tensor) + 1,) + tensor.shape[1:], dtype=np.uint8)
    strict[:-1] = FIELD[gamma, tensor]
    strict[1:] ^= FIELD[1 ^ gamma, tensor]
    while len(strict) > 1 and not np.any(strict[-1]):
        strict = strict[:-1]
    return strict.transpose(0, 2, 1).reshape(len(strict), *([dimension] * 4))
