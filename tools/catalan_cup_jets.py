#!/usr/bin/env python3
"""Hasse jets at q=1 of the full Catalan crossing, over F2[u]/u^order."""

from __future__ import annotations

from math import comb

import numpy as np

from catalan_cup_crossing import (
    binary_product,
    crossing_layout,
    cup_basis,
    field_kronecker,
    triangular_dual,
)
from quantum_block_braid import braid_word


def cup_series(mates, words, order):
    words = tuple(words)
    words = np.array(words, dtype=np.uint64 if max(words) < 1 << 64 else object)
    result = np.zeros((order, len(words), len(mates)), dtype=np.uint8)
    binomials = np.array(
        [
            [comb(k, j) & 1 if j <= k else 0 for k in range(len(mates[0]) // 2 + 1)]
            for j in range(order)
        ],
        dtype=np.uint8,
    )
    for column, mate in enumerate(mates):
        valid = np.ones(len(words), dtype=bool)
        flips = np.zeros(len(words), dtype=int)
        for left, right in enumerate(mate):
            if left < right:
                a, b = (words >> left) & 1, (words >> right) & 1
                valid &= a != b
                flips += a == 0
        result[:, :, column] = binomials[:, flips] & valid[None, :]
    return result


def inverse_series(series):
    result = np.zeros_like(series)
    result[0] = triangular_dual(series[0])
    for degree in range(1, len(series)):
        rhs = np.zeros_like(series[0])
        for j in range(1, degree + 1):
            rhs ^= binary_product(series[j], result[degree - j])
        result[degree] = binary_product(result[0], rhs)
    return result


def decoder_series(length, order, selected=None):
    mates, words = cup_basis(length)
    global_mates, _ = cup_basis(2 * length)
    labels = tuple(range(len(mates))) if selected is None else selected
    inverse = inverse_series(cup_series(mates, words, order))[:, list(labels)]
    products = np.zeros((order, len(labels) ** 2, len(mates) ** 2), dtype=np.uint8)
    for degree in range(order):
        for j in range(degree + 1):
            products[degree] ^= field_kronecker(inverse[j], inverse[degree - j])
    physical_words = [p | (r << length) for p in words for r in words]
    frame = cup_series(global_mates, physical_words, order)
    result = np.zeros((order, len(labels) ** 2, len(global_mates)), dtype=np.uint8)
    for degree in range(order):
        for j in range(degree + 1):
            result[degree] ^= binary_product(products[j], frame[degree - j])
    return result


def crossing_jets(length, order, selected=None):
    """Axes u-degree,t-degree,P,Q,R,S; no coefficient ring lookup table."""
    assert order >= 1
    d, count, boundary, _, routes = crossing_layout(length, 1, selected)
    output = decoder_series(length, order, selected)
    states = np.zeros((order, count, d * d, 1), dtype=np.uint8)
    states[0, boundary, np.arange(d * d), 0] = 1
    for site in braid_word(length, length):
        permutation, targets, starts, nonloops = routes[site]
        following = np.zeros(
            states.shape[:-1] + (states.shape[-1] + 1,), dtype=np.uint8
        )
        following[..., :-1] = states
        following[1:, ..., :-1] ^= states[:-1]
        for degree in range(order):
            action = states[degree] & nonloops[:, None, None]
            following[degree, targets, :, 1:] ^= np.bitwise_xor.reduceat(
                action[permutation], starts, axis=0
            )
            if degree >= 2:
                paired = nonloops == 0
                following[degree, paired, :, 1:] ^= np.bitwise_xor.reduce(
                    states[: degree - 1, paired], axis=0
                )
        states = following
    decoded = np.zeros((order, states.shape[-1], d * d, d * d), dtype=np.uint8)
    for degree in range(order):
        for j in range(degree + 1):
            for t in range(states.shape[-1]):
                decoded[degree, t] ^= binary_product(
                    output[j], states[degree - j, ..., t]
                )
    h = length * length // 2
    gamma = [comb(h + j - 1, j) & 1 for j in range(order)]
    square = [gamma[j // 2] if j % 2 == 0 else 0 for j in range(order)]
    linear = [a ^ b for a, b in zip(gamma, square)]
    result = np.zeros((order, decoded.shape[1] + 1, d * d, d * d), dtype=np.uint8)
    for degree in range(order):
        for j in range(degree + 1):
            if square[j]:
                result[degree, :-1] ^= decoded[degree - j]
            if linear[j]:
                result[degree, 1:] ^= decoded[degree - j]
    while result.shape[1] > 1 and not np.any(result[:, -1]):
        result = result[:, :-1]
    return result.transpose(0, 1, 3, 2).reshape(order, result.shape[1], *([d] * 4))
