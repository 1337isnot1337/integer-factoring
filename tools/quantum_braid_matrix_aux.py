#!/usr/bin/env python3
"""A fixed two-state noncommuting auxiliary deformation of the F4 braid."""

from __future__ import annotations

import numpy as np

from bounded_spin_crossing import paths
from quantum_block_braid import braid_word
from quantum_braid_path import MULTIPLY, local_action

FIELD = np.array(MULTIPLY, dtype=np.uint8)
Q = np.array([[0, 1], [1, 2]], dtype=np.uint8)
N = np.array([[0, 1], [0, 0]], dtype=np.uint8)
IDENTITY = np.eye(2, dtype=np.uint8)


def normalization_row(length):
    phase = -length * length // 2 % 5
    normalization = IDENTITY.copy()
    for _ in range(phase):
        normalization = matrix_product(normalization, Q)
    tau = (0, 1, 0, 2, 3)[phase]
    return np.array(
        [normalization[0, j] ^ FIELD[tau, normalization[1, j]] for j in range(2)],
        dtype=np.uint8,
    )


def leading_entry(four_paths):
    p, q, r, s = four_paths
    zigzag = (1, -1) * (len(p) // 2)
    height, maximum = 0, 0
    for step in s:
        height += step
        maximum = max(maximum, height)
    return int(p == zigzag and q == r and maximum <= 2)


def matrix_product(a, b):
    return np.array(
        [
            [FIELD[a[i, 0], b[0, j]] ^ FIELD[a[i, 1], b[1, j]] for j in range(2)]
            for i in range(2)
        ],
        dtype=np.uint8,
    )


def auxiliary_action(a, values):
    result = np.zeros_like(values)
    for i in range(2):
        for j in range(2):
            result[i] ^= FIELD[a[i, j], values[j]]
    return result


def coefficient_tensor(length):
    """Local gate (I+t(Q-I))*I + (t*I+t*(1-t)*N)*e_i."""
    assert length > 0 and length % 2 == 0
    ps, qs = paths(length), paths(2 * length)
    indices = {p: i for i, p in enumerate(qs)}
    d = len(ps)
    boundary = np.array([indices[p + q] for p in ps for q in ps])
    states = np.zeros((2, len(qs), d * d, 1), dtype=np.uint8)
    states[0, boundary, np.arange(d * d), 0] = 1
    for site in braid_word(length, length):
        diagonal = np.zeros(len(qs), dtype=np.uint8)
        source, target = [], []
        for i, p in enumerate(qs):
            for q, c in local_action(p, site):
                if p == q:
                    diagonal[i] = c
                else:
                    source.append(i)
                    target.append(indices[q])
        source, target = np.array(source, dtype=np.int64), np.array(
            target, dtype=np.int64
        )
        action = FIELD[diagonal[None, :, None, None], states]
        action[:, target] ^= FIELD[
            diagonal[source][None, :, None, None], states[:, source]
        ]
        following = np.zeros(
            states.shape[:-1] + (states.shape[-1] + 2,), dtype=np.uint8
        )
        following[..., :-2] = states
        following[..., 1:-1] ^= auxiliary_action(
            Q ^ IDENTITY, states
        ) ^ auxiliary_action(IDENTITY ^ N, action)
        following[..., 2:] ^= auxiliary_action(N, action)
        states = following
    left = normalization_row(length)
    selected = FIELD[left[0], states[0, boundary]] ^ FIELD[left[1], states[1, boundary]]
    result = selected.transpose(2, 1, 0).reshape(2 * length * length + 1, d, d, d, d)
    while len(result) > 1 and not np.any(result[-1]):
        result = result[:-1]
    assert len(result) - 1 <= 3 * length * length // 2
    return result
