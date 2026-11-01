#!/usr/bin/env python3
"""A four-state nilpotent auxiliary on the central seam of a block braid."""

from __future__ import annotations

import numpy as np

from bounded_spin_crossing import paths
from quantum_braid_path import MULTIPLY, local_action

FIELD = np.array(MULTIPLY, dtype=np.uint8)
IDENTITY = np.eye(4, dtype=np.uint8)
Q = np.array(((0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 2, 0), (1, 0, 0, 2)), dtype=np.uint8)
N = np.eye(4, k=1, dtype=np.uint8)


def matrix_product(a, b):
    return np.bitwise_xor.reduce(FIELD[a[:, :, None], b[None, :, :]], axis=1)


def auxiliary_action(a, values):
    result = np.zeros_like(values)
    for i, j in zip(*np.nonzero(a)):
        result[i] ^= FIELD[a[i, j], values[j]]
    return result


def normalization_row(length):
    phase = -length * length // 2 % 5
    normal = IDENTITY.copy()
    for _ in range(phase):
        normal = matrix_product(normal, Q)
    return normal[0] ^ FIELD[(0, 1, 0, 2, 3)[phase], normal[3]]


def coefficient_tensor(length):
    """Axes a,b,P,Q,R,S; retain the theoretical top a even if it is zero.

    Central gates: (I+x(Q-I))*I+(x*I+x*(1+x)*N)*e_i.
    Cross-cut gates: Q*I+y*e_i. All other gates: Q*I+e_i.
    Cells are ordered by anti-diagonal, then row. The two seams are disjoint.
    """
    assert length > 0 and length % 2 == 0
    ps, qs = paths(length), paths(2 * length)
    indices = {p: i for i, p in enumerate(qs)}
    d = len(ps)
    boundary = np.array([indices[p + q] for p in ps for q in ps])
    states = np.zeros((4, len(qs), d * d, 1, 1), dtype=np.uint8)
    states[0, boundary, np.arange(d * d), 0, 0] = 1
    for _, row, column in sorted(
        (j + k, j, k) for j in range(length) for k in range(length)
    ):
        site = length + row - 1 - column
        central, cross_cut = row + column == length - 1, row == column
        assert not (central and cross_cut)
        diagonal = np.zeros(len(qs), dtype=np.uint8)
        source, target = [], []
        for i, p in enumerate(qs):
            for q, scalar in local_action(p, site):
                if q == p:
                    diagonal[i] = scalar
                else:
                    source.append(i)
                    target.append(indices[q])
        source, target = np.array(source, dtype=int), np.array(target, dtype=int)
        action = FIELD[diagonal[None, :, None, None, None], states]
        action[:, target] ^= FIELD[
            diagonal[source][None, :, None, None, None], states[:, source]
        ]
        if central:
            shape = list(states.shape)
            shape[3] += 2
            following = np.zeros(shape, dtype=np.uint8)
            following[..., :-2, :] = states
            following[..., 1:-1, :] ^= auxiliary_action(
                Q ^ IDENTITY, states
            ) ^ auxiliary_action(IDENTITY ^ N, action)
            following[..., 2:, :] ^= auxiliary_action(N, action)
        elif cross_cut:
            shape = list(states.shape)
            shape[4] += 1
            following = np.zeros(shape, dtype=np.uint8)
            following[..., :-1] = auxiliary_action(Q, states)
            following[..., 1:] ^= action
        else:
            following = auxiliary_action(Q, states) ^ action
        states = following
    left = normalization_row(length)
    selected = np.zeros(states.shape[1:], dtype=np.uint8)
    for i in range(4):
        selected ^= FIELD[left[i], states[i]]
    tensor = (
        selected[boundary]
        .transpose(2, 3, 1, 0)
        .reshape(2 * length + 1, length + 1, d, d, d, d)
    )
    degree_bound = 2 * length - length // 4
    assert not np.any(tensor[degree_bound + 1 :])
    return tensor[: degree_bound + 1]
