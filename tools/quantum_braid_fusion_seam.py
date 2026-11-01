#!/usr/bin/env python3
"""Triplet-projector middle channels for the height-three braid quotient."""

from __future__ import annotations

import numpy as np

from bounded_spin_crossing import paths
from paired_path_crossing import MULTIPLY, power
from quantum_braid_path import EMBED, local_action

FIELD = np.array(MULTIPLY, dtype=np.uint8)
Q = 8
DELTA = 7
INVERSE_DELTA = power(DELTA, 14)


def coefficient_tensor(length):
    """Axes b,P,Q,R,S; central f2, boundary qI+y*e, other gates qI+e.

    The scalar q^(n-n^2/2) agrees with the qI+x*e deformation at
    the triplet value x=q/delta. These are candidate constant tensors;
    their whole polynomial is not assigned strict crossing endpoints.
    """
    assert length > 0 and length % 2 == 0
    ps, qs = paths(length), paths(2 * length)
    indices = {p: i for i, p in enumerate(qs)}
    d = len(ps)
    boundary = np.array([indices[p + q] for p in ps for q in ps])
    states = np.zeros((len(qs), d * d, 1), dtype=np.uint8)
    states[boundary, np.arange(d * d), 0] = 1
    for _, row, column in sorted(
        (j + k, j, k) for j in range(length) for k in range(length)
    ):
        site = length + row - 1 - column
        central, cross_cut = row + column == length - 1, row == column
        diagonal = np.zeros(len(qs), dtype=np.uint8)
        source, target = [], []
        for i, p in enumerate(qs):
            for q, scalar in local_action(p, site):
                if q == p:
                    diagonal[i] = EMBED[scalar]
                else:
                    source.append(i)
                    target.append(indices[q])
        source, target = np.array(source, dtype=int), np.array(target, dtype=int)
        action = FIELD[diagonal[:, None, None], states]
        action[target] ^= FIELD[diagonal[source, None, None], states[source]]
        if central:
            states = states ^ FIELD[INVERSE_DELTA, action]
        elif cross_cut:
            following = np.zeros((len(qs), d * d, states.shape[-1] + 1), dtype=np.uint8)
            following[..., :-1] = FIELD[Q, states]
            following[..., 1:] ^= action
            states = following
        else:
            states = FIELD[Q, states] ^ action
    tensor = states[boundary].transpose(2, 1, 0).reshape(length + 1, d, d, d, d)
    return FIELD[power(Q, length - length * length // 2), tensor]
