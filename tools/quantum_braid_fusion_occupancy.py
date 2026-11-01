#!/usr/bin/env python3
"""A complete crossing resolution by orthogonal central fusion occupancies."""

from __future__ import annotations

from math import comb

import numpy as np

from paired_path_crossing import power
from quantum_braid_fusion_seam import FIELD, INVERSE_DELTA, Q
from quantum_braid_graded_seam import _apply_e, _layout, fusion_rank


def central_rank(length, triplets):
    assert 0 <= triplets <= length
    return comb(length, triplets) * fusion_rank(triplets)


def central_projectors(length):
    """Q_k=[u^k] product(p0+u*p2), on the full 2n-strand vacuum quotient."""
    assert length > 0 and length % 2 == 0
    _, states_count, _, routes = _layout(length)
    projectors = np.eye(states_count, dtype=np.uint8)[..., None]
    for site in range(0, 2 * length, 2):
        cup = FIELD[INVERSE_DELTA, _apply_e(projectors, routes[site])]
        following = np.zeros(
            projectors.shape[:-1] + (projectors.shape[-1] + 1,), dtype=np.uint8
        )
        following[..., :-1] = cup
        following[..., 1:] ^= projectors ^ cup
        projectors = following
    return projectors.transpose(2, 0, 1)


def coefficient_tensor(length):
    """Axes k,b,P,Q,R,S; complete normalized two-parameter central resolution.

    At u=q^2,y=1 the polynomial is the crossing E. Central cells are
    p0+u*p2; boundary cells are qI+y*e; all other cells are qI+e.
    """
    assert length > 0 and length % 2 == 0
    d, states_count, boundary, routes = _layout(length)
    states = np.zeros((states_count, d * d, 1, 1), dtype=np.uint8)
    states[boundary, np.arange(d * d), 0, 0] = 1
    for diagonal in range(2 * length - 1):
        for row in range(length):
            column = diagonal - row
            if not 0 <= column < length:
                continue
            action = _apply_e(states, routes[length + row - 1 - column])
            if diagonal == length - 1:
                cup = FIELD[INVERSE_DELTA, action]
                following = np.zeros(
                    states.shape[:-2] + (states.shape[-2] + 1, states.shape[-1]),
                    dtype=np.uint8,
                )
                following[..., :-1, :] = cup
                following[..., 1:, :] ^= states ^ cup
            elif row == column:
                following = np.zeros(
                    states.shape[:-1] + (states.shape[-1] + 1,), dtype=np.uint8
                )
                following[..., :-1] = FIELD[Q, states]
                following[..., 1:] ^= action
            else:
                following = FIELD[Q, states] ^ action
            states = following
    tensor = (
        states[boundary]
        .transpose(2, 3, 1, 0)
        .reshape(length + 1, length + 1, d, d, d, d)
    )
    return FIELD[power(Q, -length - length * length // 2), tensor]


def rank_bounds(length, triplets, boundary_degree):
    """Both actual cut bounds, including occupancy convolution channels."""
    h = min(boundary_degree + 1, length - boundary_degree + 1)
    lo = max(0, triplets - length // 2)
    hi = min(length // 2, triplets)
    # Individual block tangles need not be vacuum projectors at the interface.
    # Keep every feasible left/right occupancy split, including a single triplet.
    channels = hi - lo + 1
    return (
        h * central_rank(length, triplets),
        channels * comb(length, boundary_degree) * fusion_rank(2 * boundary_degree),
    )
