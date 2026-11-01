#!/usr/bin/env python3
"""Graded braid channels from invertible three-strand fusion projectors."""

from __future__ import annotations

from functools import lru_cache
from math import comb

import numpy as np

from bounded_spin_crossing import paths
from paired_path_crossing import power
from quantum_braid_fusion_seam import FIELD, INVERSE_DELTA, Q
from quantum_braid_path import EMBED, local_action


def fusion_rank(free_strands):
    """Vacuum multiplicity after the invertible charge is canceled: F_(k-1)."""
    if free_strands == 0:
        return 1
    a, b = 0, 1
    for _ in range(free_strands - 1):
        a, b = b, a + b
    return a


def free_counts(length):
    return tuple(k for k in range(length + 1) if (length - k) % 3 == 0)


def projector_sites(length, free_left, free_right):
    assert free_left in free_counts(length) and free_right in free_counts(length)
    return tuple(range(free_left, length, 3)) + tuple(
        range(length + free_right, 2 * length, 3)
    )


@lru_cache(maxsize=1)
def _layout(length):
    ps, qs = paths(length), paths(2 * length)
    indices = {p: i for i, p in enumerate(qs)}
    boundary = np.array([indices[p + q] for p in ps for q in ps])
    routes = []
    for site in range(2 * length - 1):
        diagonal = np.zeros(len(qs), dtype=np.uint8)
        source, target = [], []
        for i, p in enumerate(qs):
            for q, c in local_action(p, site):
                if q == p:
                    diagonal[i] = EMBED[c]
                else:
                    source.append(i)
                    target.append(indices[q])
        routes.append(
            (diagonal, np.array(source, dtype=int), np.array(target, dtype=int))
        )
    return len(ps), len(qs), boundary, routes


def _apply_e(values, route):
    diagonal, source, target = route
    shape = (-1,) + (1,) * (values.ndim - 1)
    result = FIELD[diagonal.reshape(shape), values]
    result[target] ^= FIELD[diagonal[source].reshape(shape), values[source]]
    return result


def _apply_f3(values, start, routes):
    a, b = _apply_e(values, routes[start]), _apply_e(values, routes[start + 1])
    return (
        values
        ^ a
        ^ b
        ^ FIELD[
            INVERSE_DELTA,
            _apply_e(b, routes[start]) ^ _apply_e(a, routes[start + 1]),
        ]
    )


def coefficient_tensor(length, free_left, free_right):
    """Axes b,P,Q,R,S; free prefixes followed by f3 triples in each block.

    The old central layer is replaced as a whole. Boundary gates are qI+y*e;
    other gates are qI+e. These are candidate pieces, not a strict endpoint pair.
    """
    assert length > 0 and length % 2 == 0
    sites = projector_sites(length, free_left, free_right)
    d, states_count, boundary, routes = _layout(length)
    states = np.zeros((states_count, d * d, 1), dtype=np.uint8)
    states[boundary, np.arange(d * d), 0] = 1
    for diagonal in range(2 * length - 1):
        if diagonal == length - 1:
            for site in sites:
                states = _apply_f3(states, site, routes)
            continue
        for row in range(length):
            column = diagonal - row
            if not 0 <= column < length:
                continue
            action = _apply_e(states, routes[length + row - 1 - column])
            if row == column:
                following = np.zeros(
                    states.shape[:-1] + (states.shape[-1] + 1,), dtype=np.uint8
                )
                following[..., :-1] = FIELD[Q, states]
                following[..., 1:] ^= action
                states = following
            else:
                states = FIELD[Q, states] ^ action
    tensor = states[boundary].transpose(2, 1, 0).reshape(length + 1, d, d, d, d)
    return FIELD[power(Q, length - length * length // 2), tensor]


def rank_bounds(length, free_left, free_right, boundary_degree):
    return (
        min(boundary_degree + 1, length - boundary_degree + 1)
        * fusion_rank(free_left + free_right),
        comb(length, boundary_degree) * fusion_rank(2 * boundary_degree),
    )
