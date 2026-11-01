#!/usr/bin/env python3
"""Exact all-degree rank envelopes for the four-state braid seam."""

from __future__ import annotations

from math import comb


def vacuum_dimension(pairs):
    """Number of height-three vacuum paths on twice this many strands."""
    current, following = 1, 1
    for _ in range(pairs):
        current, following = following, 3 * following - current
    return current


def central_weights(length):
    """Vacuum-weighted ternary degree patterns; eight integer counters per degree.

    Degree0 leaves the mandatory-N run unchanged; degree1 resets it; degree2
    increments it and forces a cup. Four mandatory N's without a reset vanish.
    Each state holds sums of (v_(2f),v_(2f+2)) for its free-pair counts f.
    """
    states = {(0, 0): (1, 1)}
    for _ in range(length):
        following = {}

        def add(key, pair):
            old = following.get(key, (0, 0))
            following[key] = old[0] + pair[0], old[1] + pair[1]

        for (run, degree), (x, y) in states.items():
            free = y, 3 * y - x
            add((run, degree), free)
            add((0, degree + 1), free)
            if run < 3:
                add((run + 1, degree + 2), (x, y))
        states = following
    result = [0] * (2 * length + 1)
    for (_, degree), (value, _) in states.items():
        result[degree] += value
    return result


def rank_envelope(length, central_degree, boundary_degree, weights=None):
    """Bounds before the trivial L^2 cap; weights can be reused across a grid."""
    if weights is None:
        weights = central_weights(length)
    return (
        4
        * min(boundary_degree + 1, length - boundary_degree + 1)
        * weights[central_degree],
        16
        * (central_degree + 1)
        * comb(length, boundary_degree)
        * vacuum_dimension(boundary_degree),
    )
