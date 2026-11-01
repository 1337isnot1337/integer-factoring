#!/usr/bin/env python3
"""Exact leading central auxiliary ranks, proved by a finite fusion automaton."""

from __future__ import annotations

from nilpotent_seam_rank_bounds import vacuum_dimension


def four_block_rank(blocks):
    """[z^m](4-24z+55z^2-52z^3+22z^4)/((1-z)(1-3z+z^2)^2)."""
    assert blocks >= 0
    values = [4, 4, 15, 53, 178]
    for _ in range(5, blocks + 1):
        values.append(
            7 * values[-1]
            - 17 * values[-2]
            + 17 * values[-3]
            - 7 * values[-4]
            + values[-5]
        )
    return values[blocks]


def top_auxiliary_rank(length):
    """Rank of [x^(2n-floor(n/4))] on auxiliary tensor full vacuum space.

    This precedes the boundary contractions and is an upper bound for their
    row factors; it does not assert equality for a full crossing coefficient.
    """
    assert length > 0 and length % 2 == 0
    blocks, tail = divmod(length, 4)
    if tail == 0:
        return four_block_rank(blocks)
    return 2 * vacuum_dimension(blocks)
