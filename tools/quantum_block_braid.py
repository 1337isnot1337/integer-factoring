#!/usr/bin/env python3
"""Exact F16 block-braid polynomial crossing; exponential width construction."""

from __future__ import annotations

from bounded_spin_crossing import binary_path
from paired_path_crossing import MULTIPLY as M, path_entry, power, splitting_blocks


def braid_word(first, second):
    return tuple(i for j in range(second) for i in range(first + j - 1, j - 1, -1))


def polynomial_entry(four_paths, q=8):
    """Both strict endpoints; primal inputs P,Q and inverse-dual outputs R,S."""
    n = len(four_paths[0])
    assert n > 0 and n % 2 == 0 and all(len(p) == n and sum(p) == 0 for p in four_paths)
    splitting_blocks(q)
    frames = [
        {w: value for w in range(1 << n) if (value := path_entry(p, w, q, dual=i >= 2))}
        for i, p in enumerate(four_paths)
    ]
    states = {
        x | (y << n): [M[a][b]]
        for x, a in frames[0].items()
        for y, b in frames[1].items()
    }
    degree = 0
    for site in braid_word(n, n):
        degree += 1
        following = {w: [M[q][v] for v in values] + [0] for w, values in states.items()}
        for w, values in states.items():
            bits = (w >> site) & 3
            if bits not in (1, 2):
                continue
            diagonal = power(q, -1) if bits == 1 else q
            own = following[w]
            exchanged = following.setdefault(w ^ (3 << site), [0] * (degree + 1))
            for j, value in enumerate(values):
                own[j + 1] ^= M[diagonal][value]
                exchanged[j + 1] ^= value
        states = following
    mask = (1 << n) - 1
    result = [0] * (degree + 1)
    for word, values in states.items():
        scalar = M[frames[2].get(word & mask, 0)][frames[3].get(word >> n, 0)]
        if scalar:
            for j, value in enumerate(values):
                result[j] ^= M[scalar][value]
    normalization = power(q, -n * n // 2)
    result = [M[normalization][value] for value in result]
    # This affine scalar makes the zero endpoint strict and preserves the span.
    strict = [0] * (degree + 2)
    for j, value in enumerate(result):
        strict[j] ^= M[normalization][value]
        strict[j + 1] ^= M[1 ^ normalization][value]
    return strict


def binary_entry(width, labels, q=8):
    return polynomial_entry(tuple(binary_path(width, label) for label in labels), q)
