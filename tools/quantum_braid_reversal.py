#!/usr/bin/env python3
"""Reversal-stable binary alphabets, leading factors, and nilpotent braid dressings."""

from __future__ import annotations

from functools import lru_cache

from bounded_spin_crossing import suffix_counts, unrank_path
from quantum_braid_path import MULTIPLY, crossing_weights


def reverse_path(path):
    return tuple(-step for step in path[::-1])


@lru_cache(None)
def _pair_data(length):
    half = length // 2
    multiplicities = tuple(suffix_counts(half, h)[half][0] for h in range(4))
    return multiplicities, sum(m * (m - 1) // 2 for m in multiplicities)


def _unordered_pair(size, index):
    assert 0 <= index < size * (size - 1) // 2
    low, high = 0, size - 1
    while low + 1 < high:
        middle = (low + high) // 2
        if middle * (2 * size - middle - 1) // 2 <= index:
            low = middle
        else:
            high = middle
    return low, low + 1 + index - low * (2 * size - low - 1) // 2


def reversal_pair(length, index):
    multiplicities, total = _pair_data(length)
    assert 0 <= index < total
    for h, size in enumerate(multiplicities):
        count = size * (size - 1) // 2
        if index < count:
            i, j = _unordered_pair(size, index)
            first, second = unrank_path(length // 2, i, h), unrank_path(
                length // 2, j, h
            )
            path = first + reverse_path(second)
            return path, reverse_path(path)
        index -= count
    raise AssertionError("Pair index was not consumed.")


def fixed_path(length, index):
    half = length // 2
    counts = [[1] * 4]
    for _ in range(half):
        counts.append(
            [
                sum(counts[-1][j] for j in (h - 1, h + 1) if 0 <= j <= 3)
                for h in range(4)
            ]
        )
    assert 0 <= index < counts[half][0]
    prefix, height = [], 0
    for remaining in range(half - 1, -1, -1):
        up = counts[remaining][height + 1] if height < 3 else 0
        if index < up:
            prefix.append(1)
            height += 1
        else:
            index -= up
            prefix.append(-1)
            height -= 1
    prefix = tuple(prefix)
    return prefix + reverse_path(prefix)


@lru_cache(None)
def binary_layout(width):
    assert width >= 1
    length = 4
    while suffix_counts(length)[length][0] < 1 << width:
        length += 2
    _, pairs = _pair_data(length)
    return length, min(pairs, 1 << (width - 1))


def binary_path(width, label):
    length, pairs = binary_layout(width)
    assert 0 <= label < 1 << width
    if label < 2 * pairs:
        return reversal_pair(length, label // 2)[label & 1]
    return fixed_path(length, label - 2 * pairs)


def reverse_label(width, label):
    _, pairs = binary_layout(width)
    assert 0 <= label < 1 << width
    return label ^ 1 if label < 2 * pairs else label


def diagram_norm(path):
    height, value = 0, 1
    for step in path:
        if step == -1:
            value = MULTIPLY[value][(0, 2, 3, 2)[height]]
        height += step
    return value


def leading_entry(four_paths):
    p, q, r, s = four_paths
    if q != reverse_path(p) or s != reverse_path(r):
        return 0
    return MULTIPLY[diagram_norm(p)][(0, 1, 3, 2)[diagram_norm(r)]]


def crossing_leading_entry(four_paths):
    n = len(four_paths[0])
    return MULTIPLY[crossing_weights(n)[-1]][leading_entry(four_paths)]


def binary_polynomial_entry(width, labels):
    from quantum_braid_path import polynomial_entry

    return polynomial_entry(tuple(binary_path(width, label) for label in labels))


def nilpotent_dressing(tensor, scalar):
    """Conjugate the normalized row operator by I+scalar*z*(1+z)*P."""
    import numpy as np

    def multiply(a, b):
        # uint8 dot products preserve parity exactly, including wraparound.
        x, y, u, v = a & 1, a >> 1, b & 1, b >> 1
        first, second, mixed = (x @ u) & 1, (y @ v) & 1, ((x ^ y) @ (u ^ v)) & 1
        return (first ^ second) + 2 * (mixed ^ first)

    data = np.asarray(tensor, dtype=np.uint8)
    d = data.shape[1]
    row = data.transpose(0, 1, 2, 4, 3).reshape(len(data), d * d, d * d)
    p = row[-1]
    if np.any(multiply(p, p)):
        raise ValueError("The leading row operator must be square-zero.")
    field = np.array(MULTIPLY, dtype=np.uint8)
    commutators = np.array([multiply(p, a) ^ multiply(a, p) for a in row])
    sandwiches = np.array([multiply(multiply(p, a), p) for a in row])
    result = np.zeros((len(row) + 4, d * d, d * d), dtype=np.uint8)
    result[: len(row)] = row
    result[1 : len(row) + 1] ^= field[scalar, commutators]
    result[2 : len(row) + 2] ^= field[scalar, commutators]
    square = MULTIPLY[scalar][scalar]
    result[2 : len(row) + 2] ^= field[square, sandwiches]
    result[4 : len(row) + 4] ^= field[square, sandwiches]
    while len(result) > 1 and not np.any(result[-1]):
        result = result[:-1]
    return result.reshape(len(result), d, d, d, d).transpose(0, 1, 2, 4, 3)
