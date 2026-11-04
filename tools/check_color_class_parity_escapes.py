#!/usr/bin/env python3
"""Finite controls for parity-breaking variants of the color-class cube.

The one-hot color-class tensor preserves individual color parity, so the
sector obstruction applies at large q.  This checker tests the natural
parity-breaking alternatives at q=4: support forgetfulness, grouped colors,
an endpoint-restoring antipodal delta, pair-weighted support, and a
colour-averaged mixture of one-hot tensors.

The point is the endpoint and the finite Smith spectrum, not a new pressure
theorem.  A variant that fails the individual endpoint is not a candidate.
An endpoint-correct variant with many zero orders or more low orders than
the one-hot control is not an improvement.
"""
from __future__ import annotations

from collections import Counter
from itertools import product
from math import log2
import json
from pathlib import Path

import numpy as np

from check_color_class_cube_profile import (
    GAUSSIAN, evaluate, LEFT, RIGHT, MULTIPLY, INVERSE,
    color_matrix, rank_over_f256, block_toeplitz, specialized_profile)

Q = 4
L = 7
Y, Z = 2, 3
N = Q ** 4
COLUMNS = np.array(list(product(range(Q), repeat=4)), dtype=np.uint8)
EVALUATED = np.array([evaluate(GAUSSIAN[m], Y, Z, L)
                      for m in range(256)], dtype=np.uint8)
SHIFT_INDEX = []
for c in range(Q):
    shifted = COLUMNS ^ c
    keys = {tuple(row): i for i, row in enumerate(COLUMNS)}
    SHIFT_INDEX.append(np.array([keys[tuple(row)] for row in shifted]))


def block_profile(blocks):
    ranks = [0]
    for cut in range(1, L + 1):
        ranks.append(rank_over_f256(block_toeplitz(blocks, cut)))
    cumulative = {order: ranks[order + 1] - ranks[order]
                  for order in range(L)}
    profile = {}
    previous = 0
    for order in range(L):
        profile[order] = cumulative[order] - previous
        previous = cumulative[order]
    return {order: count for order, count in profile.items() if count}


def group_blocks(groups):
    result = np.zeros((N, N, L), dtype=np.uint8)
    result[:, :, 0] = 1
    for group in groups:
        mask = np.zeros((N, N), dtype=np.uint16)
        for i, vertex in enumerate(LEFT):
            member = np.isin(COLUMNS[:, i], list(group))
            mask |= member.astype(np.uint16)[:, None] << vertex
        for i, vertex in enumerate(RIGHT):
            member = np.isin(COLUMNS[:, i], list(group))
            mask |= member.astype(np.uint16)[None, :] << vertex
        factor = EVALUATED[mask]
        updated = np.zeros_like(result)
        for i in range(L):
            for j in range(L - i):
                updated[:, :, i + j] ^= MULTIPLY[
                    result[:, :, i], factor[:, :, j]]
        result = updated
    return [result[:, :, k].copy() for k in range(L)]


def pair_weight_blocks(weights):
    blocks = group_blocks([set(range(1, Q))])
    mismatch = (COLUMNS[:, None, :] != COLUMNS[None, :, :]).sum(axis=2)
    pair_weight = np.zeros_like(mismatch, dtype=np.uint8)
    power = 1
    for count in range(5):
        pair_weight[mismatch == count] = power
        power = MULTIPLY[power, weights]
    return [MULTIPLY[pair_weight, block] for block in blocks]


def color_average_blocks(weights):
    base = color_matrix(Q, Y, Z, L)
    blocks = [np.zeros_like(base[0]) for _ in range(L)]
    for c, weight in enumerate(weights):
        permutation = SHIFT_INDEX[c]
        for k in range(L):
            blocks[k] ^= MULTIPLY[
                weight, base[k][np.ix_(permutation, permutation)]]
    return blocks


def support_diagonal_profile():
    profile = Counter()
    for word in COLUMNS:
        mask = 0
        for vertex, color in zip(LEFT, word):
            if color:
                mask |= 1 << vertex
        for vertex, color in zip(RIGHT, word):
            if color:
                mask |= 1 << vertex
        polynomial = evaluate(GAUSSIAN[mask], Y, Z, L)
        nonzero = np.nonzero(polynomial)[0]
        assert nonzero.size
        profile[int(nonzero[0])] += 1
    return dict(sorted(profile.items()))


def grouped_diagonal_profile(groups):
    profile = Counter()
    for word in COLUMNS:
        total = 0
        for group in groups:
            mask = 0
            for vertex, color in zip(LEFT, word):
                if color in group:
                    mask |= 1 << vertex
            for vertex, color in zip(RIGHT, word):
                if color in group:
                    mask |= 1 << vertex
            polynomial = evaluate(GAUSSIAN[mask], Y, Z, L)
            nonzero = np.nonzero(polynomial)[0]
            assert nonzero.size
            total += int(nonzero[0])
        profile[total] += 1
    return dict(sorted(profile.items()))


def normalize(weights):
    total = 0
    for weight in weights:
        total ^= weight
    assert total
    inverse = int(np.nonzero(MULTIPLY[total] == 1)[0][0])
    return tuple(int(MULTIPLY[inverse, weight]) for weight in weights)


def pressure_floor(profile):
    orders = [order for order, count in sorted(profile.items())
              for _ in range(count)]
    best = (10.0, None)
    for i in range(1, 240):
        t = i / 240
        log_z = log2(sum(t ** order for order in orders))
        for j in range(1, 240):
            p = j / 240
            occupied = 1 - p
            entropy = (-p * log2(p) if p else 0.0)
            entropy -= (1 - p) * log2(1 - p) if p < 1 else 0.0
            entropy += occupied * log2(Q - 1)
            if entropy <= 0:
                continue
            value = (log_z - 8 * occupied * log2(t)) / (4 * entropy)
            if value < best[0]:
                best = (value, (t, p))
    return best


def endpoint_spurious_count(groups):
    """Full colourings with every group mask antipodal-closed but not
    individually antipodal-equal."""
    pairs = [(v, v ^ 7) for v in range(4)]
    count = 0
    for coloring in product(range(Q), repeat=8):
        individual = all(coloring[i] == coloring[j] for i, j in pairs)
        if individual:
            continue
        coarse = True
        for group in groups:
            for i, j in pairs:
                if (coloring[i] in group) != (coloring[j] in group):
                    coarse = False
                    break
            if not coarse:
                break
        if coarse:
            count += 1
    return count


def main():
    one_hot = specialized_profile(Q, Y, Z, L - 1)[0]
    support = block_profile(group_blocks([set(range(1, Q))]))
    grouped = block_profile(group_blocks([{1, 2}, {3}]))
    pair_two = block_profile(pair_weight_blocks(2))
    average = block_profile(color_average_blocks(
        normalize((203, 141, 70, 1))))
    support_diagonal = support_diagonal_profile()
    grouped_diagonal = grouped_diagonal_profile([{1, 2}, {3}])

    assert one_hot == {0: 1, 2: 3, 4: 6, 6: 10}
    assert support == {0: 1, 2: 1, 4: 1, 6: 1}
    assert grouped == {0: 1, 2: 2, 4: 3, 6: 4}
    assert pair_two == {0: 136, 2: 120}
    assert average == {0: 4, 2: 12, 4: 24, 6: 32}
    assert support_diagonal == {0: 136, 2: 120}
    assert grouped_diagonal == {0: 72, 2: 120, 4: 64}
    assert endpoint_spurious_count([set(range(1, Q))]) > 0
    assert endpoint_spurious_count([{1, 2}, {3}]) == 1040

    receipt = json.loads(Path("evidence/color-class-cube-profile.json").read_text())
    full_one_hot = {int(order): count
                    for order, count in receipt["q4_full_profile"]}
    full_floor = pressure_floor(full_one_hot)
    assert full_floor[0] < 0.7

    records = {
        "one_hot_through6": one_hot,
        "support_through6": support,
        "grouped_two_through6": grouped,
        "pair_weight_two_through6": pair_two,
        "color_average_through6": average,
        "support_diagonal_through6": support_diagonal,
        "grouped_diagonal_through6": grouped_diagonal,
        "one_hot_full_pressure_floor": list(full_floor[0:1]),
        "low_order_pressure_floors": {
            "one_hot": list(pressure_floor(one_hot)[0:1]),
            "support": list(pressure_floor(support)[0:1]),
            "grouped_two": list(pressure_floor(grouped)[0:1]),
            "pair_weight_two": list(pressure_floor(pair_two)[0:1]),
            "color_average": list(pressure_floor(average)[0:1]),
            "support_diagonal": list(pressure_floor(support_diagonal)[0:1]),
            "grouped_diagonal": list(pressure_floor(grouped_diagonal)[0:1]),
        },
        "finite_order_counts_through6": {
            "one_hot": sum(one_hot.values()),
            "support": sum(support.values()),
            "grouped_two": sum(grouped.values()),
            "pair_weight_two": sum(pair_two.values()),
            "color_average": sum(average.values()),
            "support_diagonal": sum(support_diagonal.values()),
            "grouped_diagonal": sum(grouped_diagonal.values()),
        },
        "zero_orders_through6": {
            "one_hot": one_hot.get(0, 0),
            "support": support.get(0, 0),
            "grouped_two": grouped.get(0, 0),
            "pair_weight_two": pair_two.get(0, 0),
            "color_average": average.get(0, 0),
            "support_diagonal": support_diagonal.get(0, 0),
            "grouped_diagonal": grouped_diagonal.get(0, 0),
        },
        "endpoint_spurious_counts_q4": {
            "support": endpoint_spurious_count([set(range(1, Q))]),
            "grouped_two": endpoint_spurious_count([{1, 2}, {3}]),
        },
        "scope": "q=4 finite controls through order six; no uniform theorem",
        "status": "passed",
    }
    print(json.dumps(records, sort_keys=True))
    return records


if __name__ == "__main__":
    main()
