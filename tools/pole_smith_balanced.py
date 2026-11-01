#!/usr/bin/env python3
"""A one-pass Smith histogram when every residue imbalance is at most one."""

from __future__ import annotations

from collections import Counter

from pole_smith_flow import residue_tree


def dyadic_balanced_padding(n: int, left: tuple[int, int], right: tuple[int, int]) -> tuple:
    """Pad two odd intervals to the same power-of-two number of nodes.

    The returned domains contain each odd class modulo the third return
    value once. Every new node exceeds N. Numerical interval lengths are
    not enumerated.
    """
    if n < 1 or n % 2 == 0:
        raise ValueError("N must be positive and odd")
    for first, count in (left, right):
        if first < 1 or first % 2 == 0 or count < 1 or first + 2 * (count - 1) > n:
            raise ValueError("Inputs must be nonempty odd progressions contained in [1,N]")
    count = 1 << ((max(left[1], right[1]) - 1).bit_length())
    modulus = 2 * count

    def pad(interval: tuple[int, int]) -> tuple[tuple[int, int], ...]:
        first, length = interval
        if length == count:
            return (interval,)
        tail = first + 2 * length
        multiple = max(0, (n - tail) // modulus + 1)
        start = tail + multiple * modulus
        assert start > n
        return (interval, (start, count - length))

    return pad(left), pad(right), modulus


def balanced_smith_histogram(n: int, rows: tuple[int, ...], cols: tuple[int, ...],
                             k: int | None = None) -> dict:
    if n < 1 or n % 2 == 0:
        raise ValueError("N must be positive and odd")
    if len(rows) != len(cols):
        raise ValueError("Balanced histogram requires a square matrix")
    if len(set(rows)) != len(rows) or len(set(cols)) != len(cols):
        raise ValueError("Require distinct node sets")
    if any(value < 1 or value % 2 == 0 for value in rows + cols):
        raise ValueError("Require positive odd nodes")
    if not rows:
        return {"histogram": [], "dimension": 0, "base_precision": 1,
                "k": 1 if k is None else k, "factor_count": 0,
                "overlap_counts": [], "tree_vertices": 0, "postorder_visits": 0}
    bound = max(abs(n - min(rows) * min(cols)), abs(n - max(rows) * max(cols)))
    h = max(1, bound.bit_length())
    if k is None:
        k = h
    if k < h:
        raise ValueError("Displacement must separate all nonzero product differences")
    tree = residue_tree(rows, cols, n, h)
    size = len(tree.parent)
    children: list[list[int]] = [[] for _ in range(size)]
    depth = [0] * size
    for vertex in range(1, size):
        parent = tree.parent[vertex]
        children[parent].append(vertex)
        depth[vertex] = depth[parent] + tree.length[vertex]
    row_count = tree.row_capacity[:]
    col_count = tree.col_capacity[:]
    histogram = Counter()
    matching_depths = Counter()
    for vertex in range(size - 1, -1, -1):
        row_count[vertex] += sum(row_count[child] for child in children[vertex])
        col_count[vertex] += sum(col_count[child] for child in children[vertex])
        if abs(row_count[vertex] - col_count[vertex]) > 1:
            raise ValueError("The residue tree has imbalance greater than one")
        matches = (min(row_count[vertex], col_count[vertex])
                   - sum(min(row_count[child], col_count[child]) for child in children[vertex]))
        assert matches >= 0
        if matches:
            base_order = h - depth[vertex]
            order = base_order + k - h if base_order else 0
            histogram[order] += matches
            matching_depths[depth[vertex]] += matches
    overlap = [0] * h
    count = 0
    for level in range(h, 0, -1):
        count += matching_depths[level]
        overlap[level - 1] = count
    assert sum(histogram.values()) == len(rows)
    assert overlap[0] == len(rows)
    return {"histogram": [{"order": order, "multiplicity": count}
                           for order, count in sorted(histogram.items())],
            "dimension": len(rows), "base_precision": h, "k": k,
            "factor_count": histogram.get(0, 0), "overlap_counts": overlap,
            "tree_vertices": size, "postorder_visits": size}
