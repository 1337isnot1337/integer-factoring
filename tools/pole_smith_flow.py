#!/usr/bin/env python3
"""Exact Cauchy pole Smith orders by convex flow on a residue trie.

Input lists are explicit. Complexity is polynomial in their numerical lengths;
this does not supply a compressed-interval factoring algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass
from bisect import bisect_left
import heapq


@dataclass
class ResidueTree:
    parent: list[int]
    length: list[int]
    row_capacity: list[int]
    col_capacity: list[int]
    uncompressed_nodes: int


def residue_tree(rows: tuple[int, ...], cols: tuple[int, ...], n: int, k: int) -> ResidueTree:
    """Build the compressed LSB-first trie directly, without its long paths."""
    if k < 1 or n % 2 == 0 or any(value % 2 == 0 for value in rows + cols):
        raise ValueError("N, row nodes, and column nodes must be odd, and k positive")
    modulus = 1 << k
    capacities: dict[int, list[int]] = {}
    for color, points in enumerate((rows, tuple(n * pow(b, -1, modulus) % modulus for b in cols))):
        for point in points:
            residue = point % modulus
            if residue not in capacities:
                capacities[residue] = [0, 0]
            capacities[residue][color] += 1
    reversed_items = sorted((int(format(value, f"0{k}b")[::-1], 2), value)
                            for value in capacities)
    keys = [item[0] for item in reversed_items]
    values = [item[1] for item in reversed_items]
    parent = [-1]
    depth = [0]
    row_capacity = [0]
    col_capacity = [0]
    stack = [(0, len(values), 0)]
    while stack:
        lower, upper, ancestor = stack.pop()
        if upper - lower == 1:
            parent.append(ancestor)
            depth.append(k)
            row_capacity.append(capacities[values[lower]][0])
            col_capacity.append(capacities[values[lower]][1])
            continue
        difference = values[lower] ^ values[upper - 1]
        common = (difference & -difference).bit_length() - 1
        if common == depth[ancestor]:
            vertex = ancestor
        else:
            vertex = len(parent)
            parent.append(ancestor)
            depth.append(common)
            row_capacity.append(0)
            col_capacity.append(0)
        prefix = (keys[lower] >> (k - common)) << (k - common)
        boundary = prefix + (1 << (k - common - 1))
        split = bisect_left(keys, boundary, lower, upper)
        assert lower < split < upper
        stack.append((split, upper, vertex))
        stack.append((lower, split, vertex))
    length = [0] + [depth[v] - depth[parent[v]] for v in range(1, len(parent))]
    assert len(parent) <= 2 * (len(rows) + len(cols))
    return ResidueTree(parent, length, row_capacity, col_capacity, 1 + sum(length))


def smith_profile_heap(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int | None = None) -> dict:
    """Return all Smith orders of 2^k*w/(N+2^k*w-a*b) for odd nodes.

    The only large list in the algorithm is the explicitly supplied node sets.
    No pair matrix, matrix entries, rational fractions, or minors are built.
    """
    if n < 1 or n % 2 == 0:
        raise ValueError("N must be a positive odd integer")
    if len(set(rows)) != len(rows) or len(set(cols)) != len(cols):
        raise ValueError("Node sets must be distinct")
    if any(value < 1 or value % 2 == 0 for value in rows + cols):
        raise ValueError("Node sets must contain positive odd integers")
    if not rows or not cols:
        return {"orders": [], "contents": [0], "k": 1 if k is None else k,
                "compressed_nodes": 0, "uncompressed_nodes": 0,
                "augmentations": 0}
    max_delta = max(abs(n - min(rows) * min(cols)), abs(n - max(rows) * max(cols)))
    h = max(1, max_delta.bit_length())
    if k is None:
        k = h
    if k < h:
        raise ValueError("k must be at least the bit length of every |N-a*b|")
    if k > h:
        result = smith_profile_heap(n, rows, cols, h)
        shift = k - h
        result["orders"] = [value + shift if value else 0 for value in result["orders"]]
        result["contents"] = [0]
        for value in result["orders"]:
            result["contents"].append(result["contents"][-1] + value)
        result["k"] = k
        result["base_precision"] = h
        return result
    tree = residue_tree(rows, cols, n, k)
    size = len(tree.parent)
    source, sink = size, size + 1
    adjacency: list[list[tuple[int, str, int, int]]] = [[] for _ in range(size + 2)]
    for child in range(1, size):
        ancestor = tree.parent[child]
        adjacency[child].append((ancestor, "tree", child, 1))
        adjacency[ancestor].append((child, "tree", child, -1))
    for vertex in range(size):
        if tree.row_capacity[vertex]:
            adjacency[source].append((vertex, "row", vertex, 1))
            adjacency[vertex].append((source, "row", vertex, -1))
        if tree.col_capacity[vertex]:
            adjacency[vertex].append((sink, "col", vertex, 1))
            adjacency[sink].append((vertex, "col", vertex, -1))
    flow = [0] * size
    rows_used = [0] * size
    cols_used = [0] * size
    potential = [0] * (size + 2)
    doubled_cost = 0
    contents = [0]
    orders = []
    for _ in range(min(len(rows), len(cols))):
        distance: list[int | None] = [None] * (size + 2)
        predecessor: list[tuple[int, str, int, int, int] | None] = [None] * (size + 2)
        distance[source] = 0
        queue = [(0, source)]
        while queue:
            current_distance, vertex = heapq.heappop(queue)
            if distance[vertex] != current_distance:
                continue
            for neighbor, kind, edge, direction in adjacency[vertex]:
                if kind == "tree":
                    cost = tree.length[edge] * (2 * direction * flow[edge] + 1)
                else:
                    used = rows_used if kind == "row" else cols_used
                    capacity = tree.row_capacity if kind == "row" else tree.col_capacity
                    if direction == 1 and used[edge] == capacity[edge]:
                        continue
                    if direction == -1 and used[edge] == 0:
                        continue
                    cost = 0
                reduced_cost = cost + potential[vertex] - potential[neighbor]
                if reduced_cost < 0:
                    raise AssertionError("Negative reduced cost invalidates the shortest-path invariant")
                candidate = current_distance + reduced_cost
                if distance[neighbor] is None or candidate < distance[neighbor]:
                    distance[neighbor] = candidate
                    predecessor[neighbor] = (vertex, kind, edge, direction, cost)
                    heapq.heappush(queue, (candidate, neighbor))
        if distance[sink] is None:
            raise AssertionError("The requested matching cardinality must be feasible")
        for vertex, delta in enumerate(distance):
            if delta is not None:
                potential[vertex] += delta
        vertex = sink
        while vertex != source:
            previous = predecessor[vertex]
            assert previous is not None
            ancestor, kind, edge, direction, cost = previous
            if kind == "tree":
                flow[edge] += direction
            elif kind == "row":
                rows_used[edge] += direction
            else:
                cols_used[edge] += direction
            doubled_cost += cost
            vertex = ancestor
        assert doubled_cost == sum(tree.length[v] * flow[v] ** 2 for v in range(1, size))
        assert doubled_cost % 2 == 0
        contents.append(doubled_cost // 2)
        orders.append(contents[-1] - contents[-2])
    assert orders == sorted(orders)
    return {"orders": orders, "contents": contents, "k": k,
            "base_precision": h,
            "compressed_nodes": size, "uncompressed_nodes": tree.uncompressed_nodes,
            "augmentations": len(orders)}


def smith_profile(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int | None = None) -> dict:
    """All Smith orders with updates confined to two residue-tree paths.

    With b-bit explicit input nodes, h=O(b) and each augmentation changes at
    most 2h+1 tree vertices. The implementation does not expand the matrix.
    """
    if n < 1 or n % 2 == 0:
        raise ValueError("N must be a positive odd integer")
    if len(set(rows)) != len(rows) or len(set(cols)) != len(cols):
        raise ValueError("Node sets must be distinct")
    if any(value < 1 or value % 2 == 0 for value in rows + cols):
        raise ValueError("Node sets must contain positive odd integers")
    if not rows or not cols:
        return {"orders": [], "contents": [0], "k": 1 if k is None else k,
                "base_precision": 1, "compressed_nodes": 0,
                "uncompressed_nodes": 0, "augmentations": 0,
                "vertex_updates": 0}
    max_delta = max(abs(n - min(rows) * min(cols)), abs(n - max(rows) * max(cols)))
    h = max(1, max_delta.bit_length())
    if k is None:
        k = h
    if k < h:
        raise ValueError("k must be at least the bit length of every |N-a*b|")
    tree = residue_tree(rows, cols, n, h)
    size = len(tree.parent)
    children: list[list[int]] = [[] for _ in range(size)]
    for child in range(1, size):
        children[tree.parent[child]].append(child)
    flow = [0] * size
    row_remaining = tree.row_capacity.copy()
    col_remaining = tree.col_capacity.copy()
    row_path: list[tuple[int, int] | None] = [None] * size
    col_path: list[tuple[int, int] | None] = [None] * size
    pair_path: list[tuple[int, int, int] | None] = [None] * size
    updates = 0

    def refresh(vertex: int) -> None:
        nonlocal updates
        updates += 1
        rr = [(0, vertex)] if row_remaining[vertex] else []
        cc = [(0, vertex)] if col_remaining[vertex] else []
        pairs = [(0, vertex, vertex)] if rr and cc else []
        child_rows = []
        child_cols = []
        for child in children[vertex]:
            row = row_path[child]
            col = col_path[child]
            pair = pair_path[child]
            if row is not None:
                item = (row[0] + tree.length[child] * (2 * flow[child] + 1), row[1])
                rr.append(item)
                child_rows.append((child, item))
            if col is not None:
                item = (col[0] + tree.length[child] * (-2 * flow[child] + 1), col[1])
                cc.append(item)
                child_cols.append((child, item))
            if pair is not None:
                pairs.append(pair)
        for row_child, row in child_rows:
            for col_child, col in child_cols:
                if row_child != col_child:
                    pairs.append((row[0] + col[0], row[1], col[1]))
        row_path[vertex] = min(rr) if rr else None
        col_path[vertex] = min(cc) if cc else None
        pair_path[vertex] = min(pairs) if pairs else None

    for vertex in range(size - 1, -1, -1):
        refresh(vertex)
    contents = [0]
    orders = []
    for _ in range(min(len(rows), len(cols))):
        best = pair_path[0]
        if best is None:
            raise AssertionError("A residual row-column path must exist")
        cost, row_leaf, col_leaf = best
        assert cost >= 0 and cost % 2 == 0
        orders.append(cost // 2)
        contents.append(contents[-1] + cost // 2)
        row_remaining[row_leaf] -= 1
        col_remaining[col_leaf] -= 1
        row_ancestors = []
        col_ancestors = []
        for leaf, path in ((row_leaf, row_ancestors), (col_leaf, col_ancestors)):
            vertex = leaf
            while vertex >= 0:
                path.append(vertex)
                vertex = tree.parent[vertex]
        shared = 0
        while (shared < min(len(row_ancestors), len(col_ancestors)) and
               row_ancestors[-1 - shared] == col_ancestors[-1 - shared]):
            shared += 1
        exclusive_rows = row_ancestors[:-shared]
        exclusive_cols = col_ancestors[:-shared]
        for vertex in exclusive_rows:
            flow[vertex] += 1
        for vertex in exclusive_cols:
            flow[vertex] -= 1
        for vertex in exclusive_rows + exclusive_cols + row_ancestors[-shared:]:
            refresh(vertex)
    assert orders == sorted(orders)
    assert 2 * contents[-1] == sum(tree.length[v] * flow[v] ** 2 for v in range(1, size))
    if k > h:
        orders = [value + k - h if value else 0 for value in orders]
        contents = [0]
        for value in orders:
            contents.append(contents[-1] + value)
    assert updates <= size + min(len(rows), len(cols)) * (2 * h + 1)
    return {"orders": orders, "contents": contents, "k": k,
            "base_precision": h, "compressed_nodes": size,
            "uncompressed_nodes": tree.uncompressed_nodes,
            "augmentations": len(orders), "vertex_updates": updates}
