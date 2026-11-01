#!/usr/bin/env python3
"""Disjointness determinants between equal-size partial matchings over F_2.

An edge is a two-bit nonnegative integer.  Within each matching the edges
are vertex-disjoint.  The component algorithm performs O(k) graph work;
it does not construct a k-by-k matrix or perform Gaussian elimination.
"""


def intersection_components(first, second):
    """Return (left size, right size, is_cycle) for each simple component."""
    if len(first) != len(second):
        raise ValueError("The two partial matchings must have equal sizes.")
    size = len(first)
    owners = {}
    for index, edge in enumerate(second):
        if edge < 0 or edge.bit_count() != 2:
            raise ValueError("Each edge must contain exactly two vertices.")
        remaining = edge
        while remaining:
            vertex = remaining & -remaining
            remaining ^= vertex
            if vertex in owners:
                raise ValueError("Edges within a matching must be disjoint.")
            owners[vertex] = size + index
    adjacency = [[] for _ in range(2 * size)]
    used = 0
    for index, edge in enumerate(first):
        if edge < 0 or edge.bit_count() != 2:
            raise ValueError("Each edge must contain exactly two vertices.")
        if used & edge:
            raise ValueError("Edges within a matching must be disjoint.")
        used |= edge
        neighbors = set()
        remaining = edge
        while remaining:
            vertex = remaining & -remaining
            remaining ^= vertex
            if vertex in owners:
                neighbors.add(owners[vertex])
        for neighbor in neighbors:
            adjacency[index].append(neighbor)
            adjacency[neighbor].append(index)
    seen = set()
    components = []
    for start in range(2 * size):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        left = right = 0
        cycle = True
        while stack:
            vertex = stack.pop()
            left += vertex < size
            right += vertex >= size
            cycle &= len(adjacency[vertex]) == 2
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append((left, right, cycle))
    return components


def disjointness_determinant(first, second, *, with_case=False):
    """Return det([e & f == 0]) modulo two by component classification."""
    components = intersection_components(first, second)
    cycles = []
    left_major = []
    right_major = []
    balanced = []
    for left, right, cycle in components:
        if cycle:
            assert left == right and left >= 2
            cycles.append(left)
        elif left == right:
            balanced.append(left)
        elif left == right + 1:
            left_major.append(left)
        else:
            assert right == left + 1
            right_major.append(right)
    assert len(left_major) == len(right_major)
    nullity = len(cycles) + len(left_major)
    if nullity >= 2:
        result, case = 0, "nullity_at_least_two"
    elif cycles:
        result, case = cycles[0] % 2, "one_cycle"
    elif left_major:
        result = left_major[0] * right_major[0] % 2
        case = "one_unbalanced_pair"
    else:
        result = (1 + sum((length + 1) // 2 for length in balanced)) % 2
        case = "balanced_paths_only"
    return (result, case) if with_case else result
