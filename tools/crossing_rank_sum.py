#!/usr/bin/env python3
"""Exact rank-sum optimization in a finite F2 crossing coefficient span.

Tensor entries are packed in row-major order with row (a,y), column
(b,u). The second cut has row (a,u), column (b,y). Thus the default
crossing [a=b][u=y] is the identity matrix in either cut.
"""

from dataclasses import dataclass
import heapq
from math import isqrt


def binary_rank(rows):
    basis = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def unpack_tensor(tensor, labels):
    dimension = labels * labels
    mask = (1 << dimension) - 1
    return tuple((tensor >> (dimension * row)) & mask for row in range(dimension))


def pack_tensor(rows):
    dimension = len(rows)
    if any(row < 0 or row.bit_length() > dimension for row in rows):
        raise ValueError("Rows must encode a square binary matrix.")
    return sum(row << (dimension * i) for i, row in enumerate(rows))


def other_cut(rows, labels):
    result = [0] * (labels * labels)
    for row, value in enumerate(rows):
        a, y = divmod(row, labels)
        while value:
            bit = value & -value
            value ^= bit
            b, u = divmod(bit.bit_length() - 1, labels)
            result[a * labels + u] ^= 1 << (b * labels + y)
    return tuple(result)


def crossing_tensor(labels):
    dimension = labels * labels
    return sum(1 << (dimension * i + i) for i in range(dimension))


def independent_generators(tensors):
    basis = {}
    generators = []
    for value in tensors:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot][0]
            else:
                basis[pivot] = value, len(generators)
                generators.append(value)
                break
    return tuple(generators), basis


def span_coordinates(tensor, basis):
    result = 0
    while tensor:
        pivot = tensor.bit_length() - 1
        if pivot not in basis:
            return None
        value, index = basis[pivot]
        tensor ^= value
        result ^= 1 << index
    return result


def xor_shortest_path(costs, target):
    """Exact integer shortest path; edge v moves x to x XOR v."""
    if not costs or costs[0] != 0 or len(costs) & (len(costs) - 1):
        raise ValueError("Costs must index a binary vector space, with cost(0)=0.")
    if any(cost <= 0 for cost in costs[1:]):
        raise ValueError("Every nonzero vector must have positive cost.")
    infinity = sum(costs) + 1
    distance = [infinity] * len(costs)
    distance[0] = 0
    parents = [None] * len(costs)
    heap = [(0, 0)]
    edges = sorted((cost, vector) for vector, cost in enumerate(costs) if vector)
    while heap:
        cost, node = heapq.heappop(heap)
        if cost != distance[node]:
            continue
        if node == target:
            witness = []
            while node:
                previous, edge = parents[node]
                witness.append(edge)
                node = previous
            return cost, tuple(witness)
        for edge_cost, edge in edges:
            candidate = cost + edge_cost
            if candidate >= distance[target]:
                break
            neighbor = node ^ edge
            if candidate < distance[neighbor]:
                distance[neighbor] = candidate
                parents[neighbor] = node, edge
                heapq.heappush(heap, (candidate, neighbor))
    raise AssertionError("A nonzero target is itself an available edge.")


@dataclass(frozen=True)
class RankCost:
    row_rank: int
    column_rank: int
    lower: int
    upper: int


class RankOracle:
    def __init__(self, labels, scale=1_000_000, *, coordinate_generators=None):
        if labels < 2 or scale < 1:
            raise ValueError("At least two labels and a positive integer scale are required.")
        self.labels = labels
        self.scale = scale
        self.coordinate_basis = None
        if coordinate_generators is not None:
            _, self.coordinate_basis = independent_generators(coordinate_generators)
        self.cache = {0: RankCost(0, 0, 0, 0)}

    def __call__(self, tensor):
        if tensor < 0 or tensor.bit_length() > self.labels**4:
            raise ValueError("The packed tensor has the wrong size.")
        key = tensor
        if self.coordinate_basis is not None:
            key = span_coordinates(tensor, self.coordinate_basis)
            if key is None:
                raise ValueError("The tensor is outside the oracle's declared coordinate span.")
        if key not in self.cache:
            rows = unpack_tensor(tensor, self.labels)
            row_rank = binary_rank(rows)
            column_rank = binary_rank(other_cut(rows, self.labels))
            square = row_rank * column_rank * self.scale**2
            lower = isqrt(square)
            upper = lower + int(lower**2 < square)
            self.cache[key] = RankCost(row_rank, column_rank, lower, upper)
        return self.cache[key]


def optimize_span(tensors, oracle, *, target=None, maximum_dimension=16):
    """Enclose the least sum sqrt(row_rank*column_rank) in the full span.

    The lower bound uses floor costs for every nonzero span vector.
    The upper bound is the cost of a concrete XOR decomposition using
    ceiling costs. Both bounds are integers divided by oracle.scale.
    """
    if target is None:
        target = crossing_tensor(oracle.labels)
    generators, basis = independent_generators(tensors)
    if len(generators) > maximum_dimension:
        raise ValueError("The coefficient span exceeds the declared exhaustive-search limit.")
    target_index = span_coordinates(target, basis)
    if target_index is None:
        raise ValueError("The target tensor does not lie in the supplied span.")
    values = [0] * (1 << len(generators))
    lower = [0] * len(values)
    upper = [0] * len(values)
    for index in range(1, len(values)):
        bit = index & -index
        values[index] = values[index ^ bit] ^ generators[bit.bit_length() - 1]
        bound = oracle(values[index])
        lower[index], upper[index] = bound.lower, bound.upper
    lower_value, _ = xor_shortest_path(lower, target_index)
    upper_value, witness = xor_shortest_path(upper, target_index)
    reconstruction = 0
    entries = []
    for index in witness:
        reconstruction ^= values[index]
        bound = oracle(values[index])
        entries.append({
            "span_vector": index,
            "row_rank": bound.row_rank,
            "column_rank": bound.column_rank,
            "lower_scaled": bound.lower,
            "upper_scaled": bound.upper,
        })
    assert reconstruction == target
    return {
        "labels": oracle.labels,
        "scale": oracle.scale,
        "span_dimension": len(generators),
        "target_coordinates": target_index,
        "lower_scaled": lower_value,
        "upper_scaled": upper_value,
        "witness": entries,
    }


def certifies_upper_power(result, numerator, denominator):
    """Certify S < labels^(numerator/denominator), using integer powers."""
    if numerator < 0 or denominator < 1:
        raise ValueError("A nonnegative rational power is required.")
    return result["upper_scaled"]**denominator < (
        result["scale"]**denominator * result["labels"]**numerator)


def certifies_lower_power(result, numerator, denominator):
    """Certify S >= labels^(numerator/denominator), using integer powers."""
    if numerator < 0 or denominator < 1:
        raise ValueError("A nonnegative rational power is required.")
    return result["lower_scaled"]**denominator >= (
        result["scale"]**denominator * result["labels"]**numerator)


def pair_vectors(first, second, old_labels):
    """Kronecker columns for two binary label embeddings."""
    result = []
    for first_vector in first:
        for second_vector in second:
            value = 0
            remaining = first_vector
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                value ^= second_vector << (old_labels * (bit.bit_length() - 1))
            result.append(value)
    return tuple(result)


def dual_pair(first, second):
    return len(first) == len(second) and all(
        (a & b).bit_count() % 2 == int(i == j)
        for i, a in enumerate(first) for j, b in enumerate(second))


def restrict_tensor(tensor, old_labels, a, b, u, y):
    """Apply four possibly different, d-column binary label embeddings.

    The crossing is preserved when a.T*b=I and u.T*y=I.
    This function performs the linear restriction; the caller may also
    use it when those conditions are not desired.
    """
    if len({len(a), len(b), len(u), len(y)}) != 1:
        raise ValueError("The four new label spaces must have the same dimension.")
    if any(v < 0 or v.bit_length() > old_labels for family in (a, b, u, y) for v in family):
        raise ValueError("A label embedding has coordinates outside the old space.")
    rows = unpack_tensor(tensor, old_labels)
    left = pair_vectors(a, y, old_labels)
    right = pair_vectors(b, u, old_labels)
    result = []
    for vector in left:
        value = 0
        remaining = vector
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            value ^= rows[bit.bit_length() - 1]
        result.append(sum(((value & column).bit_count() % 2) << j
                          for j, column in enumerate(right)))
    return pack_tensor(result)
