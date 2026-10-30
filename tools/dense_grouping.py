"""Exact helpers for the dense-width-graph grouping theorem.

These helpers check finite certificates. They are not a factoring evaluator
and do not establish the asymptotic theorem by testing.

The two parameter choices and their mathematical proofs are documented in
research/dense-grouping.md and research/bounded-edge-grouping.md. Their
asymptotically degree-two factoring schedules remain conditional. These
helpers also do not check the Appendix B alpha=1/2 machine bridge. The
separate certified k6 alpha=1/30 theorem is documented in
research/certified-k6.md and proved in variants/k6/FactoringK6/Certified.lean.
"""

from __future__ import annotations

from itertools import product
from typing import Sequence


def dense_group_size(groups: int, accuracy: int) -> int:
    if groups < 2 or accuracy < 2:
        raise ValueError("groups and accuracy must both be at least 2")
    length = groups.bit_length() + accuracy.bit_length() + 8
    return 576 * accuracy**2 * groups * length


def bounded_edge_group_size(groups: int, accuracy: int) -> int:
    """The stronger grid-relevant theorem, using degree bound g*B."""
    if groups < 2 or accuracy < 2:
        raise ValueError("groups and accuracy must both be at least 2")
    length = groups.bit_length() + accuracy.bit_length() + 8
    return 144 * accuracy**2 * groups * length


def validate_weights(weights: Sequence[Sequence[int]]) -> tuple[list[int], int]:
    size = len(weights)
    if not size or any(len(row) != size for row in weights):
        raise ValueError("weights must be a nonempty square matrix")
    for u, row in enumerate(weights):
        for v, value in enumerate(row):
            if not isinstance(value, int) or value < 0:
                raise ValueError("weights must be nonnegative integers")
            if u == v and value:
                raise ValueError("diagonal weights must be zero")
            if value != weights[v][u]:
                raise ValueError("weights must be symmetric")
    degrees = [sum(row) for row in weights]
    return degrees, sum(degrees)


def has_dense_widths(weights: Sequence[Sequence[int]], bound: int) -> bool:
    validate_weights(weights)
    if bound <= 0:
        raise ValueError("positive-width density needs a positive bound")
    return all(
        2 * bound <= 3 * weights[u][v] <= 3 * bound
        for u in range(len(weights))
        for v in range(u + 1, len(weights))
    )


def grouping_statistics(
    weights: Sequence[Sequence[int]], coloring: Sequence[int], groups: int
) -> tuple[list[int], list[list[int]], list[list[int]], list[int], int]:
    degrees, total = validate_weights(weights)
    if groups < 2 or len(coloring) != len(weights):
        raise ValueError("invalid coloring dimensions")
    if any(not isinstance(a, int) or not 0 <= a < groups for a in coloring):
        raise ValueError("color labels must be integers in range(groups)")
    sizes = [0] * groups
    degree_sums = [0] * groups
    vertex_sums = [[0] * groups for _ in weights]
    pairs = [[0] * groups for _ in range(groups)]
    for u, a in enumerate(coloring):
        sizes[a] += 1
        degree_sums[a] += degrees[u]
        for v, b in enumerate(coloring):
            vertex_sums[u][b] += weights[u][v]
            # For a != b, this directed sum counts each undirected
            # edge once in pairs[a][b] and once in pairs[b][a].
            if a != b:
                pairs[a][b] += weights[u][v]
    return sizes, pairs, vertex_sums, degree_sums, total


def good_grouping(
    weights: Sequence[Sequence[int]],
    coloring: Sequence[int],
    groups: int,
    accuracy: int,
    degree_bound: int | None = None,
) -> bool:
    if accuracy < 2:
        raise ValueError("accuracy must be at least 2")
    degrees, total = validate_weights(weights)
    bound = max(degrees) if degree_bound is None else degree_bound
    if bound < max(degrees):
        raise ValueError("degree_bound is smaller than an actual degree")
    sizes, pairs, _, _, _ = grouping_statistics(weights, coloring, groups)
    size = len(weights)
    return all(groups * value <= 2 * size for value in sizes) and all(
        accuracy * groups**2 * pairs[a][b] <= accuracy * total + size * bound
        for a in range(groups)
        for b in range(a + 1, groups)
    )


def strong_balance_certificate(
    weights: Sequence[Sequence[int]],
    coloring: Sequence[int],
    groups: int,
    accuracy: int,
) -> bool:
    """The exact intermediate inequalities used in the existence proof."""
    if accuracy < 2:
        raise ValueError("accuracy must be at least 2")
    degrees, _ = validate_weights(weights)
    sizes, _, vertex_sums, degree_sums, total = grouping_statistics(
        weights, coloring, groups
    )
    scale = 4 * accuracy * groups
    return (
        all(groups * value <= 2 * len(weights) for value in sizes)
        and all(
            scale * vertex_sums[u][b] <= (4 * accuracy + 1) * degrees[u]
            for u in range(len(weights))
            for b in range(groups)
        )
        and all(scale * value <= (4 * accuracy + 1) * total for value in degree_sums)
    )


def first_good_grouping(
    weights: Sequence[Sequence[int]], groups: int, accuracy: int
) -> tuple[int, ...] | None:
    """Exact finite search; exhaustive and intended only for small examples."""
    validate_weights(weights)
    for coloring in product(range(groups), repeat=len(weights)):
        if good_grouping(weights, coloring, groups, accuracy):
            return coloring
    return None


def additive_balance_certificate(
    weights: Sequence[Sequence[int]],
    coloring: Sequence[int],
    groups: int,
    accuracy: int,
    edge_bound: int,
) -> bool:
    """Intermediate balances for the bounded-edge form; all tests are integral."""
    if accuracy < 2 or edge_bound < 0:
        raise ValueError("accuracy must be at least 2 and edge_bound nonnegative")
    degrees, _ = validate_weights(weights)
    if any(value > edge_bound for row in weights for value in row):
        raise ValueError("edge_bound is smaller than an actual edge weight")
    sizes, _, vertex_sums, degree_sums, total = grouping_statistics(
        weights, coloring, groups
    )
    size = len(weights)
    scale = 4 * accuracy * groups
    return (
        all(groups * value <= 2 * size for value in sizes)
        and all(
            scale * vertex_sums[u][b] <= 4 * accuracy * degrees[u] + size * edge_bound
            for u in range(size)
            for b in range(groups)
        )
        and all(
            scale * value <= 4 * accuracy * total + size**2 * edge_bound
            for value in degree_sums
        )
    )


def union_budget_certificate(groups: int, accuracy: int) -> dict[str, int | bool]:
    """Integer inequalities sufficient for the analytic union bound.

    The theorem proves failure <= events * exp(-exponent), with
    exponent >= 4*length and events < 2**(3*length-5) < exp(3*length).
    This checks the integer premises; it never approximates an exponential.
    """
    size = dense_group_size(groups, accuracy)
    length = groups.bit_length() + accuracy.bit_length() + 8
    events = (size + 2) * groups
    tail_scale = 72 * accuracy**2 * groups
    return {
        "groups": groups,
        "accuracy": accuracy,
        "length": length,
        "vertices": size,
        "events": events,
        "tail_exponent_ge_4L": size - 1 >= 4 * length * tail_scale,
        "events_lt_two_pow_3L_minus_5": events < 2 ** (3 * length - 5),
    }


def bounded_edge_union_budget(groups: int, accuracy: int) -> dict[str, int | bool]:
    size = bounded_edge_group_size(groups, accuracy)
    length = groups.bit_length() + accuracy.bit_length() + 8
    events = (size + 2) * groups
    tail_scale = 36 * accuracy**2 * groups
    return {
        "groups": groups,
        "accuracy": accuracy,
        "length": length,
        "vertices": size,
        "events": events,
        "tail_exponent_ge_4L": size >= 4 * length * tail_scale,
        "events_lt_two_pow_3L_minus_7": events < 2 ** (3 * length - 7),
    }


def dense_schedule(depth: int, initial: int = 48) -> list[tuple[int, int, int]]:
    if depth < 0 or initial < 48:
        raise ValueError("depth must be nonnegative and initial at least 48")
    result = []
    groups = initial
    for level in range(depth):
        accuracy = (level + 48) ** 2
        size = dense_group_size(groups, accuracy)
        result.append((groups, accuracy, size))
        groups = size**2
    return result
