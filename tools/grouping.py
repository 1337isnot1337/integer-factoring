"""Exact, bounded certificate search for the generalized grouping lemma.

This is a finite research checker, not the proposed asymptotic factoring
algorithm. It uses only Python's standard library and integer arithmetic.
Repeated undirected edges are combined by adding their weights.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from itertools import product
import json
from typing import Iterable


Edge = tuple[int, int, int]


def _integer(name: str, value: int, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return value


@dataclass(frozen=True)
class WeightedGraph:
    """A loopless graph with canonical, combined nonzero integer weights."""

    vertex_count: int
    edges: tuple[Edge, ...]
    degrees: tuple[int, ...] = field(init=False)
    max_degree: int = field(init=False)
    ordered_total: int = field(init=False)

    def __post_init__(self) -> None:
        n = _integer("vertex_count", self.vertex_count)
        combined: dict[tuple[int, int], int] = {}
        for edge in self.edges:
            if not isinstance(edge, (tuple, list)) or len(edge) != 3:
                raise ValueError("each edge must be a three-item (u, v, weight) tuple")
            u, v, weight = edge
            _integer("edge endpoint", u)
            _integer("edge endpoint", v)
            _integer("edge weight", weight)
            if u >= n or v >= n:
                raise ValueError("edge endpoint is outside the graph")
            if u == v:
                raise ValueError("self loops are not permitted, including zero-weight loops")
            if weight:
                pair = (min(u, v), max(u, v))
                combined[pair] = combined.get(pair, 0) + weight
        canonical = tuple((u, v, weight) for (u, v), weight in sorted(combined.items()))
        degrees = [0] * n
        for u, v, weight in canonical:
            degrees[u] += weight
            degrees[v] += weight
        object.__setattr__(self, "edges", canonical)
        object.__setattr__(self, "degrees", tuple(degrees))
        object.__setattr__(self, "max_degree", max(degrees, default=0))
        object.__setattr__(self, "ordered_total", sum(degrees))


def sufficient_node_count(groups: int, discrepancy: int) -> int:
    """Integer-only sufficient threshold, stronger than 12 D^2 K^2 ln(K)."""

    k = _integer("groups", groups, 2)
    d = _integer("discrepancy", discrepancy, 2)
    return 12 * d * d * k * k * k.bit_length()


@dataclass(frozen=True)
class GroupingCertificate:
    vertex_count: int
    groups: int
    discrepancy: int
    degree_bound: int
    ordered_total: int
    colors: tuple[int, ...]
    # Sparse records omit zero-sized classes and zero-weight pairs.
    group_sizes: tuple[tuple[int, int], ...]
    pair_weights: tuple[Edge, ...]
    sufficient_threshold_met: bool
    valid: bool
    violations: tuple[str, ...]


def verify_grouping(
    graph: WeightedGraph,
    groups: int,
    discrepancy: int,
    colors: Iterable[int],
    *,
    degree_bound: int | None = None,
) -> GroupingCertificate:
    """Recompute a certificate from the graph; never trust supplied summaries.

    Structural input errors raise ValueError. A well-formed coloring that
    violates the two mathematical inequalities returns valid=False.
    """

    k = _integer("groups", groups, 2)
    d = _integer("discrepancy", discrepancy, 2)
    p = graph.max_degree if degree_bound is None else _integer("degree_bound", degree_bound)
    if p < graph.max_degree:
        raise ValueError("degree_bound is smaller than the actual maximum weighted degree")
    assignment = tuple(colors)
    if len(assignment) != graph.vertex_count:
        raise ValueError("the coloring must contain exactly one color per vertex")
    sizes: dict[int, int] = {}
    for color in assignment:
        _integer("color", color)
        if color >= k:
            raise ValueError("color is outside the group range")
        sizes[color] = sizes.get(color, 0) + 1
    pairs: dict[tuple[int, int], int] = {}
    for u, v, weight in graph.edges:
        a, b = assignment[u], assignment[v]
        if a != b:
            pair = (min(a, b), max(a, b))
            pairs[pair] = pairs.get(pair, 0) + weight
    violations: list[str] = []
    for color, size in sorted(sizes.items()):
        if k * size > 2 * graph.vertex_count:
            violations.append(f"group {color}: K * size = {k * size} > 2g = {2 * graph.vertex_count}")
    pair_bound = d * graph.ordered_total + graph.vertex_count * p
    for (a, b), weight in sorted(pairs.items()):
        if d * k * k * weight > pair_bound:
            violations.append(
                f"pair ({a}, {b}): D * K^2 * Y = {d * k * k * weight} > D * U + g * P = {pair_bound}"
            )
    return GroupingCertificate(
        vertex_count=graph.vertex_count,
        groups=k,
        discrepancy=d,
        degree_bound=p,
        ordered_total=graph.ordered_total,
        colors=assignment,
        group_sizes=tuple(sorted(sizes.items())),
        pair_weights=tuple((a, b, weight) for (a, b), weight in sorted(pairs.items())),
        sufficient_threshold_met=graph.vertex_count >= sufficient_node_count(k, d),
        valid=not violations,
        violations=tuple(violations),
    )


class SearchLimitExceeded(RuntimeError):
    """The requested exhaustive search exceeds its explicit candidate cap."""


class NoGroupingFound(RuntimeError):
    """Every coloring was checked and none meets the requested inequalities."""


def _bounded_power(base: int, exponent: int, cap: int) -> int:
    """Return base**exponent if <= cap, otherwise cap+1, without huge integers."""

    result = 1
    factor = min(base, cap + 1)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = min(result * factor, cap + 1)
            if result > cap:
                return cap + 1
        remaining >>= 1
        if remaining:
            factor = min(factor * factor, cap + 1)
    return result


@dataclass(frozen=True)
class SearchResult:
    certificate: GroupingCertificate
    candidates_examined: int
    total_colorings: int
    candidate_cap: int


def find_grouping(
    graph: WeightedGraph,
    groups: int,
    discrepancy: int,
    *,
    degree_bound: int | None = None,
    max_candidates: int = 100_000,
) -> SearchResult:
    """Find the lexicographically first certificate, subject to a preflight cap.

    The cap applies to the complete K**g search space, even when an early
    valid coloring might exist. A failed cap never starts the enumeration.
    If the space fits the cap but no coloring works, NoGroupingFound is
    evidence of finite nonexistence for those parameters.
    """

    k = _integer("groups", groups, 2)
    d = _integer("discrepancy", discrepancy, 2)
    cap = _integer("max_candidates", max_candidates, 1)
    p = graph.max_degree if degree_bound is None else _integer("degree_bound", degree_bound)
    if p < graph.max_degree:
        raise ValueError("degree_bound is smaller than the actual maximum weighted degree")
    total = _bounded_power(k, graph.vertex_count, cap)
    if total > cap:
        raise SearchLimitExceeded(
            f"K^g exceeds the explicit cap of {cap} candidates (K={k}, g={graph.vertex_count}); "
            "no colorings were enumerated"
        )
    if graph.vertex_count == 0:
        # itertools.product materializes its input pools even when repeat=0.
        # Avoid constructing range(k) for an empty graph with an enormous k.
        return SearchResult(verify_grouping(graph, k, d, (), degree_bound=p), 1, 1, cap)
    for examined, colors in enumerate(product(range(k), repeat=graph.vertex_count), 1):
        certificate = verify_grouping(graph, k, d, colors, degree_bound=p)
        if certificate.valid:
            return SearchResult(certificate, examined, total, cap)
    raise NoGroupingFound(f"none of the {total} colorings satisfies the grouping inequalities")


def demonstration() -> dict:
    """Small exact examples, including a certified nonexistence case."""

    examples = {
        "matching": WeightedGraph(6, ((0, 1, 5), (2, 3, 5), (4, 5, 5))),
        "star": WeightedGraph(5, tuple((0, v, v) for v in range(1, 5))),
        "duplicate_edges": WeightedGraph(4, ((0, 1, 2), (1, 0, 3), (2, 3, 4), (0, 2, 0))),
        "zero_graph": WeightedGraph(4, ()),
    }
    results = {
        name: {"edges": graph.edges, "result": asdict(find_grouping(graph, 3, 2))}
        for name, graph in examples.items()
    }
    try:
        find_grouping(WeightedGraph(2, ((0, 1, 1),)), 3, 2)
    except NoGroupingFound as error:
        results["nonexistence"] = {"vertices": 2, "edges": ((0, 1, 1),), "K": 3, "D": 2, "result": str(error)}
    return {
        "purpose": "finite exact certificates; these examples do not prove the asymptotic theorem",
        "integer_sufficient_threshold_K3_D2": sufficient_node_count(3, 2),
        "examples": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo", action="store_true", help="print the reproducible finite certificate examples")
    args = parser.parse_args()
    if not args.demo:
        parser.error("pass --demo; import find_grouping for explicitly bounded custom searches")
    print(json.dumps(demonstration(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
