"""Exact path-graph obstruction to subquadratic universal grouping.

The all-size proof is in research/grouping-obstruction.md. Exhaustive search
below the threshold and explicit witnesses at it independently check small
instances using the existing grouping verifier.
"""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.grouping import NoGroupingFound, WeightedGraph, find_grouping, verify_grouping


def path_graph(vertices: int) -> WeightedGraph:
    if not isinstance(vertices, int) or isinstance(vertices, bool) or vertices < 3:
        raise ValueError("path requires an integer vertex count at least 3")
    return WeightedGraph(vertices, tuple((u, u + 1, 1) for u in range(vertices - 1)))


def path_threshold(groups: int, discrepancy: int) -> int:
    """Exact minimum g >= 3 admitting a path grouping, for K >= 3, D >= 2."""
    if (not isinstance(groups, int) or isinstance(groups, bool) or groups < 3
            or not isinstance(discrepancy, int) or isinstance(discrepancy, bool) or discrepancy < 2):
        raise ValueError("require integer K >= 3 and integer D >= 2")
    numerator = discrepancy * (groups * groups + 2)
    denominator = 2 * (discrepancy + 1)
    return (numerator + denominator - 1) // denominator


def run_suite() -> dict:
    cases = []
    for groups, discrepancy in ((3, 2), (3, 3), (4, 2), (4, 10)):
        threshold = path_threshold(groups, discrepancy)
        below = path_graph(threshold - 1)
        margin = discrepancy * groups ** 2 - (discrepancy * below.ordered_total
                                                + below.vertex_count * below.max_degree)
        assert margin > 0
        try:
            find_grouping(below, groups, discrepancy, max_candidates=100_000)
        except NoGroupingFound:
            pass
        else:
            raise AssertionError("a grouping exists below the proved path threshold")
        graph = path_graph(threshold)
        block_size = (threshold + groups - 1) // groups
        witness = verify_grouping(graph, groups, discrepancy,
                                  tuple(u // block_size for u in range(threshold)))
        assert witness.valid
        cases.append({
            "K": groups, "D": discrepancy, "threshold": threshold,
            "excluded_vertices": below.vertex_count, "excluded_edges": below.edges,
            "ordered_total": below.ordered_total, "degree_bound": below.max_degree,
            "strict_obstruction_margin": margin,
            "exhausted_colorings": groups ** below.vertex_count,
            "boundary_witness": asdict(witness),
        })
    return {"scope": "finite checks of the path obstruction and its sharp boundary; not a factoring bound",
            "cases": cases}


if __name__ == "__main__":
    print(json.dumps(run_suite(), indent=2))
