#!/usr/bin/env python3
"""Verify component-path refinements of exact channel-walk bounds.

The graph example certifies a stronger counting method. It is not a
crossing tensor and supplies no improved factoring exponent.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from check_channel_flow_duals import (
    Graph, coupled_certificate, upper_beta, walk_types,
)


def component_path_cover(graph):
    """Return SCCs and induced subgraphs on maximal condensation paths."""
    n = graph.order
    reachable = [[s == t for t in range(n)] for s in range(n)]
    for s, t, _ in graph.edges:
        reachable[s][t] = True
    for middle in range(n):
        for s in range(n):
            for t in range(n):
                reachable[s][t] |= reachable[s][middle] and reachable[middle][t]
    unseen = set(range(n))
    components = []
    while unseen:
        first = min(unseen)
        part = tuple(t for t in sorted(unseen)
                     if reachable[first][t] and reachable[t][first])
        components.append(part)
        unseen.difference_update(part)
    owner = {v: i for i, part in enumerate(components) for v in part}
    outgoing = [set() for _ in components]
    incoming = [set() for _ in components]
    for s, t, _ in graph.edges:
        a, b = owner[s], owner[t]
        if a != b:
            outgoing[a].add(b)
            incoming[b].add(a)

    paths = []

    def extend(path):
        tail = path[-1]
        if not outgoing[tail]:
            paths.append(path)
        for following in sorted(outgoing[tail]):
            assert following not in path, "condensation graph has a cycle"
            extend(path + (following,))

    for start in range(len(components)):
        if not incoming[start]:
            extend((start,))
    subgraphs = []
    for path in paths:
        vertices = tuple(sorted(v for i in path for v in components[i]))
        remap = {v: i for i, v in enumerate(vertices)}
        edges = tuple((remap[s], remap[t], label)
                      for s, t, label in graph.edges if s in remap and t in remap)
        subgraphs.append(Graph(len(vertices), graph.pieces, edges))
    return tuple(components), tuple(paths), tuple(subgraphs)


def check_cover(graph, lengths):
    components, paths, subgraphs = component_path_cover(graph)
    comparisons = 0
    for length in lengths:
        original = walk_types(graph, length)
        covered = Counter()
        for subgraph in subgraphs:
            covered.update(walk_types(subgraph, length))
        assert all(covered[kind] >= count for kind, count in original.items())
        comparisons += len(original)
    return components, paths, subgraphs, comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    row = Graph(3, 3, ((0, 0, 0), (0, 0, 0),
                       (1, 1, 1), (1, 1, 1), (2, 2, 2)))
    col = Graph(3, 3, ((0, 1, 0), (0, 1, 0),
                       (1, 0, 1), (1, 0, 1), (2, 2, 2)))
    rc, rp, rows, comparisons_r = check_cover(row, range(13))
    cc, cp, cols, comparisons_c = check_cover(col, range(13))
    assert rc == ((0,), (1,), (2,)) and rp == ((0,), (1,), (2,))
    assert cc == ((0, 1), (2,)) and cp == ((0,), (1,))

    total, dual = coupled_certificate(row, dict.fromkeys(range(4), 1),
                                      col, dict.fromkeys(range(4), 1))
    assert total == 4 and dual == 2**8
    assert upper_beta(row, [F(1)] * 3, [F(1)] * 3) == 2
    assert upper_beta(col, [F(1)] * 3, [F(1)] * 3) == 2

    # row path, column path, piece weights, row beta, column beta, column h.
    certificates = (
        (0, 0, (F(1, 4), F(4), F(1)), F(1, 2), F(2), (F(4), F(1))),
        (1, 0, (F(4), F(1, 4), F(1)), F(1, 2), F(2), (F(1), F(4))),
        (2, 0, (F(4), F(4), F(1)), F(1), F(1, 2), (F(1), F(1))),
        (0, 1, (F(1, 4), F(1), F(1)), F(1, 2), F(1), (F(1),)),
        (1, 1, (F(1), F(1, 4), F(1)), F(1, 2), F(1), (F(1),)),
        (2, 1, (F(1), F(1), F(1)), F(1), F(1), (F(1),)),
    )
    report_certificates = []
    for r, c, weights, rb, cb, ch in certificates:
        assert upper_beta(rows[r], weights, [F(1)] * rows[r].order) == rb
        assert upper_beta(cols[c], [1 / z for z in weights], ch) == cb
        assert rb * cb <= 1
        report_certificates.append({
            "row_path": r, "column_path": c,
            "piece_weights": list(map(str, weights)),
            "row_beta": str(rb), "column_beta": str(cb),
            "column_positive_vector": list(map(str, ch)),
        })
    final_total, final_dual = coupled_certificate(rows[2], {0: 1}, cols[1], {0: 1})
    assert final_total == 1 and final_dual == 1
    for length in range(2, 13):
        rt, ct = walk_types(row, length), walk_types(col, length)
        common = rt.keys() & ct.keys()
        assert common == {(0, 0, length)}
        assert rt[0, 0, length] == ct[0, 0, length] == 1

    # A walk can move from one component to the next. Splitting them
    # separately would discard legitimate words and is not a cover.
    connected = Graph(2, 2, ((0, 0, 0), (0, 1, 1), (1, 1, 1)))
    nc, np, ng, comparisons_connected = check_cover(connected, range(8))
    assert nc == ((0,), (1,)) and np == ((0, 1),)
    assert ng == (connected,)
    invalid_split = (Graph(1, 2, ((0, 0, 0),)), Graph(1, 2, ((0, 0, 1),)))
    bad_count = sum((walk_types(g, 2)[1, 1] for g in invalid_split))
    assert bad_count == 0 < walk_types(connected, 2)[1, 1]

    # Branching paths, an internal SCC, and a skip edge check the full
    # induced-path convention rather than just disconnected components.
    branching = Graph(6, 2, ((0, 1, 0), (1, 0, 1), (1, 2, 0),
                             (0, 3, 1), (2, 4, 0), (3, 4, 1),
                             (0, 4, 0), (4, 5, 1), (5, 5, 0)))
    bc, bp, _, comparisons_branching = check_cover(branching, range(9))
    assert bc == ((0, 1), (2,), (3,), (4,), (5,))
    assert bp == ((0, 1, 3, 4), (0, 2, 3, 4), (0, 3, 4))
    report = {
        "status": "passed",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "original_exact_optimized_base": 2,
        "refined_exact_optimized_base": 1,
        "row_components": rc, "column_components": cc,
        "rational_path_pair_certificates": report_certificates,
        "integer_dual_certificates": [
            {"total": total, "product": str(dual)},
            {"total": final_total, "product": str(final_dual)},
        ],
        "common_type_sum_one_lengths": list(range(2, 13)),
        "cover_type_comparisons": (comparisons_r + comparisons_c
                                   + comparisons_connected + comparisons_branching),
        "connected_component_negative_control": "passed",
        "branching_and_skip_edge_control_paths": bp,
        "scope": "Abstract counting-graph refinement; no improving crossing or factoring exponent.",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
