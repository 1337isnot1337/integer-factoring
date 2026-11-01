#!/usr/bin/env python3
"""Exact finite checks for connected matching-type flow growth.

The general limsup theorem is proved in the companion note. These checks
use integer walk counts and rational entropy ratios, without logarithms.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import factorial, prod
from pathlib import Path

from check_channel_flow_duals import Graph, coupled_certificate, validate_flow


def totals(graph, counts):
    outgoing = [0] * graph.order
    incoming = [0] * graph.order
    kind = [0] * graph.pieces
    for (s, t, label), a in zip(graph.edges, counts):
        outgoing[s] += a
        incoming[t] += a
        kind[label] += a
    return outgoing, incoming, tuple(kind)


def validate_base(graph, counts, start, finish):
    if len(counts) != len(graph.edges) or not counts:
        raise ValueError("base count length or empty support")
    if any(not isinstance(a, int) or a <= 0 for a in counts):
        raise ValueError("base counts must be positive on the whole support")
    outgoing, incoming, kind = totals(graph, counts)
    expected = [int(v == start) - int(v == finish) for v in range(graph.order)]
    if [a - b for a, b in zip(outgoing, incoming)] != expected:
        raise ValueError("base does not have the stated Euler balance")
    used = {v for edge in graph.edges for v in edge[:2]}
    if start not in used or finish not in used:
        raise ValueError("endpoint outside support")
    found = {start}
    while True:
        following = found | {t for s, t, _ in graph.edges if s in found}
        following |= {s for s, t, _ in graph.edges if t in found}
        if following == found:
            break
        found = following
    if found != used:
        raise ValueError("support is not weakly connected")
    return kind


def validate_base_pair(row, ra, rs, rt, col, ca, cs, ct):
    first = validate_base(row, ra, rs, rt)
    second = validate_base(col, ca, cs, ct)
    if first != second:
        raise ValueError("base piece totals do not agree")
    return sum(ra)


def edge_count_words(graph, start, finish, length):
    states = Counter({(start, (0,) * len(graph.edges)): 1})
    for _ in range(length):
        following = Counter()
        for (vertex, counts), multiplicity in states.items():
            for e, (s, t, _) in enumerate(graph.edges):
                if s != vertex:
                    continue
                changed = list(counts)
                changed[e] += 1
                following[t, tuple(changed)] += multiplicity
        states = following
    return Counter({counts: multiplicity for (vertex, counts), multiplicity
                    in states.items() if vertex == finish})


def type_counts(graph, start, finish, length):
    result = Counter()
    for counts, multiplicity in edge_count_words(graph, start, finish, length).items():
        result[totals(graph, counts)[2]] += multiplicity
    return result


def entropy_ratio(graph, counts):
    outgoing, _, _ = totals(graph, counts)
    return Fraction(prod(a**a for a in outgoing if a),
                    prod(a**a for a in counts if a))


def check_count_classes(graph, lengths):
    checks = 0
    for length in lengths:
        for start in range(graph.order):
            for finish in range(graph.order):
                for counts, words in edge_count_words(graph, start, finish, length).items():
                    support = [(edge, a) for edge, a in zip(graph.edges, counts) if a]
                    restricted = Graph(graph.order, graph.pieces, tuple(e for e, _ in support))
                    positive = tuple(a for _, a in support)
                    validate_base(restricted, positive, start, finish)
                    outgoing, _, _ = totals(restricted, positive)
                    vertices = len({v for e in restricted.edges for v in e[:2]})
                    numerator = prod(factorial(a) for a in outgoing)
                    denominator = prod(factorial(a) for a in positive)
                    assert numerator % denominator == 0
                    queues = numerator // denominator
                    ratio = entropy_ratio(restricted, positive)
                    assert words <= queues <= ratio
                    assert words * (length + 1)**vertices >= queues
                    assert queues * (length + 1)**(2 * len(positive)) >= ratio
                    checks += 1
    return checks


def expect_rejection(function):
    try:
        function()
    except ValueError:
        return
    raise AssertionError("invalid base certificate accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    parallel_row = Graph(1, 1, ((0, 0, 0),) * 2)
    parallel_col = Graph(1, 1, ((0, 0, 0),) * 3)
    assert validate_base_pair(parallel_row, (3, 3), 0, 0,
                              parallel_col, (2, 2, 2), 0, 0) == 6
    d6, dual6 = coupled_certificate(parallel_row, {0: 3, 1: 3},
                                    parallel_col, {0: 2, 1: 2, 2: 2})
    assert d6 == 6 and dual6 == 6**6
    for k in range(1, 13):
        rr = type_counts(parallel_row, 0, 0, k)
        cc = type_counts(parallel_col, 0, 0, k)
        assert rr == {(k,): 2**k} and cc == {(k,): 3**k}
        assert rr[k,] * cc[k,] == 6**k

    # A positive bridge is used once. The stationary direction is zero
    # on it, and has half of each piece type in both graphs.
    bridge = Graph(2, 2, ((0, 0, 0), (0, 0, 0),
                          (1, 1, 1), (1, 1, 1), (0, 1, 1)))
    cycle = Graph(2, 2, ((0, 1, 0), (0, 1, 0),
                         (1, 0, 1), (1, 0, 1)))
    ra, ca = (1, 1, 1, 1, 1), (1, 1, 1, 2)
    assert validate_base_pair(bridge, ra, 0, 1, cycle, ca, 1, 0) == 5
    row_flow, col_flow = (1, 1, 1, 1, 0), (1, 1, 1, 1)
    d2, dual2 = coupled_certificate(bridge, dict(enumerate(row_flow)),
                                    cycle, dict(enumerate(col_flow)))
    assert d2 == 4 and dual2 == 2**8
    repeated_base_checks = 0
    for n in range(7):
        rb = tuple(a + n * u for a, u in zip(ra, row_flow))
        cb = tuple(a + n * u for a, u in zip(ca, col_flow))
        assert validate_base_pair(bridge, rb, 0, 1, cycle, cb, 1, 0) == 5 + 4*n
        repeated_base_checks += 1
    # At each odd length 2m+1, the cycle has 2**k words of its one
    # type. The matching bridge walks have the unique bridge position
    # after m A-loops, with 2**(k-1) choices for all remaining loops.
    for k in range(1, 12):
        rr, cc = type_counts(bridge, 0, 1, k), type_counts(cycle, 1, 0, k)
        if k % 2 == 0:
            assert not cc
        else:
            kind = (k // 2, k // 2 + 1)
            assert cc == {kind: 2**k}
            assert rr[kind] == 2**(k-1)

    one_cycle = Graph(2, 1, ((0, 1, 0), (1, 0, 0)))
    parity_rows = []
    for k in range(1, 15):
        closed = type_counts(one_cycle, 0, 0, k)
        opened = type_counts(one_cycle, 0, 1, k)
        assert bool(closed) == (k % 2 == 0)
        assert bool(opened) == (k % 2 == 1)
        assert not closed.keys() & opened.keys()
        parity_rows.append({"length": k, "closed_count": sum(closed.values()),
                            "open_count": sum(opened.values())})
    # The exact balance equations, not this finite list, prove parity:
    # a0=a1 for closed walks, b0=b1+1 for open walks.
    for a1 in range(30):
        for b1 in range(30):
            assert (2*a1) % 2 != (2*b1+1) % 2

    disconnected = Graph(2, 1, ((0, 0, 0), (1, 1, 0)))
    # Balance and a stationary entropy ratio exist, yet no one walk
    # can use both loops. Connectivity must be checked independently.
    validate_flow(disconnected, {0: 1, 1: 1})
    expect_rejection(lambda: validate_base(disconnected, (1, 1), 0, 0))
    expect_rejection(lambda: validate_base(one_cycle, (1, 1), 0, 1))
    expect_rejection(lambda: validate_base(bridge, (1, 1, 1, 1, 0), 0, 1))
    expect_rejection(lambda: validate_base_pair(parallel_row, (1, 1), 0, 0,
                                               parallel_col, (1, 1, 1), 0, 0))
    branching = Graph(3, 2, ((0, 1, 0), (1, 0, 1), (1, 2, 0),
                             (2, 1, 1), (1, 1, 0)))
    classes = sum(check_count_classes(g, range(1, 9)) for g in
                  (parallel_row, parallel_col, bridge, cycle, one_cycle, branching))
    report = {
        "status": "passed",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "exact_edge_count_classes": classes,
        "queue_and_entropy_inequalities_per_class": 4,
        "connected_integer_base_repetitions": repeated_base_checks,
        "integer_growth_certificates": [
            {"base": "sqrt(6)", "flow_total": d6, "entropy_product": str(dual6)},
            {"base": "2", "flow_total": d2, "entropy_product": str(dual2),
             "positive_base_bridge_has_zero_stationary_flow": True},
        ],
        "parity_controls": parity_rows,
        "invalid_base_rejections": 4,
        "scope": "Finite exact counting checks for a proved limsup theorem; no tensor-rank lower bound or factoring improvement.",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
