#!/usr/bin/env python3
"""Check integer stationary-flow lower certificates for channel bounds.

The lower bound concerns the reciprocal-weight positive-supervector
method. It is not a lower bound on the true common-type counts or ranks.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
import random


@dataclass(frozen=True)
class Graph:
    order: int
    pieces: int
    edges: tuple[tuple[int, int, int], ...]


def validate_flow(graph, counts):
    if any(not isinstance(e, int) or not 0 <= e < len(graph.edges)
           for e in counts):
        raise ValueError("flow uses an edge outside the graph")
    if any(not isinstance(a, int) or a < 0 for a in counts.values()):
        raise ValueError("flow counts must be nonnegative integers")
    outgoing = [0] * graph.order
    incoming = [0] * graph.order
    piece_totals = [0] * graph.pieces
    for e, (s, t, label) in enumerate(graph.edges):
        a = counts.get(e, 0)
        outgoing[s] += a
        incoming[t] += a
        piece_totals[label] += a
    if outgoing != incoming:
        raise ValueError("flow is not stationary")
    total = sum(outgoing)
    if not total:
        raise ValueError("flow has zero total")
    numerator = prod(a**a for a in outgoing if a)
    denominator = prod(a**a for a in counts.values() if a)
    return total, tuple(piece_totals), Fraction(numerator, denominator)


def coupled_certificate(row_graph, row_flow, col_graph, col_flow):
    rd, rt, rk = validate_flow(row_graph, row_flow)
    cd, ct, ck = validate_flow(col_graph, col_flow)
    if rd != cd or rt != ct:
        raise ValueError("row and column piece totals do not match")
    return rd, rk * ck


def upper_beta(graph, weights, vector):
    assert all(z > 0 for z in weights)
    assert all(h > 0 for h in vector)
    image = [Fraction(0) for _ in range(graph.order)]
    for s, t, label in graph.edges:
        image[s] += weights[label] * vector[t]
    return max(a / h for a, h in zip(image, vector))


def complete_graph(piece_labels):
    n = len(piece_labels)
    return Graph(n, max(piece_labels) + 1,
                 tuple((s, t, piece_labels[t])
                       for s in range(n) for t in range(n)))


def product_flow(graph, probabilities, denominator):
    counts = {}
    for e, (s, t, _) in enumerate(graph.edges):
        value = denominator * probabilities[s] * probabilities[t]
        assert value.denominator == 1
        counts[e] = value.numerator
    return counts


def walk_types(graph, length):
    states = Counter((s, (0,) * graph.pieces) for s in range(graph.order))
    for _ in range(length):
        following = Counter()
        for (s, counts), multiplicity in states.items():
            for source, target, label in graph.edges:
                if source == s:
                    changed = list(counts)
                    changed[label] += 1
                    following[target, tuple(changed)] += multiplicity
        states = following
    result = Counter()
    for (_, counts), multiplicity in states.items():
        result[counts] += multiplicity
    return result


def expect_rejection(function):
    try:
        function()
    except ValueError:
        return
    raise AssertionError("invalid flow certificate was accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    labels = tuple(i for i in range(3) for _ in range(4))
    nilpotent = Graph(12, 3, tuple((s, t, labels[t])
                                 for s in range(12) for t in range(12)
                                 if labels[s] != labels[t]))
    unit_flow = {e: 1 for e in range(len(nilpotent.edges))}
    d8, bound8 = coupled_certificate(nilpotent, unit_flow, nilpotent, unit_flow)
    assert d8 == 96 and bound8 == 8**192
    ones = [Fraction(1)] * 12
    assert upper_beta(nilpotent, [Fraction(1)] * 3, ones) == 8

    row = complete_graph((0, 1, 1, 1, 1))
    col = complete_graph((0, 0, 0, 0, 1))
    row_prob = [Fraction(1, 2)] + [Fraction(1, 8)] * 4
    col_prob = [Fraction(1, 8)] * 4 + [Fraction(1, 2)]
    row_flow = product_flow(row, row_prob, 64)
    col_flow = product_flow(col, col_prob, 64)
    d4, bound4 = coupled_certificate(row, row_flow, col, col_flow)
    assert d4 == 64 and bound4 == 4**128
    assert upper_beta(row, [Fraction(2), Fraction(1, 2)], [Fraction(1)] * 5) == 4
    assert upper_beta(col, [Fraction(1, 2), Fraction(2)], [Fraction(1)] * 5) == 4

    rng = random.Random(20261030)
    random_checks = 0
    for rg, rf, cg, cf in ((nilpotent, unit_flow, nilpotent, unit_flow),
                          (row, row_flow, col, col_flow)):
        total, lower = coupled_certificate(rg, rf, cg, cf)
        for _ in range(48):
            weights = [Fraction(rng.randrange(1, 9), rng.randrange(1, 9))
                       for _ in range(rg.pieces)]
            rh = [Fraction(rng.randrange(1, 9), rng.randrange(1, 9))
                  for _ in range(rg.order)]
            ch = [Fraction(rng.randrange(1, 9), rng.randrange(1, 9))
                  for _ in range(cg.order)]
            rb = upper_beta(rg, weights, rh)
            cb = upper_beta(cg, [1 / z for z in weights], ch)
            assert (rb * cb)**total >= lower
            random_checks += 1

    invalid = dict(unit_flow)
    invalid[0] = 0
    expect_rejection(lambda: validate_flow(nilpotent, invalid))
    expect_rejection(lambda: validate_flow(nilpotent, {len(nilpotent.edges): 1}))
    expect_rejection(lambda: validate_flow(nilpotent, {}))
    loops = Graph(2, 2, ((0, 0, 0), (1, 1, 1)))
    cycle = Graph(2, 2, ((0, 1, 0), (1, 0, 1)))
    expect_rejection(lambda: coupled_certificate(loops, {0: 2}, cycle, {0: 1, 1: 1}))

    dc, kc = coupled_certificate(loops, {0: 1, 1: 1}, cycle, {0: 1, 1: 1})
    assert dc == 2 and kc == 1
    assert upper_beta(loops, [Fraction(1)] * 2, [Fraction(1)] * 2) == 1
    assert upper_beta(cycle, [Fraction(1)] * 2, [Fraction(1)] * 2) == 1
    for length in range(2, 13):
        row_types = walk_types(loops, length)
        col_types = walk_types(cycle, length)
        assert not (row_types.keys() & col_types.keys())

    report = {
        "status": "passed",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "exact_optimized_certificate_bases": [8, 4],
        "square_zero_flow": {"vertices": 12, "edges": 96, "total": d8,
                             "piece_totals": [32, 32, 32],
                             "dual_product": "8^192"},
        "complete_graph_flow": {"total": d4, "piece_totals": [32, 32],
                                "dual_product": "4^128"},
        "independent_rational_supervector_comparisons": random_checks,
        "invalid_certificate_rejections": 4,
        "spectral_vs_common_type_negative_control_lengths": list(range(2, 13)),
        "scope": "Optimality of the weighted-supervector bound only; no tensor-rank or factoring lower bound.",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
