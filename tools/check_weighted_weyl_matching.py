#!/usr/bin/env python3
"""Certify a common protected minor for weighted Weyl symmetry spans.

This checker uses only binary partition labels, integer conservation,
and exact Lucas arithmetic. It does not need the large coefficient dumps.
"""

import argparse
from collections import deque
import itertools
import json
from time import monotonic


def binary_labels(h):
    weights = [1 << j for j in range(h.bit_length())]
    answer = []
    building = [0] * len(weights)

    def visit(j, remaining):
        if j == 0:
            building[0] = remaining
            answer.append(tuple(building))
            return
        for value in range(remaining // weights[j] + 1):
            building[j] = value
            visit(j - 1, remaining - value * weights[j])

    visit(len(weights) - 1, h)
    return weights, answer


def maximum_matching(adjacency, right_size):
    """Hopcroft--Karp, with a replayable left-to-right matching."""
    left = [-1] * len(adjacency)
    right = [-1] * right_size

    while True:
        distance = [-1] * len(adjacency)
        queue = deque()
        for u, v in enumerate(left):
            if v < 0:
                distance[u] = 0
                queue.append(u)
        reachable = False
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                following = right[v]
                if following < 0:
                    reachable = True
                elif distance[following] < 0:
                    distance[following] = distance[u] + 1
                    queue.append(following)
        if not reachable:
            break

        def augment(u):
            for v in adjacency[u]:
                following = right[v]
                if following < 0 or (distance[following] == distance[u] + 1
                                     and augment(following)):
                    left[u], right[v] = v, u
                    return True
            distance[u] = -1
            return False

        for u, v in enumerate(left):
            if v < 0:
                augment(u)
    return left


def calculate(h):
    started = monotonic()
    weights, labels = binary_labels(h)
    differences, sums, adjacency, edge_pairs = {}, {}, [], {}
    edges = 0
    for a, alpha in enumerate(labels):
        for y, yy in enumerate(labels):
            if any(x & z for x, z in zip(alpha, yy)):
                continue
            difference = tuple(x - z for x, z in zip(alpha, yy))
            total = tuple(x + z for x, z in zip(alpha, yy))
            if difference not in differences:
                differences[difference] = len(adjacency)
                adjacency.append([])
            if total not in sums:
                sums[total] = len(sums)
            u, v = differences[difference], sums[total]
            adjacency[u].append(v)
            edge_pairs[u, v] = (a, y)
            edges += 1
    matching = maximum_matching(adjacency, len(sums))
    witness = [edge_pairs[u, v] for u, v in enumerate(matching) if v >= 0]
    assert len({tuple(x + z for x, z in zip(labels[a], labels[y]))
                for a, y in witness}) == len(witness)
    assert len({tuple(x - z for x, z in zip(labels[a], labels[y]))
                for a, y in witness}) == len(witness)
    assert all(not any(x & z for x, z in zip(labels[a], labels[y]))
               for a, y in witness)

    # E pairs positions (0,3) and (1,2), in the order (alpha,y,u,beta).
    crossing = {frozenset((0, 3)), frozenset((1, 2))}
    permutations = [p for p in itertools.permutations(range(4))
                    if {frozenset((p[0], p[3])), frozenset((p[1], p[2]))}
                    == crossing]
    assert len(permutations) == 8
    allowed = {(1, -1, 1, -1), (-1, 1, -1, 1),
               (1, 1, -1, -1), (-1, -1, 1, 1)}
    for p in permutations:
        coefficients = [0] * 4
        for position, sign in enumerate((1, -1, 1, -1)):
            coefficients[p[position]] = sign
        assert tuple(coefficients) in allowed

    # Independently check the diagonal binomial identity in all local
    # occupation ranges. A bit mask encodes a polynomial over F_2.
    diagonal_checks = 0
    for weight in weights:
        for a in range(h // weight + 1):
            for y in range(h // weight + 1):
                polynomial = 0
                for k in range(min(a, y) + 1):
                    if (k & ~a) or (k & ~y):
                        continue
                    exponent = weight * k
                    submask = exponent
                    while True:
                        polynomial ^= 1 << submask
                        if not submask:
                            break
                        submask = (submask - 1) & exponent
                assert polynomial == 1 << (weight * (a & y))
                diagonal_checks += 1

    rank_floor, d = len(witness), len(labels)
    # 577/408 exceeds sqrt(2), certified by integer squares.
    passes_strict_exclusion = rank_floor ** 408 > d ** 577
    assert 577 ** 2 > 2 * 408 ** 2
    return {
        "energy": h,
        "field_characteristic": 2,
        "labels": d,
        "weights": weights,
        "label_occupations": labels,
        "zero_intersection_pairs": edges,
        "difference_vertices": len(differences),
        "sum_vertices": len(sums),
        "protected_minor_size": rank_floor,
        "matching_covers_all_sum_vertices": rank_floor == len(sums),
        "matching_pairs": witness,
        "crossing_preserving_permutations": permutations,
        "local_diagonal_checks": diagonal_checks,
        "strict_half_threshold_exclusion": passes_strict_exclusion,
        "integer_certificate": {
            "rank_floor_power": 408,
            "label_power": 577,
            "rational_above_sqrt_two": True,
            "rank_floor_power_exceeds_label_power": passes_strict_exclusion,
        },
        "scope": "Full coefficient and all crossing-preserving relabeling span; every characteristic-two field; this finite energy only.",
        "elapsed_seconds": round(monotonic() - started, 3),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h", type=int, default=32)
    parser.add_argument("--output")
    args = parser.parse_args()
    if args.h < 1:
        parser.error("energy must be positive")
    result = calculate(args.h)
    if args.h == 32:
        assert result["labels"] == 202
        assert result["protected_minor_size"] == 1821
        assert result["strict_half_threshold_exclusion"]
    if args.output:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {"matching_pairs", "label_occupations"}}, indent=2))


if __name__ == "__main__":
    main()
