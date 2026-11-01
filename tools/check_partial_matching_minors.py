#!/usr/bin/env python3
"""Check the component formula against independent exact determinants."""

import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path
from random import Random

from partial_matching_minors import disjointness_determinant


def partial_matchings(vertices):
    if not vertices:
        yield ()
        return
    first, *rest = vertices
    yield from partial_matchings(rest)
    for index, second in enumerate(rest):
        edge = (1 << first) | (1 << second)
        for tail in partial_matchings(rest[:index] + rest[index + 1:]):
            yield (edge,) + tail


def direct_determinant(first, second):
    """Gaussian elimination, with no component or graph information."""
    rows = [sum(int(edge & other == 0) << j
                for j, other in enumerate(second)) for edge in first]
    size = len(rows)
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i] >> column & 1), None)
        if pivot is None:
            return 0
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(column + 1, size):
            if rows[row] >> column & 1:
                rows[row] ^= rows[column]
    return 1


def verify():
    cases = Counter()
    profiles = []
    for vertices in range(9):
        by_size = defaultdict(list)
        for matching in partial_matchings(list(range(vertices))):
            by_size[len(matching)].append(matching)
        checked = 0
        for family in by_size.values():
            for first, second in product(family, repeat=2):
                answer, case = disjointness_determinant(first, second, with_case=True)
                assert answer == direct_determinant(first, second), (first, second, case)
                cases[(case, answer)] += 1
                checked += 1
        profiles.append({"vertices": vertices, "partial_matchings": sum(map(len, by_size.values())),
                         "equal_size_pairs_checked": checked})
    rng = Random(20261030)
    randomized = 0
    for vertices in (10, 12, 16, 24, 40):
        for _ in range(200):
            size = rng.randrange(vertices // 2 + 1)
            pair = []
            for _ in range(2):
                order = list(range(vertices))
                rng.shuffle(order)
                pair.append(tuple((1 << order[2 * i]) | (1 << order[2 * i + 1])
                                  for i in range(size)))
            answer, case = disjointness_determinant(*pair, with_case=True)
            assert answer == direct_determinant(*pair)
            cases[(case, answer)] += 1
            randomized += 1
    assert all(cases[(case, answer)] for case in
               ("balanced_paths_only", "one_cycle", "one_unbalanced_pair") for answer in (0, 1))
    assert cases[("nullity_at_least_two", 0)]
    invalid = 0
    for first, second in [((3,), ()), ((1,), (3,)), ((3, 6), (3, 12)), ((3, 12), (3, 6))]:
        try:
            disjointness_determinant(first, second)
        except ValueError:
            invalid += 1
        else:
            raise AssertionError("Invalid matching accepted")
    return {"scope": "exact partial-matching determinant queries; no crossing rank bound",
            "exhaustive_profiles": profiles,
            "exhaustive_pairs": sum(item["equal_size_pairs_checked"] for item in profiles),
            "deterministic_random_pairs": randomized,
            "component_cases": [{"case": case, "determinant": value, "count": count}
                                for (case, value), count in sorted(cases.items())],
            "invalid_input_controls": invalid}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    report = verify()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(rendered)
    print(rendered, end="")
