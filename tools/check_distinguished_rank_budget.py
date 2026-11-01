#!/usr/bin/env python3
"""Exact verification of the distinguished-piece rank-budget relaxation."""

import argparse
from collections import Counter
from functools import lru_cache
from math import isqrt
import json
from pathlib import Path


@lru_cache(None)
def root_form(value):
    """Return a,s with sqrt(value)=a*sqrt(s), s square-free."""
    coefficient = 1
    prime = 2
    while prime * prime <= value:
        while value % (prime * prime) == 0:
            value //= prime * prime
            coefficient *= prime
        prime += 1
    return coefficient, value


def surds(values):
    result = Counter()
    for value in values:
        coefficient, square_free = root_form(value)
        result[square_free] += coefficient
    return tuple(sorted(result.items()))


def enclosure(value, bits):
    scale = 1 << bits
    low = high = 0
    for radicand, coefficient in value:
        numerator = radicand * scale * scale
        floor = isqrt(numerator)
        ceiling = floor if floor * floor == numerator else floor + 1
        low += coefficient * floor
        high += coefficient * ceiling
    return low, high


def compare(left, right):
    if left == right:
        return 0
    for bits in (8, 16, 32, 64, 128):
        lo_left, hi_left = enclosure(left, bits)
        lo_right, hi_right = enclosure(right, bits)
        if hi_left < lo_right:
            return -1
        if hi_right < lo_left:
            return 1
    raise AssertionError("square-root intervals did not separate")


@lru_cache(None)
def compositions(total, count):
    if count == 1:
        return ((total,),) if total >= 1 else ()
    return tuple((first,) + tail
                 for first in range(1, total - count + 2)
                 for tail in compositions(total - first, count - 1))


def formula(n, m):
    if m >= n - 1:
        return surds((n * n,))
    candidates = (surds((n * n,)),
                  surds((m * (n - 1), n - m)),
                  surds((m * m, n - m - 1, n - m - 1)))
    best = candidates[0]
    for candidate in candidates[1:]:
        if compare(candidate, best) < 0:
            best = candidate
    return best


def verify():
    parameters = vector_pairs = 0
    for n in range(2, 10):
        for m in range(1, n + 1):
            best = None
            for count in range(1, n - m + 2):
                vectors = tuple(vector for vector in compositions(n, count) if vector[0] >= m)
                for row in vectors:
                    for column in vectors:
                        value = surds(tuple(a * b for a, b in zip(row, column)))
                        if best is None or compare(value, best) < 0:
                            best = value
                        vector_pairs += 1
            assert best == formula(n, m), (n, m, best, formula(n, m))
            parameters += 1

    # The three displayed vertex types realize the three candidates.
    attainment_cases = 0
    for n in range(3, 101):
        for m in range(1, n - 1):
            vectors = (((n,), (n,)),
                       ((n - 1, 1), (m, n - m)),
                       ((m, n - m - 1, 1), (m, 1, n - m - 1)))
            attained = []
            for row, column in vectors:
                assert sum(row) == sum(column) == n
                assert min(row) >= 1 and min(column) >= 1
                assert row[0] >= m and column[0] >= m
                attained.append(surds(tuple(a * b for a, b in zip(row, column))))
            assert any(value == formula(n, m) for value in attained)
            attainment_cases += 1

    special = surds((72 * 72, 503, 503))
    assert special == formula(576, 72)
    assert 22 * 22 < 503 < 23 * 23
    assert 118 < 203 and 203 * 203 < 72 * 575
    assert 116 ** 7 > 24 ** 10
    assert formula(36, 6) == surds((36, 29, 29))
    return {
        "status": "pass",
        "exhaustive_parameters": parameters,
        "integer_rank_vector_pairs": vector_pairs,
        "attainment_parameters": attainment_cases,
        "permutation_b4_lower_bound": "72 + 2*sqrt(503)",
        "permutation_b4_lower_bound_strictly_exceeds": 116,
        "scope": "Sharp numerical rank-budget relaxation, not tensor-rank realization.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
