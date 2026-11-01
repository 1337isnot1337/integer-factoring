#!/usr/bin/env python3
"""Exact controls for the all-field mixed-cut crossing rank budget."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import isqrt
from pathlib import Path


def rank_binary(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def cuts(mask, q):
    n = q * q
    row = [(mask >> (n * i)) & ((1 << n) - 1) for i in range(n)]
    column = [0] * n
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    if (row[q * x + y] >> (q * u + v)) & 1:
                        column[q * x + u] |= 1 << (q * v + y)
    return rank_binary(row), rank_binary(column)


def crossing(q):
    n = q * q
    return sum(1 << (n * (q * x + y) + q * y + x)
               for x in range(q) for y in range(q))


def identity(q):
    n = q * q
    return sum(1 << (n * i + i) for i in range(n))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    q = 2
    target = crossing(q)
    profiles = [cuts(mask, q) for mask in range(1 << 16)]
    mixed_histogram = Counter()
    cheap_histogram = Counter()
    for mask, (r, c) in enumerate(profiles):
        other_r, other_c = profiles[mask ^ target]
        assert r + other_c >= q
        mixed_histogram[r + other_c] += 1
        cheap = min(r, c) + min(other_r, other_c)
        assert cheap >= q
        cheap_histogram[cheap] += 1

    unit = identity(q)
    negative = [profiles[unit], profiles[unit ^ target]]
    assert negative == [(4, 1), (1, 4)]
    assert sum(min(pair) for pair in negative) == 2 < q * q

    scale = 10 ** 6
    cutoff = []
    for q in range(2, 13):
        n = q * q
        radicand = n - q + 1
        root_lower = isqrt(radicand * scale * scale)
        lower = (q - 2) * scale + 2 * root_lower
        upper = lower + 2
        assert root_lower ** 2 <= radicand * scale ** 2 < (root_lower + 1) ** 2
        # Every comparison below is between integers, with positive denominators.
        excludes = lower ** 70 > q ** 99 * scale ** 70
        if q <= 11:
            assert excludes
        else:
            assert upper ** 99 < q ** 140 * scale ** 99
            assert not excludes
        abstract = [(radicand, 1), (1, radicand)] + [(1, 1)] * (q - 2)
        assert sum(r for r, _ in abstract) == n
        assert sum(c for _, c in abstract) == n
        assert sum(min(r, c) for r, c in abstract) == q
        cutoff.append({'labels': q, 'radicand': radicand,
                       'lower_numerator': lower, 'upper_numerator': upper,
                       'denominator': scale,
                       'excluded_by_exact_power_comparison': excludes})
    assert 140 ** 2 < 2 * 99 ** 2
    assert 99 ** 2 > 2 * 70 ** 2

    evidence = {
        'scope': 'All-field theorem; finite binary controls and exact rational cutoff certificates.',
        'binary_two_label_tensors': len(profiles),
        'mixed_budget_histogram': dict(sorted(mixed_histogram.items())),
        'sum_of_cheaper_cuts_histogram': dict(sorted(cheap_histogram.items())),
        'negative_control_against_q_squared_budget': negative,
        'cutoff_certificates': cutoff,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = json.dumps(evidence, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
