#!/usr/bin/env python3
"""Exact convolution bounds and field-free envelopes for Pascal complements.

The checker uses integer comparisons and modular elimination only. It does
not assert a universal entropy floor or an improved factoring exponent.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from check_schur_pascal_complement import greedy_complement


def complement_envelope(indices):
    """Largest disjoint complement satisfying all translated-basis deadlines."""
    indices = tuple(sorted(indices))
    n = len(indices)
    if not n or len(set(indices)) != n or indices[0] < 0:
        raise ValueError("a nonempty set of nonnegative integers is required")
    largest = 2 * indices[-1] + 1
    deadlines = [largest] * n
    for i, a in enumerate(indices, 1):
        for j, b in enumerate(indices, 1):
            bound = a + b + 1
            required = i + j - bisect_right(indices, bound)
            if required > 0:
                deadlines[required - 1] = min(deadlines[required - 1], bound)
    for k in range(n - 2, -1, -1):
        deadlines[k] = min(deadlines[k], deadlines[k + 1])
    excluded = set(indices)
    result = [0] * n
    following = largest + 1
    for k in range(n - 1, -1, -1):
        candidate = min(deadlines[k], following - 1)
        while candidate in excluded:
            candidate -= 1
        if candidate < 0:
            raise ArithmeticError("the proven complement envelope should exist")
        result[k] = candidate
        following = candidate
    return tuple(result), tuple(deadlines)


def modular_complement(indices, prime):
    """Independent direct Pascal-column scan with modular Gaussian elimination."""
    indices = tuple(indices)
    excluded = set(indices)
    basis = {}
    result = []
    for c in range(2 * max(indices) + 2):
        if c in excluded:
            continue
        column = [comb(c, a) % prime if a <= c else 0 for a in indices]
        for pivot in sorted(basis):
            scale = column[pivot]
            if scale:
                column = [(a - scale * b) % prime for a, b in zip(column, basis[pivot])]
        pivot = next((i for i, value in enumerate(column) if value), None)
        if pivot is not None:
            inverse = pow(column[pivot], -1, prime)
            basis[pivot] = [value * inverse % prime for value in column]
            result.append(c)
        if len(result) == len(indices):
            return tuple(result)
    raise AssertionError("translated consecutive columns must span")


def check_bounds(indices, complement):
    indices, complement = tuple(indices), tuple(complement)
    n = len(indices)
    assert len(complement) == n
    assert len(set(indices + complement)) == 2 * n
    union = tuple(sorted(indices + complement))
    for i, a in enumerate(indices, 1):
        assert union[i - 1] <= a
        assert union[2 * i - 2] <= 2 * a
        assert union[2 * i - 1] <= 2 * a + 1
        for j, b in enumerate(indices, 1):
            assert union[i + j - 1] <= a + b + 1
    envelope, deadlines = complement_envelope(indices)
    assert all(c <= d <= h for c, d, h in zip(complement, envelope, deadlines))
    assert len(set(indices + envelope)) == 2 * n
    assert all(envelope[k] < envelope[k + 1] for k in range(n - 1))
    # Verify the equivalent prefix constraint independently on small ranges.
    prefix_checks = 0
    if max(indices) <= 64:
        for x in range(2 * max(indices) + 2):
            actual = bisect_right(union, x)
            required = max(bisect_right(indices, a) + bisect_right(indices, x - a - 1)
                           for a in range(-1, x + 1))
            assert actual >= required
            prefix_checks += 1
    return envelope, n * n, prefix_checks


def check_consecutive_factorization():
    entries = 0
    for degree in range(13):
        for start in range(18):
            for a in range(degree + 1):
                for j in range(degree + 1):
                    left = comb(start + j, a) if a <= start + j else 0
                    right = sum((comb(start, a - h) if a - h <= start else 0)
                                * (comb(j, h) if h <= j else 0) for h in range(a + 1))
                    assert left == right
                    entries += 1
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--ambient", type=int, default=12,
                        help="enumerate every nonempty subset of this small interval")
    args = parser.parse_args()
    if not 1 <= args.ambient <= 18:
        raise ValueError("finite exhaustive ambient must lie between 1 and 18")
    pair_checks = prefix_checks = 0
    differing_envelopes = 0
    for mask in range(1, 1 << args.ambient):
        indices = tuple(i for i in range(args.ambient) if mask >> i & 1)
        chosen = greedy_complement(indices)
        envelope, pairs, prefixes = check_bounds(indices, chosen)
        pair_checks += pairs
        prefix_checks += prefixes
        differing_envelopes += int(envelope != chosen)
        if args.ambient <= 10 or mask < 256:
            assert chosen == modular_complement(indices, 2)
    other_fields = []
    for prime in (3, 5, 7):
        count = 0
        for mask in range(1, 1 << 9):
            indices = tuple(i for i in range(9) if mask >> i & 1)
            chosen = modular_complement(indices, prime)
            _, pairs, prefixes = check_bounds(indices, chosen)
            pair_checks += pairs
            prefix_checks += prefixes
            count += 1
        other_fields.append({"characteristic": prime, "mode_sets": count})
    intervals = []
    for n in range(1, 65):
        indices = tuple(range(n))
        chosen = greedy_complement(indices)
        envelope, _ = complement_envelope(indices)
        expected = tuple(range(n, 2 * n))
        assert chosen == envelope == expected
        union = tuple(range(2 * n))
        assert all(union[i + j - 1] == a + b + 1
                   for i, a in enumerate(indices, 1) for j, b in enumerate(indices, 1))
        intervals.append(n)
    large = []
    for bits in (20, 64, 128, 256):
        power = 1 << bits
        for indices in ((power - 1,), (0, power), (power + 1, power + 3, power + 7),
                        (0, power - 1, power, 2 * power - 1)):
            chosen = greedy_complement(indices)
            envelope, pairs, prefixes = check_bounds(indices, chosen)
            pair_checks += pairs
            prefix_checks += prefixes
            large.append({"index_bits": bits, "indices": indices,
                          "complement": chosen, "envelope": envelope})
    sample = (0, 1, 2, 4)
    chosen = greedy_complement(sample)
    envelope, deadlines = complement_envelope(sample)
    assert chosen == (3, 5, 6, 7) and envelope == (3, 5, 7, 9)
    report = {"status": "exact finite checks of all-size convolution bounds; no entropy-floor claim",
              "integer_factorization_entries": check_consecutive_factorization(),
              "binary_exhaustive_ambient": args.ambient,
              "binary_mode_sets": (1 << args.ambient) - 1,
              "binary_sets_with_strict_envelope": differing_envelopes,
              "pair_inequalities_checked": pair_checks,
              "prefix_inequalities_checked": prefix_checks,
              "other_fields": other_fields, "sharp_interval_sizes": intervals,
              "large_index_cases": large,
              "strict_example": {"indices": sample, "complement": chosen,
                                 "envelope": envelope, "deadlines": deadlines},
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
