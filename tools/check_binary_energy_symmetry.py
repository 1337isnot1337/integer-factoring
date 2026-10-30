#!/usr/bin/env python3
"""Exact embedding checks for the all-field binary-energy shell barrier."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb, prod
from pathlib import Path
import random


def parameters(r):
    h = 1 << r
    ceil_log = (r + 1).bit_length()
    cutoff = r - 2 * ceil_log - 2
    free = [4 * (s // 4) for s in range(max(0, cutoff + 1))]
    energy = sum((1 << s) * q // 2 for s, q in enumerate(free))
    return h, cutoff, free, energy


def shell_count(r):
    h = 1 << r
    counts = [0] * (h + 1)
    counts[0] = 1
    for s in range(r + 1):
        for _ in range(s + 1):
            for target in range(h, (1 << s) - 1, -1):
                counts[target] += counts[target - (1 << s)]
    return counts[h]


def partition_count(r):
    h = 1 << r
    counts = [0] * (h + 1)
    counts[0] = 1
    for s in range(r + 1):
        for target in range(1 << s, h + 1):
            counts[target] += counts[target - (1 << s)]
    return counts[h]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    exact_counts = []
    for r in range(1, 13):
        count = shell_count(r)
        assert count == partition_count(r)
        h, cutoff, free, energy = parameters(r)
        assert prod(comb(q, q // 2) for q in free) <= count
        assert count <= 1 << ((r + 1) * (r + 2) // 2)
        exact_counts.append({"R": r, "shell_size": count})

    envelopes = 0
    for r in range(1, 1025):
        h, cutoff, free, energy = parameters(r)
        assert 0 <= energy < h
        assert h - energy <= h
        assert all(q % 4 == 0 and q <= s for s, q in enumerate(free))
        total = sum(free)
        if cutoff >= 4:
            assert energy < cutoff * (1 << cutoff)
            assert (1 << (r - cutoff)) >= 4 * (r + 2)**2
            assert total >= cutoff * (cutoff + 1) // 2 - 3 * (cutoff + 1)
            assert total <= cutoff * (cutoff + 1) // 2
        for q in free:
            assert comb(q, q // 2) * (q + 1) >= 1 << q
        envelopes += 1

    rng = random.Random(20261030)
    embeddings = 0
    for r in (14, 16, 20, 24, 32, 48, 64):
        h, cutoff, free, energy = parameters(r)
        filler = h - energy
        seen = {}
        for _ in range(128):
            coordinates = []
            choices = []
            for s in range(r + 1):
                q = free[s] if s < len(free) else 0
                selected = frozenset(rng.sample(range(q), q // 2))
                choices.append(tuple(sorted(selected)))
                bits = [(filler >> s) & 1]
                bits += [int(i in selected) for i in range(q)]
                bits += [0] * (s - q)
                assert len(bits) == s + 1
                coordinates.append(tuple(bits))
            assert sum((1 << s) * sum(bits) for s, bits in enumerate(coordinates)) == h
            encoded = tuple(coordinates)
            original = tuple(choices)
            assert encoded not in seen or seen[encoded] == original
            seen[encoded] = original
            embeddings += 1

    result = {
        "status": "passed",
        "scope": "Exact shell counts, integer envelopes, and reserved-bit embeddings; the all-field rank theorem is proved separately.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "independent_shell_counts": exact_counts,
        "integer_resource_envelopes": envelopes,
        "reserved_bit_embedding_checks": embeddings,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
