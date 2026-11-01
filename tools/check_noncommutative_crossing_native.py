#!/usr/bin/env python3
"""Independent Python audit of the native central-correction rank screen.

The native program must already be compiled.  This checker uses the direct
four-coordinate matrix formula and Python integer elimination.  The default
audits bounded blocks; --full replays all 32768 normalized F4 quadratics.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time

from check_noncommutative_crossing_orbit import Heisenberg, binary_rank


MASK64 = (1 << 64) - 1
QUADRATIC_MINIMIZERS_F4 = (
    9553, 9556, 9562, 9567, 24480, 24485, 24491, 24494,
    30001, 30004, 30010, 30015,
)


def quadratic_terms(k: int) -> list[tuple[int, int]]:
    return [(-1, -1)] + [(j, -1) for j in range(2 * k)] + list(
        itertools.combinations(range(2 * k), 2))


def family_terms(k: int, family: str) -> list[tuple[int, int]]:
    if family == "quadratic":
        return quadratic_terms(k)
    if family == "bilinear":
        return [(i, k + j) for i in range(k) for j in range(k)]
    if family == "lookup":
        return list(itertools.product(range(1 << k), repeat=2))
    raise ValueError(family)


def normalized_free_bits(k: int, family: str = "quadratic") -> list[tuple[int, int]]:
    result = []
    for index, (a, b) in enumerate(family_terms(k, family)):
        for bit in range(k):
            if family == "lookup":
                if (a, b) in ((0, 0), (1, 0), (0, 1)) or (a, b, bit) == (1, 1, 0):
                    continue
            else:
                if family == "quadratic" and ((a, b) == (-1, -1) or b == -1 and a in (0, k)):
                    continue
                if (a, b, bit) == (0, k, 0):
                    continue
            result.append((index, bit))
    return result


def correction_from_index(k: int, index: int, family: str = "quadratic") -> tuple[int, ...]:
    q = 1 << k
    terms = family_terms(k, family)
    coefficients = [0] * len(terms)
    for place, (term, bit) in enumerate(normalized_free_bits(k, family)):
        coefficients[term] ^= ((index >> place) & 1) << bit
    if family == "lookup":
        return tuple(coefficients)
    table = []
    for a in range(q):
        for b in range(q):
            word = a | (b << k)
            value = 0
            for (i, j), coefficient in zip(terms, coefficients):
                if (i == -1 or (word >> i) & 1) and (j == -1 or (word >> j) & 1):
                    value ^= coefficient
            table.append(value)
    assert table[0] == table[q] == table[1] == 0
    return tuple(table)


def correction_matrix_blocks(group: Heisenberg) -> list[list[list[tuple[int, int]]]]:
    """Precompute direct entries for each value f(a+A,b+B), independently."""
    q = group.q
    blocks = [[[] for _ in range(q)] for _ in range(q * q)]
    for a, b, A, B in itertools.product(range(q), repeat=4):
        cross = group.field_mul(a, B) ^ group.field_mul(A, b)
        square_sum = group.field_mul(a, b) ^ group.field_mul(A, B)
        for value in range(q):
            row = (a * q + b) * q + (cross ^ value)
            col = (A * q + B) * q + (square_sum ^ value)
            blocks[(a ^ A) * q + (b ^ B)][value].append((row, 1 << col))
    return blocks


def lookup_index(k: int, correction: tuple[int, ...]) -> int:
    result = sum(((correction[term] >> bit) & 1) << place
                 for place, (term, bit) in enumerate(normalized_free_bits(k, "lookup")))
    assert correction_from_index(k, result, "lookup") == correction
    return result


def rank_from_blocks(group: Heisenberg, blocks, correction: tuple[int, ...]) -> int:
    matrix = [0] * group.d
    for s, value in enumerate(correction):
        for row, mask in blocks[s][value]:
            matrix[row] ^= mask
    return binary_rank(matrix)


def audit_block(binary: Path, k: int, start: int, count: int,
                family: str = "quadratic") -> dict:
    began = time.monotonic()
    command = [str(binary), "--k", str(k), "--family", family,
               "--normalized", "--start", str(start), "--count", str(count), "--exact"]
    native = [json.loads(line) for line in subprocess.check_output(command, text=True).splitlines()]
    assert native[0]["free_bits"] == len(normalized_free_bits(k, family))
    group = Heisenberg(k)
    blocks = correction_matrix_blocks(group)
    histogram = Counter()
    digest = 14695981039346656037
    best = group.q + (group.q - 1) * (3 ** k)
    minimizers = []
    for index in range(start, start + count):
        correction = correction_from_index(k, index, family)
        rank = rank_from_blocks(group, blocks, correction)
        histogram[rank] += 1
        if rank < best:
            best = rank
            minimizers = [index]
        elif rank == best:
            minimizers.append(index)
        digest = ((digest ^ index) * 1099511628211) & MASK64
        digest = ((digest ^ rank) * 1099511628211) & MASK64
    result = native[-1]
    assert result["event"] == "complete"
    assert result["histogram"] == [[rank, number] for rank, number in sorted(histogram.items())]
    assert result["digest_fnv64"] == str(digest)
    assert result["best_rank"] == best
    return {"field_degree": k, "family": family, "start": start, "count": count,
            "best_rank": best, "minimizing_indices": minimizers,
            "digest_fnv64": str(digest),
            "histogram": result["histogram"],
            "elapsed_seconds": round(time.monotonic() - began, 4)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--lookup-summary", type=Path,
                        help="also validate the recorded exhaustive lookup summary")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert len(normalized_free_bits(2)) == 15
    assert len(normalized_free_bits(3)) == 56
    assert len(normalized_free_bits(2, "lookup")) == 25
    assert len(normalized_free_bits(3, "bilinear")) == 26
    blocks = [audit_block(args.binary, 2, 0, 32768 if args.full else 256),
              audit_block(args.binary, 3, 0, 32),
              audit_block(args.binary, 3, 9553, 32),
              audit_block(args.binary, 3, (1 << 55) + 479, 32),
              audit_block(args.binary, 2, 0, 256, "lookup"),
              audit_block(args.binary, 2, (1 << 20) + 937, 64, "lookup"),
              audit_block(args.binary, 2, (1 << 25) - 64, 64, "lookup"),
              audit_block(args.binary, 3, 0, 32, "bilinear"),
              audit_block(args.binary, 3, (1 << 25) + 937, 32, "bilinear"),
              audit_block(args.binary, 3, (1 << 26) - 32, 32, "bilinear")]
    if args.full:
        assert blocks[0]["minimizing_indices"] == list(QUADRATIC_MINIMIZERS_F4)
    lookup_winners = []
    for quadratic_index in QUADRATIC_MINIMIZERS_F4:
        table = correction_from_index(2, quadratic_index)
        index = lookup_index(2, table)
        checked = audit_block(args.binary, 2, index, 1, "lookup")
        assert checked["best_rank"] == 28
        lookup_winners.append({"quadratic_index": quadratic_index,
                               "lookup_index": index, "correction": list(table)})
    lookup_summary = None
    if args.lookup_summary:
        recorded = json.loads(args.lookup_summary.read_text())
        start, complete = recorded["start"], recorded["complete"]
        assert start["family"] == "lookup" and start["k"] == 2
        assert start["normalized"] and start["exact_ranks"] and not start["random"]
        assert start["start"] == 0 and start["count"] == 1 << 25
        assert complete["tested"] == 1 << 25
        assert complete["best_rank"] == 28 and complete["histogram"][0] == [28, 12]
        assert sum(count for rank, count in complete["histogram"]) == 1 << 25
        assert len({entry["lookup_index"] for entry in lookup_winners}) == 12
        assert recorded["minimizers"] == lookup_winners
        assert recorded["native_source_sha256"] == hashlib.sha256(
            Path(__file__).with_name("explore_group_correction_rank.cpp").read_bytes()).hexdigest()
        lookup_summary = {"file": args.lookup_summary.name,
                          "sha256": hashlib.sha256(args.lookup_summary.read_bytes()).hexdigest(),
                          "scope": "full native run; bounded independent Python replay and all minimizers"}
    result = {"status": "native ranks independently match exact Python ranks",
              "complete_normalized_F4_quadratic_replay": args.full,
              "blocks": blocks,
              "lookup_minimizers": lookup_winners,
              "exhaustive_lookup_summary": lookup_summary,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "native_source_sha256": hashlib.sha256(Path(__file__).with_name(
                  "explore_group_correction_rank.cpp").read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
