#!/usr/bin/env python3
"""Check exact quadratic-pivot rank reduction in a binary square-zero algebra."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random


def terms(polynomial):
    while polynomial:
        bit = polynomial & -polynomial
        yield bit.bit_length() - 1
        polynomial ^= bit


def multiply(a, b):
    answer = 0
    for x in terms(a):
        for y in terms(b):
            if not x & y:
                answer ^= 1 << (x | y)
    return answer


def direct_rank(f, m):
    basis = {}
    for s in range(1 << m):
        column = sum(1 << (a | s) for a in terms(f) if not a & s)
        while column:
            pivot = column.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = column
                break
            column ^= basis[pivot]
    return len(basis)


def pivot_decomposition(f, m, i, j):
    others = [k for k in range(m) if k not in (i, j)]
    blocks = [0] * 4
    for support in terms(f):
        index = ((support >> i) & 1) + 2 * ((support >> j) & 1)
        reduced = sum(((support >> k) & 1) << q for q, k in enumerate(others))
        blocks[index] ^= 1 << reduced
    d, b, c, a = blocks
    assert a & 1
    assert not b & 1 and not c & 1 and not d & 1
    assert multiply(a, a) == 1
    assert multiply(b, b) == multiply(c, c) == 0
    g = d ^ multiply(multiply(a, b), c)
    assert not g & 1
    return g


def reduced_rank(f, m, counts=None):
    if not f:
        if counts is not None:
            counts["zero"] += 1
        return 0
    if f & 1:
        if counts is not None:
            counts["unit"] += 1
        return 1 << m
    if any(f >> (1 << i) & 1 for i in range(m)):
        if counts is not None:
            counts["linear"] += 1
        return 1 << (m - 1)
    for i in range(m):
        for j in range(i + 1, m):
            if f >> ((1 << i) | (1 << j)) & 1:
                if counts is not None:
                    counts["quadratic"] += 1
                g = pivot_decomposition(f, m, i, j)
                return (1 << (m - 2)) + 2 * reduced_rank(g, m - 2, counts)
    if counts is not None:
        counts["gaussian"] += 1
    return direct_rank(f, m)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = Counter()
    for f in range(1 << 16):
        assert reduced_rank(f, 4, counts) == direct_rank(f, 4)
    rng = random.Random(20261030)
    random_cases = 0
    for m in range(5, 10):
        for _ in range(12):
            f = rng.getrandbits(1 << m) & ~1
            for i in range(m):
                f &= ~(1 << (1 << i))
            assert reduced_rank(f, m) == direct_rank(f, m)
            random_cases += 1
    # The block recurrence requires g^2=0. The separate unit branch is
    # necessary: for g=1, multiplication by xy+1 has rank four, not three.
    assert direct_rank((1 << 3) | 1, 2) == 4
    # With a linear g the square-zero block recurrence is still valid.
    assert direct_rank((1 << 3) | (1 << 4), 3) == 4
    result = {
        "status": "passed", "all_four_variable_polynomials": 65536,
        "random_larger_cases": random_cases,
        "branch_counts_including_recursive_calls": dict(counts),
        "constant_and_linear_residual_guards": 2,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
