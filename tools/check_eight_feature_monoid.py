#!/usr/bin/env python3
"""Check the order-dependent eight-feature monoid and its exact ranks."""
from __future__ import annotations
import argparse
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
from check_eight_feature_crossings import entry

MATRICES = [9, 2, 4, 9, 0, 0, 0, 0]


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def multiply(first, second):
    result = []
    for row in first:
        value = 0
        while row:
            low = row & -row
            row ^= low
            value ^= second[low.bit_length() - 1]
        result.append(value)
    return result


def coefficient_matrices(t):
    d = 1 << t
    rows = [[0] * (d * d) for _ in range(t + 1)]
    columns = [[0] * (d * d) for _ in range(t + 1)]
    for x, y, u in product(range(d), repeat=3):
        v = x ^ y ^ u
        polynomial = entry(t, (x, y, u, v), MATRICES, 0)
        for k in range(t + 1):
            if polynomial >> k & 1:
                rows[k][x * d + y] |= 1 << (v * d + u)
                columns[k][x * d + u] |= 1 << (y * d + v)
    return rows, columns


def check_auxiliary_moments():
    # A bit at 9*i+j is the monomial a**i b**j; degrees never exceed eight.
    count = 0
    for length in range(9):
        for word in product(range(2), repeat=length):  # 0=U, 1=V
            state = [1, 0, 0, 0]
            for letter in word:
                next_state = [0] * 4
                for index, polynomial in enumerate(state):
                    for factor, shift in ((0, 9), (1, 1)):
                        bit = (index >> factor) & 1
                        if bit == letter:
                            next_state[index ^ (1 << factor)] ^= polynomial << shift
                state = next_state
            expected = 0
            if length % 2 == 0 and all(letter == (i % 2) for i, letter in enumerate(word)):
                for i in range(length + 1):
                    if comb(length, i) % 2:
                        expected ^= 1 << (9 * i + length - i)
            assert state[0] == expected
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    products_checked = 0
    finite_profiles = []
    for t in range(1, 5):
        rows, columns = coefficient_matrices(t)
        n = 1 << (2 * t)
        for i, j in product(range(t + 1), repeat=2):
            assert multiply(rows[i], rows[j]) == (rows[i] if i == j else [0] * n)
            products_checked += 1
        for index in range(n):
            total = 0
            for matrix in rows:
                total ^= matrix[index]
            assert total == 1 << index
        expected = [(1 << t) * comb(t, k) for k in range(t + 1)]
        assert list(map(rank, rows)) == list(map(rank, columns)) == expected
        finite_profiles.append({"t": t, "row_ranks": expected, "column_ranks": expected})
    root = Path(__file__).resolve().parent
    evidence = root.parent / "evidence" / "eight-feature-extensions.json"
    if not evidence.exists():
        evidence = root / "eight-feature-extensions.json"
    saved = json.loads(evidence.read_text())
    model = next(a for a in saved["models"] if a["id"] == "model-03")
    assert model["matrices"] == MATRICES and model["off_shift"] == 0 and model["boundary"] == 0
    for profile in model["profiles"]:
        t = profile["t"]
        expected = [(1 << t) * comb(t, k) for k in range(t + 1)]
        assert profile["row_ranks"] == profile["column_ranks"] == expected
        assert profile["rank_sum_integer_lower_bound"] == profile["rank_sum_integer_upper_bound"] == 1 << (2 * t)
    result = {"status": "passed", "field": "F_2",
              "exact_coefficient_products_checked": products_checked,
              "formal_two_auxiliary_moments_checked": check_auxiliary_moments(),
              "independent_profiles": finite_profiles,
              "native_profile_lengths": [a["t"] for a in model["profiles"]],
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Finite verification of the all-length proof. The exact rank sum is 4**t; no improved factoring bound."}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
