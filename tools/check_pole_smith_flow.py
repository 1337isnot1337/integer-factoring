#!/usr/bin/env python3
"""Independent exact controls for the residue-tree Cauchy Smith algorithm."""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import random

from check_pole_smith_stabilization import matrix_at, smith_orders
from pole_smith_flow import smith_profile, smith_profile_heap


def residue_energy(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int) -> int:
    total = 0
    for level in range(1, k + 1):
        modulus = 1 << level
        balance: dict[int, int] = {}
        for a in rows:
            residue = a % modulus
            balance[residue] = balance.get(residue, 0) + 1
        for b in cols:
            residue = n * pow(b, -1, modulus) % modulus
            balance[residue] = balance.get(residue, 0) - 1
        squares = sum(value * value for value in balance.values())
        assert squares % 2 == 0
        total += squares // 2
    return total


def brute_contents(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int) -> list[int]:
    return [0] + [min(residue_energy(n, aa, bb, k)
                      for aa in combinations(rows, size) for bb in combinations(cols, size))
                  for size in range(1, min(len(rows), len(cols)) + 1)]


def v2_int(value: int) -> int:
    assert value
    value = abs(value)
    return (value & -value).bit_length() - 1


def direct_content(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int) -> int:
    assert len(rows) == len(cols)
    answer = k * len(rows)
    answer += sum(v2_int(b - a) for values in (rows, cols) for a, b in combinations(values, 2))
    answer -= sum(k if n == a * b else v2_int(n - a * b) for a in rows for b in cols)
    return answer


def odd_progression_content(n: int, first_row: int, first_col: int, size: int, k: int) -> int:
    """Independent floor-count evaluation, linear in the explicit row count."""
    assert first_row % 2 and first_col % 2 and n % 2
    factorial_content = 0
    modulus = 2
    while modulus < size:
        quotient, remainder = divmod(size, modulus)
        factorial_content += modulus * quotient * (quotient - 1) // 2 + quotient * remainder
        modulus *= 2
    vandermonde = size * (size - 1) // 2 + factorial_content
    cross = 0
    inverse_modulus = 1 << max(1, k - 1)
    for index in range(size):
        a = first_row + 2 * index
        half_delta = (n - a * first_col) // 2
        residue = half_delta * pow(a, -1, inverse_modulus) % inverse_modulus
        cross += size
        modulus = 2
        for _ in range(2, k + 1):
            root = residue % modulus
            if root < size:
                cross += 1 + (size - 1 - root) // modulus
            modulus *= 2
    return k * size + 2 * vandermonde - cross


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(810215)
    gaussian_checks = 0
    exhaustive_subset_checks = 0
    records = []
    for index in range(2000):
        n = rng.randrange(1, 250, 2)
        rows = tuple(sorted(rng.sample(range(1, 32, 2), rng.randrange(1, 7))))
        cols = tuple(sorted(rng.sample(range(1, 32, 2), rng.randrange(1, 7))))
        flow = smith_profile(n, rows, cols)
        direct = smith_orders(matrix_at(n, rows, cols, flow["k"]))
        assert flow["orders"] == direct
        gaussian_checks += 1
        if index < 500:
            assert flow["orders"] == smith_profile_heap(n, rows, cols)["orders"]
        if index < 200:
            assert flow["contents"] == brute_contents(n, rows, cols, flow["k"])
            exhaustive_subset_checks += 1
        records.append((n, rows, cols, flow["orders"]))
    large = []
    for size in (16, 64, 256, 1024, 16384, 65536):
        first = (1 << 100) + 3
        second = (1 << 100) + 7
        n = first * second
        rows = tuple(first + 2 * i for i in range(size))
        cols = tuple(second + 2 * i for i in range(size))
        flow = smith_profile(n, rows, cols)
        assert flow["orders"].count(0) == 1
        if size <= 1024:
            assert flow["contents"][-1] == direct_content(n, rows, cols, flow["k"])
        assert flow["contents"][-1] == odd_progression_content(n, first, second, size, flow["k"])
        large.append({"size": size, "input_bits": n.bit_length(), "k": flow["k"],
                      "compressed_nodes": flow["compressed_nodes"],
                      "uncompressed_nodes": flow["uncompressed_nodes"],
                      "vertex_updates": flow["vertex_updates"],
                      "vertex_update_bound": flow["compressed_nodes"] + size * (2 * flow["k"] + 1),
                      "factor_count": flow["orders"].count(0),
                      "maximal_minor_content": flow["contents"][-1],
                      "first_orders": flow["orders"][:8],
                      "last_orders": flow["orders"][-8:],
                      "profile_sha256": hashlib.sha256(json.dumps(flow["orders"]).encode()).hexdigest()})
    high_k = 10 ** 80
    base = smith_profile(105, (3, 5, 7, 11), (3, 5, 7, 15, 21))
    high = smith_profile(105, (3, 5, 7, 11), (3, 5, 7, 15, 21), high_k)
    assert high["orders"] == [value + high_k - base["k"] if value else 0 for value in base["orders"]]
    assert high["compressed_nodes"] == base["compressed_nodes"]
    assert smith_profile(15, (), (3, 5))["orders"] == []
    assert smith_profile(15, (3, 5), ())["orders"] == []
    result = {
        "claim": "Exact convex residue-tree algorithm for Cauchy pole Smith profiles",
        "dense_Qw_profile_checks": gaussian_checks,
        "heap_reference_profile_checks": 500,
        "exhaustive_subset_energy_checks": exhaustive_subset_checks,
        "random_seed": 810215,
        "small_profile_sha256": hashlib.sha256(json.dumps(records).encode()).hexdigest(),
        "large_explicit_controls": large,
        "binary_displacement_control": {"k": high_k, "base_precision": high["base_precision"],
                                        "compressed_nodes": high["compressed_nodes"],
                                        "orders": high["orders"]},
        "scope": "Explicit node lists. No compressed-interval evaluator and no improved factoring exponent.",
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
