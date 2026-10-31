#!/usr/bin/env python3
"""Check binary hook-residue identities and full Schur central ranks.

This checker compiles the native packed-determinant implementation locally.
It needs a C++17 compiler, Python 3, and no third-party Python packages.
The optional weight-sixteen run proves a finite rank, not a factoring bound.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

from check_partition_channels import partitions
from schur_fermion import frobenius, schur_mode


ROOT = Path(__file__).resolve().parent.parent


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


def support_bits(support):
    return sum(1 << i for i in support)


def query_native(executable, queries):
    def encode(partition):
        return ",".join(map(str, partition)) or "-"

    text = "".join(" ".join(map(encode, query)) + "\n" for query in queries)
    completed = subprocess.run(
        [str(executable), "--queries"], input=text, text=True,
        capture_output=True, check=True,
    )
    values = list(map(int, completed.stdout.split()))
    assert len(values) == len(queries)
    assert set(values) <= {0, 1}
    return values


def check_modes(executable):
    small = [p for weight in range(6) for p in partitions(weight)]
    queries = list(product(small, repeat=3))
    rng = random.Random(20261030)
    larger = {w: partitions(w) for w in (8, 12, 16, 24, 32)}
    for _ in range(256):
        weights = [rng.choice(tuple(larger)) for _ in range(3)]
        queries.append(tuple(rng.choice(larger[w]) for w in weights))
    actual = query_native(executable, queries)
    for query, value in zip(queries, actual):
        assert value == schur_mode(*query, prime=2), (query, value)
    return {"all_weights_through_five": len(small) ** 3,
            "larger_deterministic_samples": 256,
            "total": len(queries)}


def check_vacuum(executable, h):
    labels = partitions(h)
    queries = [(x, (), y) for x in labels for y in labels]
    queries += [(x, y, ()) for x in labels for y in labels]
    values = query_native(executable, queries)
    d = len(labels)
    rows = [sum(values[x * d + y] << y for y in range(d)) for x in range(d)]
    rank = rank_binary(rows)
    assert rank == d, (h, rank, d)
    assert values[d * d:] == [int(x == y) for x in labels for y in labels]
    return {"vacuum_gram_rank": rank, "creation_identity_entries": d * d}


def check_weight(executable, h, full_h16):
    command = [str(executable), str(h)]
    if h == 16 and not full_h16:
        command.append("13000")
    completed = subprocess.run(command, text=True, capture_output=True, check=True)
    result = json.loads(completed.stdout)
    result.update(check_vacuum(executable, h))
    upper = len(partitions(2 * h)) - 3
    key = "central_rank" if result["all_pairs"] else "central_rank_lower_bound"
    measured = result[key]
    assert measured <= upper
    expected = {2: 2, 4: 12, 8: 228, 16: 8346}[h]
    assert measured == expected, (h, measured, expected)
    result["all_weight_hook_parity_upper_bound"] = upper
    result["certified_full_central_rank"] = (
        measured if result["all_pairs"] or measured == upper else None
    )
    result["central_alone_exceeds_seed_threshold_by_rank_squared_gt_d_cubed"] = (
        measured ** 2 > result["dimension"] ** 3
    )

    middle = ((),) + partitions(2 * h)
    parity_supports = [
        [i for i, p in enumerate(middle)
         if p and len(frobenius(p)[0]) == 1 and frobenius(p)[0][0] % 2 == parity]
        for parity in (0, 1)
    ]
    parity_bits = list(map(support_bits, parity_supports))
    assert all((a & b).bit_count() % 2 == 0 for a in parity_bits for b in parity_bits)
    result["hook_parity_supports"] = parity_supports
    for name in ("left_perpendicular_supports", "right_perpendicular_supports"):
        if name not in result:
            continue
        actual = list(map(support_bits, result[name]))
        assert rank_binary(actual + parity_bits) == rank_binary(actual)
        if h >= 8:
            assert rank_binary(actual) == 2
    return result


def negative_control():
    # Residue classes modulo four are not universally annihilating.
    # For (2,2), the arm-zero residue charge is one.
    label = (2, 2)
    value = 0
    for arm in range(0, 8, 4):
        hook = (1,) * (7 - arm) + (arm + 1,)
        value ^= schur_mode(label, label, hook, prime=2)
    assert value == 1
    return {"weight": 4, "modulus": 4, "residue": 0,
            "external_label": list(label), "nonzero_hook_sum": value}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--large", action="store_true", help="include weight sixteen")
    parser.add_argument("--full-h16", action="store_true",
                        help="scan all weight-sixteen pairs, rather than a sufficient subset")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        raise SystemExit(f"C++17 compiler not found: {args.compiler}")
    compiler_version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    with tempfile.TemporaryDirectory(prefix="schur-central-") as temp:
        executable = Path(temp) / "schur-central"
        subprocess.run([compiler, "-std=c++17", "-O3",
                        str(ROOT / "tools/schur_homogeneous_central.cpp"),
                        "-o", str(executable)], check=True)
        modes = check_modes(executable)
        weights = [2, 4, 8] + ([16] if args.large or args.full_h16 else [])
        profiles = [check_weight(executable, h, args.full_h16) for h in weights]
    files = ["tools/schur_fermion.py", "tools/schur_fermion_core.hpp",
             "tools/schur_homogeneous_central.cpp",
             "tools/check_schur_homogeneous_central.py"]
    result = {
        "status": "exact finite ranks and all-weight hook-residue theorem; no new factoring exponent",
        "compiler": compiler_version,
        "compile_flags": ["-std=c++17", "-O3"],
        "mode_cross_checks": modes,
        "profiles": profiles,
        "negative_control": negative_control(),
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in files},
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
