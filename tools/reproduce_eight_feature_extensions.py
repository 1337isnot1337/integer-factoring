#!/usr/bin/env python3
"""Reproduce selected eight-feature matrix-product crossing profiles."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

MODELS = [('model-00', 0, 0, [9, 6, 9, 8, 1, 0, 1, 0]), ('model-01', 0, 0, [9, 1, 6, 0, 9, 1, 0, 0]), ('model-02', 0, 1, [9, 8, 0, 15, 0, 0, 0, 0]), ('model-03', 0, 0, [9, 2, 4, 9, 0, 0, 0, 0]), ('model-04', 0, 0, [0, 0, 0, 0, 0, 0, 2, 6]), ('model-05', 0, 0, [9, 9, 6, 8, 1, 1, 0, 0]), ('model-06', 0, 1, [15, 8, 0, 9, 0, 0, 0, 0]), ('model-07', 0, 2, [9, 2, 4, 9, 0, 0, 0, 0]), ('model-08', 0, 2, [1, 1, 6, 0, 9, 1, 0, 0]), ('model-09', 1, 0, [11, 11, 2, 2, 15, 15, 6, 6]), ('model-10', 1, 0, [9, 13, 11, 15, 4, 0, 6, 2])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lengths", nargs="+", type=int, default=[4, 6, 8])
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory / "eight_feature_crossings"
    source = root / "eight_feature_crossings.cpp"
    subprocess.run(["g++", "-O3", "-std=c++17", str(source), "-o", str(binary)], check=True)
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    results = []
    for name, boundary, shift, matrices in MODELS:
        profiles = []
        for t in sorted(set(args.lengths)):
            output = directory / f"{name}-t{t}.json"
            command = [str(binary), "profile", str(t), str(output), str(boundary),
                       *map(str, matrices), str(shift)]
            if not (args.resume and output.exists()):
                subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
            profile = json.loads(output.read_text())
            assert (profile["t"], profile["boundary"], profile["off_shift"], profile["matrices"]) == (t, boundary, shift, matrices)
            assert len(profile["row_ranks"]) == len(profile["column_ranks"]) == t * (1 + shift) + 1
            products = [r * c for r, c in zip(profile["row_ranks"], profile["column_ranks"])]
            assert max(products) == profile["maximum_rank_product"]
            profile["rank_sum_integer_lower_bound"] = sum(map(math.isqrt, products))
            profile["rank_sum_integer_upper_bound"] = sum(math.isqrt(p) + (math.isqrt(p) ** 2 < p) for p in products)
            profile["native_profile_sha256"] = hashlib.sha256(output.read_bytes()).hexdigest()
            profiles.append(profile)
            print(json.dumps({"model": name, "t": t, "maximum_product_rho": profile["maximum_product_rho"],
                              "rank_sum_rho": profile["rank_sum_rho"], "verified_native_profile": True}), flush=True)
        increments = []
        for first, second in zip(profiles, profiles[1:]):
            dt = second["t"] - first["t"]
            increments.append({"from_t": first["t"], "to_t": second["t"],
                "maximum_product_increment_rho": math.log2(second["maximum_rank_product"] / first["maximum_rank_product"]) / (4 * dt),
                "rank_sum_increment_rho": math.log2(second["rank_sum_estimate"] / first["rank_sum_estimate"]) / (2 * dt)})
        results.append({"id": name, "boundary": boundary, "off_shift": shift, "matrices": matrices,
                        "profiles": profiles, "growth_increments": increments})
        checkpoint = {"scope": "Finite exact ranks and numerical growth comparisons; no asymptotic contraction claim.",
                      "source_sha256": source_hash, "completed_models": len(results), "models": results}
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
        temporary.replace(args.output)


if __name__ == "__main__":
    main()
