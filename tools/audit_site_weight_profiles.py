#!/usr/bin/env python3
"""Audit reported site-weight profiles and retain every exact coefficient rank."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isclose, isqrt, log2, sqrt
from pathlib import Path

FILES = ["site-square-t10.json", "site-square-t12.json", "site-square-t14.json",
         "site-square-t16.json", "site-beatty-t12.json", "site-beatty-t16.json",
         "site-asymmetric-t12.json"]


def audit(record):
    t, weights = record["t"], record["weights"]
    assert len(weights) == t and all(isinstance(a, int) and a >= 0 for a in weights)
    total = sum(weights)
    reachable = [0] * (t + 1)
    reachable[0] = 1
    for i, a in enumerate(weights):
        for k in range(i + 1, 0, -1):
            reachable[k] |= reachable[k - 1] << a
    keys = {(k, u) for k, word in enumerate(reachable) for u in range(total + 1) if word >> u & 1}
    entries = {(k, u): (r, c) for k, u, r, c in record["pieces"]}
    assert len(entries) == len(record["pieces"]) == record["nonzero_pieces"]
    assert set(entries) == keys
    for (k, u), (r, c) in entries.items():
        assert isinstance(r, int) and isinstance(c, int) and 0 < r <= 4 ** t and 0 < c <= 4 ** t
        assert entries[t - k, total - u] == (c, r)
    products = [r * c for r, c in entries.values()]
    maximum, rank_sum = max(products), sum(sqrt(p) for p in products)
    assert maximum == record["maximum_rank_product"]
    assert isclose(rank_sum, record["rank_sum_estimate"], rel_tol=1e-12)
    assert isclose(log2(maximum) / (4 * t), record["maximum_product_rho"], abs_tol=1e-12)
    assert isclose(log2(rank_sum) / (2 * t), record["rank_sum_rho"], abs_tol=1e-12)
    lower = sum(isqrt(p) for p in products)
    upper = sum(isqrt(p) + (isqrt(p) ** 2 != p) for p in products)
    assert lower > 2 ** ((3 * t + 1) // 2)
    ordered = sorted(weights)
    symmetric = len({a + b for a, b in zip(ordered, reversed(ordered))}) == 1
    return {"rank_sum_integer_bounds": [lower, upper],
            "centrally_symmetric_weights": symmetric,
            "symmetry_pair_sum": ordered[0] + ordered[-1] if symmetric else None,
            "maximum_product_pieces": [row for row in record["pieces"] if row[2] * row[3] == maximum]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profiles = {}
    for name in FILES:
        path = args.directory / name
        profile = json.loads(path.read_text())
        profile["audit"] = audit(profile)
        profile["raw_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        profiles[name] = profile
    growth = []
    for a, b in zip((10, 12, 14), (12, 14, 16)):
        first, second = profiles[f"site-square-t{a}.json"], profiles[f"site-square-t{b}.json"]
        growth.append({"from_t": a, "to_t": b,
                       "incremental_max_product_rho": log2(second["maximum_rank_product"] / first["maximum_rank_product"]) / (4 * (b - a)),
                       "incremental_rank_sum_rho": log2(second["rank_sum_estimate"] / first["rank_sum_estimate"]) / (2 * (b - a))})
    result = {"status": "passed", "field": "F_2",
              "scope": "Exact finite profiles and integer rank-sum exclusions; no asymmetric asymptotic contraction claim.",
              "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "profiles": profiles, "square_weight_growth": growth}
    output = args.output or args.directory / "site-weight-profiles.json"
    # Retain complete exact ranks without expanding every four-integer row.
    output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"status": result["status"], "profiles": len(profiles),
                      "coefficient_pieces": sum(p["nonzero_pieces"] for p in profiles.values()),
                      "centrally_symmetric_profiles": [name for name, p in profiles.items() if p["audit"]["centrally_symmetric_weights"]]}, indent=2))


if __name__ == "__main__":
    main()
