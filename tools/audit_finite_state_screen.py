#!/usr/bin/env python3
"""Summarize exact finite-state profiles without discarding the Pareto frontier."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isclose, isqrt, log2, sqrt
from pathlib import Path


SCREENS = [
    ("markov-cost-screen-t5.jsonl", "cost", 6144, 2),
    ("markov-cost-screen-t6.jsonl", "cost", 26112, 5),
    ("markov-polynomial-screen-012-t6.jsonl", "polynomial", 19683, None),
    ("markov-polynomial-screen-001-t6.jsonl", "polynomial", 19683, None),
    *[(f"auxiliary-matrix-screen-{a}-b{b}-t6.jsonl", "auxiliary", 4096, None)
      for a in ("001", "012") for b in (0, 1)],
]

PROFILES = [
    "markov-run-cost-t8.json", "markov-run-cost-t10.json",
    "markov-domainwall-cost-t8.json", "markov-domainwall-cost-t10.json",
    "markov-domainwall-cost-t12.json", "markov-run-count-bigrading-t10.json",
    "markov-run-count-bigrading-t12.json", "markov-two-pair-baseline-t8.json",
]


def validate(record):
    t = record["t"]
    rows, columns = record["row_ranks"], record["column_ranks"]
    assert len(rows) == len(columns)
    assert all(isinstance(a, int) and 0 <= a <= 4 ** t for a in rows + columns)
    assert all(bool(a) == bool(b) for a, b in zip(rows, columns))
    products = [a * b for a, b in zip(rows, columns)]
    maximum = max(products)
    total = sum(sqrt(p) for p in products)
    assert maximum == record["maximum_rank_product"]
    assert isclose(total, record["rank_sum_estimate"], rel_tol=1e-12)
    assert isclose(log2(maximum) / (4 * t), record["maximum_product_rho"], abs_tol=1e-12)
    assert isclose(log2(total) / (2 * t), record["rank_sum_rho"], abs_tol=1e-12)
    return sum(isqrt(p) for p in products)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summaries, frontiers, profiles = [], {}, {}
    seeds = {
        "run-t6": ([0, 0, 2], [0] * 8 + [1]),
        "domainwall-t6": ([0, 0, 5], [0, 0, 1, 0, 0, 1, 1, 1, 0]),
    }
    for name, kind, expected, vertex_max in SCREENS:
        path = args.directory / name
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        assert len(rows) == expected
        lower = min(map(validate, rows))
        # This integer inequality alone excludes the strict finite rank-sum
        # criterion S < 2^(sqrt(2) t), since sqrt(2) < 3/2.
        assert lower > 2 ** ((3 * rows[0]["t"] + 1) // 2)
        frontier, best = [], float("inf")
        for r in sorted(rows, key=lambda x: (x["maximum_rank_product"], x["rank_sum_estimate"])):
            if r["rank_sum_estimate"] < best - 1e-8:
                frontier.append(r)
                best = r["rank_sum_estimate"]
        record = {
            "file": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "cases": len(rows), "t": rows[0]["t"], "kind": kind,
            "minimum_integer_rank_sum_lower_bound": lower,
            "max_product_best": min(r["maximum_rank_product"] for r in rows),
            "rank_sum_best_estimate": min(r["rank_sum_estimate"] for r in rows),
            "pareto_representatives": len(frontier),
        }
        if vertex_max is not None:
            record["vertex_max"] = vertex_max
        summaries.append(record)
        frontiers[name] = frontier
        if name == "markov-cost-screen-t6.jsonl":
            for label, (vertex, edge) in seeds.items():
                profiles[label] = next(r for r in rows if r["vertex"] == vertex and r["edge"] == edge)
    for name in PROFILES:
        profiles[name] = json.loads((args.directory / name).read_text())
        validate(profiles[name])
    for name in ("markov-domainwall-cost-t14.json",):
        path = args.directory / name
        if path.exists():
            profiles[name] = json.loads(path.read_text())
            validate(profiles[name])
    growth = {}
    for family, names in {
        "domainwall": ["domainwall-t6", "markov-domainwall-cost-t8.json", "markov-domainwall-cost-t10.json", "markov-domainwall-cost-t12.json", "markov-domainwall-cost-t14.json"],
        "run": ["run-t6", "markov-run-cost-t8.json", "markov-run-cost-t10.json"],
        "run_count_bigrading": ["markov-run-count-bigrading-t10.json", "markov-run-count-bigrading-t12.json"],
    }.items():
        values = [profiles[name] for name in names if name in profiles]
        growth[family] = [
            {"from_t": a["t"], "to_t": b["t"],
             "incremental_max_product_rho": log2(b["maximum_rank_product"] / a["maximum_rank_product"]) / (4 * (b["t"] - a["t"])),
             "incremental_rank_sum_rho": log2(b["rank_sum_estimate"] / a["rank_sum_estimate"]) / (2 * (b["t"] - a["t"]))}
            for a, b in zip(values, values[1:])
        ]
    result = {
        "status": "passed",
        "scope": "Exact finite model screen; no asymptotic rank claim or improved factoring exponent. Decimal square-root metrics select candidates only; finite exclusions use integer lower bounds.",
        "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "total_cases": sum(r["cases"] for r in summaries),
        "screens": summaries, "pareto_profiles": frontiers,
        "selected_larger_profiles": profiles, "growth_increments": growth,
    }
    output = args.output or args.directory / "finite-state-pairing-screen.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "cases": result["total_cases"],
                      "screens": summaries, "growth_increments": growth}, indent=2))


if __name__ == "__main__":
    main()
