#!/usr/bin/env python3
"""Audit both exhaustive four-feature matrix-product screens and retain Pareto data."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isclose, isqrt, sqrt
from pathlib import Path


def multiply(a, b):
    return sum((sum(((a >> (2 * i + k)) & 1) * ((b >> (2 * k + j)) & 1)
                    for k in range(2)) % 2) << (2 * i + j)
               for i in range(2) for j in range(2))


PRODUCT = [[multiply(a, b) for b in range(16)] for a in range(16)]


def boundary_value(a, boundary):
    return (a & 1) if boundary == 0 else ((a >> 1) ^ (a >> 2) ^ (a >> 3)) & 1


def order_witness(matrices, boundary):
    for i in range(4):
        for j in range(i + 1, 4):
            commutator = PRODUCT[matrices[i]][matrices[j]] ^ PRODUCT[matrices[j]][matrices[i]]
            if boundary_value(commutator, boundary):
                return [[i, j], [j, i]]
            for k in range(4):
                if boundary_value(PRODUCT[commutator][matrices[k]], boundary):
                    return [[i, j, k], [j, i, k]]
                if boundary_value(PRODUCT[matrices[k]][commutator], boundary):
                    return [[k, i, j], [k, j, i]]
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "status": "passed", "field": "F_2", "total_cases": 131072,
        "scope": "Exact finite screen; no asymptotic rank or factoring-score claim. Order dependence is certified by a length-two or length-three word witness when supplied; the classes overlap.",
        "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "screens": [],
    }
    for boundary in (0, 1):
        path = args.directory / f"parity-matrix-screen-b{boundary}-t6.jsonl"
        records = [json.loads(line) for line in path.read_text().splitlines()]
        assert len(records) == 65536
        assert len({tuple(record["matrices"]) for record in records}) == 65536
        classes = {name: [] for name in ("all", "nonzero_mixed_feature",
                   "witnessed_order_dependence", "mixed_and_witnessed",
                   "mixed_ordered_with_identity_padding")}
        for record in records:
            assert record["boundary"] == boundary and record["t"] == 6
            assert len(record["row_ranks"]) == len(record["column_ranks"]) == 7
            assert all(isinstance(rank, int) and 0 <= rank <= 4096
                       for rank in record["row_ranks"] + record["column_ranks"])
            assert all(bool(r) == bool(c) for r, c in zip(record["row_ranks"], record["column_ranks"]))
            products = [r * c for r, c in zip(record["row_ranks"], record["column_ranks"])]
            assert max(products) == record["maximum_rank_product"]
            assert isclose(sum(sqrt(p) for p in products), record["rank_sum_estimate"], rel_tol=1e-12)
            record["integer_rank_sum_lower_bound"] = sum(isqrt(p) for p in products)
            assert record["integer_rank_sum_lower_bound"] > 512
            matrices = record["matrices"]
            mixed = matrices[0] ^ matrices[1] ^ matrices[2] ^ matrices[3]
            witness = order_witness(matrices, boundary)
            classes["all"].append(record)
            if mixed:
                classes["nonzero_mixed_feature"].append(record)
            if witness:
                record["order_witness"] = witness
                classes["witnessed_order_dependence"].append(record)
            if mixed and witness:
                classes["mixed_and_witnessed"].append(record)
                if matrices[0] == 0 or matrices[3] == 0:
                    record["identity_padding_letter"] = 0 if matrices[0] == 0 else 3
                    classes["mixed_ordered_with_identity_padding"].append(record)
        summary = {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                   "boundary": boundary, "classes": {}}
        for name, records in classes.items():
            frontier, best = [], float("inf")
            for record in sorted(records, key=lambda r: (r["maximum_rank_product"], r["rank_sum_estimate"])):
                if record["rank_sum_estimate"] < best - 1e-8:
                    frontier.append(record)
                    best = record["rank_sum_estimate"]
            summary["classes"][name] = {
                "count": len(records),
                "minimum_integer_rank_sum_lower_bound": min(r["integer_rank_sum_lower_bound"] for r in records),
                "best_maximum_product": min(r["maximum_rank_product"] for r in records),
                "best_rank_sum_estimate": min(r["rank_sum_estimate"] for r in records),
                "pareto_profiles": frontier,
            }
        result["screens"].append(summary)
    output = args.output or args.directory / "parity-matrix-screen.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "cases": result["total_cases"],
                      "global_minimum_integer_rank_sum_lower_bound": min(s["classes"]["all"]["minimum_integer_rank_sum_lower_bound"] for s in result["screens"])}, indent=2))


if __name__ == "__main__":
    main()
