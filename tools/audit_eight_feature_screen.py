#!/usr/bin/env python3
"""Audit reproducible structural screens of full local parity matrix products."""
from __future__ import annotations
import argparse
import hashlib
from itertools import product, zip_longest
import json
import math
from pathlib import Path
from make_eight_feature_candidates import DIAGONALS, OFF_DIAGONALS, efficient_diagonals, matrices
from audit_parity_matrix_screen import PRODUCT
from reproduce_eight_feature_extensions import MODELS

SCREENS = [
    ("eight-feature-screen-b0-t4.jsonl", "broad", 0, 0),
    ("eight-feature-screen-b1-t4.jsonl", "broad", 1, 0),
    ("eight-feature-efficient-b0-shift0-t4.jsonl", "efficient", 0, 0),
    ("eight-feature-efficient-b0-shift1-t4.jsonl", "efficient", 0, 1),
    ("eight-feature-sparse-b0-shift2-t4.jsonl", "sparse", 0, 2),
]


def candidates(pool):
    diagonals = efficient_diagonals() if pool == "efficient" else DIAGONALS
    off = OFF_DIAGONALS if pool == "broad" else [1 << i for i in range(8)]
    result = []
    for parameters in product(diagonals, diagonals, off, off):
        values = matrices(parameters)
        if values[:4] == values[4:]:
            continue
        if all(PRODUCT[a][b] == PRODUCT[b][a] for a in values for b in values):
            continue
        result.append(parameters)
    assert len(result) == {"broad": 8246, "sparse": 4030, "efficient": 172480}[pool]
    return result


def normalized(profile):
    return [profile["t"], profile["boundary"], profile.get("off_shift", 0),
            profile["matrices"], profile["row_ranks"], profile["column_ranks"]]


def summary_profile(profile):
    return {k: profile[k] for k in ("t", "boundary", "matrices", "row_ranks", "column_ranks",
            "maximum_rank_product", "rank_sum_estimate", "maximum_product_rho", "rank_sum_rho")} | {"off_shift": profile.get("off_shift", 0)}


def audit(directory):
    scale = 10 ** 12
    pools = {name: candidates(name) for name in ("broad", "sparse", "efficient")}
    checked = []
    unique_models = set()
    wanted = {(boundary, shift, tuple(values)): name for name, boundary, shift, values in MODELS}
    matched = {}
    for filename, pool, boundary, shift in SCREENS:
        path = directory / filename
        canonical_hash = hashlib.sha256()
        best_product = best_sum = None
        minimum_integer_lower = None
        minimum_scaled_lower = minimum_scaled_upper = None
        count = local_crossing_only = 0
        profiles = {}
        with path.open() as stream:
            for parameters, line in zip_longest(pools[pool], stream):
                assert parameters is not None and line is not None, (filename, count, "record count")
                profile = json.loads(line)
                assert profile["t"] == 4 and profile["boundary"] == boundary and profile.get("off_shift", 0) == shift
                assert profile["matrices"] == matrices(parameters)
                r, c = profile["row_ranks"], profile["column_ranks"]
                assert len(r) == len(c) == 4 * (1 + shift) + 1
                assert all(isinstance(a, int) and 0 <= a <= 256 for a in r + c)
                assert all(bool(a) == bool(b) for a, b in zip(r, c))
                products = [a * b for a, b in zip(r, c)]
                assert max(products) == profile["maximum_rank_product"]
                roots = [math.isqrt(p * scale * scale) for p in products]
                lower = sum(roots)
                upper = lower + sum(v * v != p * scale * scale for v, p in zip(roots, products))
                estimate = profile["rank_sum_estimate"]
                assert lower / scale - 1e-10 <= estimate <= upper / scale + 1e-10
                integer_lower = sum(map(math.isqrt, products))
                minimum_integer_lower = integer_lower if minimum_integer_lower is None else min(minimum_integer_lower, integer_lower)
                minimum_scaled_lower = lower if minimum_scaled_lower is None else min(minimum_scaled_lower, lower)
                minimum_scaled_upper = upper if minimum_scaled_upper is None else min(minimum_scaled_upper, upper)
                key = (boundary, shift, tuple(profile["matrices"]))
                unique_models.add(key)
                if key in wanted:
                    matched[wanted[key]] = normalized(profile)
                local_crossing_only += not any(profile["matrices"][i] for i in (1, 2, 5, 6))
                canonical_hash.update((json.dumps(normalized(profile), separators=(",", ":")) + "\n").encode())
                if best_product is None or (profile["maximum_rank_product"], estimate) < (best_product["maximum_rank_product"], best_product["rank_sum_estimate"]):
                    best_product = profile
                if best_sum is None or estimate < best_sum["rank_sum_estimate"]:
                    best_sum = profile
                profiles.setdefault((tuple(r), tuple(c)), profile)
                count += 1
        frontier = []
        least_sum = float("inf")
        for profile in sorted(profiles.values(), key=lambda a: (a["maximum_rank_product"], a["rank_sum_estimate"])):
            if profile["rank_sum_estimate"] < least_sum - 1e-10:
                frontier.append(summary_profile(profile))
                least_sum = profile["rank_sum_estimate"]
        assert minimum_integer_lower > 64  # 64 > 2**(sqrt(2)*4), using 2*4**2 < 6**2.
        checked.append({"raw_file": filename, "candidate_pool": pool, "boundary": boundary,
            "off_shift": shift, "checked_cases": count, "distinct_exact_rank_profiles": len(profiles),
            "local_crossing_only_cases": local_crossing_only,
            "minimum_integer_rank_sum_lower_bound": minimum_integer_lower,
            "minimum_rank_sum_interval_scale": scale,
            "minimum_rank_sum_scaled_lower": minimum_scaled_lower,
            "minimum_rank_sum_scaled_upper": minimum_scaled_upper,
            "best_maximum_product_profile": summary_profile(best_product),
            "best_rank_sum_profile": summary_profile(best_sum),
            "numerical_pareto_representatives": frontier,
            "raw_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "normalized_rank_sha256": canonical_hash.hexdigest()})
        print(json.dumps({"audited": filename, "cases": count, "minimum_integer_rank_sum_lower_bound": minimum_integer_lower}), flush=True)
    assert set(matched) == {model[0] for model in MODELS}
    return {"scope": "Exact finite structural screens over F_2. Pools overlap; cases are not all distinct. Numerical Pareto choices are selection diagnostics, while rank-product and integer rank-sum bounds are exact. No asymptotic contraction claim.",
            "checked_screen_cases": sum(a["checked_cases"] for a in checked),
            "distinct_parameter_models": len(unique_models),
            "screens": checked, "selected_t4_profiles": matched,
            "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.directory)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
