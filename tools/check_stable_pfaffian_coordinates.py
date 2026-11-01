#!/usr/bin/env python3
"""Check complete stable pencil/chart bounds and reproduce their finite grids."""

from itertools import combinations, product
import argparse
import hashlib
import json
from pathlib import Path

from explore_complement_edge_pencils import all_length_labels
from pfaffian_particle_hole import rank_profile
from stable_pfaffian_charts import encode_linear, score_chart
from stable_pfaffian_pencils import score_pencil


TRIPLES = tuple((p, q, r) for (p, q), r in product(combinations((3, 5, 6), 2), (1, 2, 4, 7)))


def verify_small():
    checks = []
    for pairs, length in ((2, 2), (3, 3)):
        labels = all_length_labels(pairs, length)
        for p, q, r in TRIPLES:
            actual = rank_profile(labels, p=p, q=q, r=r)
            upper = score_pencil(pairs, length, p, q, r)
            bounds = {degree: (a, b) for degree, a, b in upper["profile"]}
            assert all(a <= bounds.get(degree, (0, 0))[0] and b <= bounds.get(degree, (0, 0))[1]
                       for degree, a, b in actual["profile"])
            stride = 2 * pairs + 3
            actual_chart = rank_profile(labels, p=encode_linear(p, stride),
                                         q=encode_linear(q, stride), r=encode_linear(r, stride))
            upper_chart = score_chart(pairs, length, p, q, r)
            chart_bounds = {(a, b): (row, col) for a, b, row, col in upper_chart["profile"]}
            for degree, row, col in actual_chart["profile"]:
                a, b = degree % stride, degree // stride
                ur, uc = chart_bounds.get((a, b), (0, 0))
                assert row <= ur and col <= uc
            checks.append({"pairs": pairs, "length": length, "parameter_bits": [p, q, r],
                           "complete_univariate_coefficients_checked": len(actual["profile"]),
                           "complete_bivariate_coefficients_checked": len(actual_chart["profile"])})
    return checks


def compact(report):
    result = dict(report)
    profile = result.pop("profile")
    result["complete_coefficient_count"] = len(profile)
    result["complete_profile_sha256"] = hashlib.sha256(
        json.dumps(profile, separators=(",", ":")).encode()).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "complete endpoint-preserving coordinate evaluations with monomial ancillary entries; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pencil_module_sha256": hashlib.sha256(Path(__file__).with_name("stable_pfaffian_pencils.py").read_bytes()).hexdigest(),
        "chart_module_sha256": hashlib.sha256(Path(__file__).with_name("stable_pfaffian_charts.py").read_bytes()).hexdigest(),
        "actual_complete_profile_controls": verify_small(),
        "pencil_grid": [], "chart_grid": [],
    }
    for pairs, length in ((4, 2), (8, 4), (16, 8), (32, 15)):
        pencils = [compact(score_pencil(pairs, length, p, q, r)) for p, q, r in TRIPLES]
        charts = [compact(score_chart(pairs, length, p, q, r)) for p, q, r in TRIPLES]
        report["pencil_grid"].extend(pencils)
        report["chart_grid"].extend(charts)
        print(json.dumps({"pairs": pairs, "length": length,
                          "best_pencil": min(pencils, key=lambda row: row["rho_of_upper_bound"]),
                          "best_chart": min(charts, key=lambda row: row["rho_of_upper_bound"])}), flush=True)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
