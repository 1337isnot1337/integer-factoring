#!/usr/bin/env python3
"""Exact abstract-profile controls; no tensor realization is asserted."""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path


def profile(d: int) -> list[tuple[int, int]]:
    assert d >= 2 and d & (d - 1) == 0
    result: list[tuple[int, int]] = []
    for j in range(d.bit_length()):
        t = 1 << j
        assert d % t == 0
        result.extend(((d * t, d // t), (d // t, d * t)))
    return result


def flagged_complement(d: int, m: int) -> int:
    if m == 0:
        return d * d
    if m <= 2:
        return d * (d - 1) // 2
    slope = m - 1
    s = d // (2 * slope)
    return max(0, s * (d - slope * s), (s + 1) * (d - slope * (s + 1)))


def flagged_subset_checks(d: int, pairs: list[tuple[int, int]]) -> int:
    """Exact knapsack: every subset with summed small cut at most d+1."""
    cap = d + 1
    checked = 0
    for exchange in (False, True):
        oriented = [(c, r) if exchange else (r, c) for r, c in pairs]
        total = sum(r for r, _ in oriented)
        removed = [-1] * (cap + 1)
        removed[0] = 0
        for r, c in oriented:
            for budget in range(cap, c - 1, -1):
                if removed[budget - c] >= 0:
                    removed[budget] = max(removed[budget], removed[budget - c] + r)
        assert min(removed) >= 0
        for m, maximum in enumerate(removed):
            assert total - maximum >= flagged_complement(d, m)
            checked += 1
    return checked


def check_dimension(d: int) -> dict[str, int]:
    pairs = profile(d)
    m, n = d.bit_length() - 1, d * d
    assert len(pairs) == 2 * (m + 1)
    assert pairs.count((d, d)) == 2
    assert all(1 <= r <= n and 1 <= c <= n and r * c == n for r, c in pairs)
    r_total = sum(r for r, _ in pairs)
    c_total = sum(c for _, c in pairs)
    assert r_total == c_total == 2 * n + d - 1
    cost = len(pairs) * d
    assert cost == 2 * d * (m + 1)
    for s in range(d + 1):
        assert sum(min(r, s * s * c) for r, c in pairs) >= s * d
        assert sum(min(s * s * r, c) for r, c in pairs) >= s * d
        if 0 < s < d:
            t = 1 << (s.bit_length() - 1)
            assert t <= s < 2 * t <= d
            assert d * t + min(2 * d * t, d * s * s // (2 * t)) >= d * s
        elif s == d:
            assert (n, 1) in pairs and min(n, s * s) == s * d
    complement_values = [flagged_complement(d, m) for m in range(d + 2)]
    for s in range(d + 1):
        budget = min(f + s * s * m for m, f in enumerate(complement_values))
        assert sum(min(r, s * s * c) for r, c in pairs) >= budget
        assert sum(min(s * s * r, c) for r, c in pairs) >= budget
    flagged_subsets = flagged_subset_checks(d, pairs)
    eligible = 0
    for exchange in (False, True):
        oriented = [(c, r) if exchange else (r, c) for r, c in pairs]
        assert sum(c == 1 for _, c in oriented) == 1
        # Every coordinate is positive; cardinalities above two cannot qualify.
        for size in range(3):
            for indices in combinations(range(len(oriented)), size):
                small_cut = sum(oriented[i][1] for i in indices)
                if small_cut > 2:
                    continue
                eligible += 1
                assert size <= 1
                remaining = sum(r for i, (r, _) in enumerate(oriented) if i not in indices)
                if not indices:
                    assert remaining >= n
                else:
                    assert remaining >= n + d - 1 >= d * (d - 1) // 2
                if small_cut == 1:
                    assert oriented[indices[0]] == (n, 1)
                    assert all(remaining >= n - h * (h + 1) // 2 for h in range(d + 1))
    return {"d": d, "pieces": len(pairs), "cost": cost,
            "total_row_rank": r_total, "total_column_rank": c_total,
            "weighted_inequalities": 2 * (d + 1),
            "small_cut_subsets": eligible,
            "flagged_weighted_inequalities": 2 * (d + 1),
            "flagged_subset_budget_classes": flagged_subsets}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    dimensions = [check_dimension(1 << m) for m in range(1, 13)]
    example = dimensions[-1]
    assert example["d"] == 4096 and example["cost"] == 106496
    assert 1393 * 1393 - 2 * 985 * 985 == -1
    assert example["cost"]**985 < example["d"]**1393
    result = {
        "status": "PASS", "tensor_realization_asserted": False,
        "scope": "limitation of scalar cut-budget obstructions only",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "weighted_inequalities": sum(row["weighted_inequalities"] for row in dimensions),
        "small_cut_subsets": sum(row["small_cut_subsets"] for row in dimensions),
        "flagged_weighted_inequalities": sum(row["flagged_weighted_inequalities"] for row in dimensions),
        "flagged_subset_budget_classes": sum(row["flagged_subset_budget_classes"] for row in dimensions),
        "dimensions": dimensions,
        "last_example_strictly_below_target": True,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "dimensions"}, indent=2))


if __name__ == "__main__":
    main()
