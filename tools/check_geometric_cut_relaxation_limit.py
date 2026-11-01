#!/usr/bin/env python3
"""Exact dyadic profile controls for geometric split constraints."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from check_weighted_cut_relaxation_limit import profile


def envelope(d: int, m: int) -> int:
    if m == 0:
        return d * d
    if m <= 2:
        return d * (d - 1) // 2
    if m >= 2 * d - 1:
        return 0
    return d * d // (4 * (m - 2))


def geometric(d: int, m: int) -> int:
    """All exact integer/half-integer split allocations, by doubled a."""
    result = 0
    for twice_a in range(2, 2 * m - 1):
        twice_b = 2 * m - twice_a
        # Odd allocations require a,b >= 3/2; the loop excludes 1/2.
        assert twice_a >= 2 and twice_b >= 2
        for s in range(1, d):
            q = min((2 * d - 2 * s) // twice_a,
                    (2 * d - twice_b * s) // 2)
            if q > 0:
                assert 2 * s + twice_a * q <= 2 * d
                assert twice_b * s + 2 * q <= 2 * d
                result = max(result, s * q)
    return result


def check_split_envelope(d: int) -> int:
    bounds = [geometric(d, m) for m in range(2 * d + 2)]
    padded = [max(bounds[m:], default=0) for m in range(len(bounds))]
    assert all(bounds[m] == 0 for m in range(2 * d - 1, len(bounds)))
    for m in range(3, len(bounds)):
        assert bounds[m] <= padded[m] <= envelope(d, m)
    for m in (1, 2):
        assert padded[m] <= envelope(d, m)
    return len(bounds)


def check_dimension(d: int) -> dict[str, int]:
    pairs = profile(d)
    cap = 2 * d - 1
    bounds = [envelope(d, m) for m in range(cap + 1)]
    checked = 0
    for exchanged in (False, True):
        oriented = [(c, r) if exchanged else (r, c) for r, c in pairs]
        total = sum(r for r, _ in oriented)
        removed = [-1] * (cap + 1)
        removed[0] = 0
        for r, c in oriented:
            for m in range(cap, c - 1, -1):
                if removed[m - c] >= 0:
                    removed[m] = max(removed[m], removed[m - c] + r)
        assert min(removed) >= 0
        for m, maximum in enumerate(removed):
            assert total - maximum >= bounds[m], (d, m, total - maximum, bounds[m])
            checked += 1
    for s in range(d + 1):
        minimum = min(value + s * s * m for m, value in enumerate(bounds))
        assert sum(min(r, s * s * c) for r, c in pairs) >= minimum
        assert sum(min(s * s * r, c) for r, c in pairs) >= minimum
    return {"d": d, "pieces": len(pairs), "cost": len(pairs) * d,
            "subset_budget_classes": checked,
            "weighted_inequalities": 2 * (d + 1)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    split_checks = sum(check_split_envelope(1 << m) for m in range(1, 7))
    dimensions = [check_dimension(1 << m) for m in range(1, 13)]
    assert dimensions[-1]["cost"] == 106496
    assert 1393 * 1393 - 2 * 985 * 985 == -1
    assert 106496**985 < 4096**1393
    result = {
        "status": "PASS", "tensor_realization_asserted": False,
        "scope": "scalar limitation including integer and half-integer geometric splits",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "all_split_envelope_checks": split_checks,
        "subset_budget_classes": sum(x["subset_budget_classes"] for x in dimensions),
        "weighted_inequalities": sum(x["weighted_inequalities"] for x in dimensions),
        "dimensions": dimensions,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "dimensions"}, indent=2))


if __name__ == "__main__":
    main()
