#!/usr/bin/env python3
"""Exact partition counts supporting the second-order channel calculation."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path


def partitions(maximum):
    values = [1] + [0] * maximum
    shifts = []
    k = 1
    while k * (3 * k - 1) // 2 <= maximum:
        sign = 1 if k % 2 else -1
        shifts.append((k * (3 * k - 1) // 2, sign))
        shifts.append((k * (3 * k + 1) // 2, sign))
        k += 1
    shifts.sort()
    for n in range(1, maximum + 1):
        total = 0
        for shift, sign in shifts:
            if shift > n:
                break
            total += sign * values[n - shift]
        values[n] = total
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-energy", type=int, default=65536)
    args = parser.parse_args()
    p = partitions(args.maximum_energy)
    small = [1] + [0] * 80
    for part in range(1, 81):
        for n in range(part, 81):
            small[n] += small[n - part]
    assert p[:81] == small
    cumulative = []
    total = 0
    for value in p:
        total += value
        cumulative.append(total)
    results = []
    for width in [8, 16, 32, 64, 128, 256, 512]:
        target = 1 << width
        h = next(n for n in range(2, len(p) // 2) if p[n] >= target)
        assert p[h - 1] < target <= p[h]
        r = cumulative[2 * h]
        exact_width = (r - 1).bit_length()
        assert (1 << (exact_width - 1)) < r <= (1 << exact_width)
        results.append({"lambda": width, "h": h, "partition_count": p[h],
            "central_channel_count": r, "central_channel_width": exact_width,
            "second_order_residual_estimate": math.log2(r) - math.sqrt(2) * width - (2 * math.sqrt(2) - 1) * math.log2(width),
            "top_energy_residual_estimate": math.log2(p[2 * h]) - math.sqrt(2) * width - (2 * math.sqrt(2) - 2) * math.log2(width)})
    result = {"status": "passed", "maximum_energy": args.maximum_energy,
        "independent_partition_dp_comparisons": 81, "profiles": results,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "Counts, minimal energies, and bit widths are exact. Residuals are numerical diagnostics; the asymptotic theorem has a separate proof."}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({"status": "passed", "profiles": [{k:a[k] for k in ("lambda", "h", "central_channel_width", "second_order_residual_estimate")} for a in results]}, indent=2))


if __name__ == "__main__":
    main()
