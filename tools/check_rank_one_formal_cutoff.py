#!/usr/bin/env python3
"""Exact factor-kernel refinement for formal crossing limits through 71."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from check_geometric_subspace_compression import geometric_complement
from check_geometric_border_cutoff import harmonic_bound, moving_case
from check_weighted_cut_rank_budget import (
    encoded_fraction, exact_above_sqrt2_power, exact_below_sqrt2_power, weights,
)
from check_weighted_cut_residual_budget import chord_coefficients


def kernel_complement(d, h):
    assert 1 <= h <= d
    return (2 * d - h) ** 2 // 4


@lru_cache(None)
def class_complement(d, h, m):
    assert 1 <= h <= d and m >= 0
    if m == 0:
        return d * d
    if m == 1:
        return kernel_complement(d, h)
    return geometric_complement(d, m)


@lru_cache(None)
def budget_without_one(d, s, minimum_products):
    end = max(minimum_products, 2 * d - 1)
    return min((d * d if m == 0 else geometric_complement(d, m)) + s * s * m
               for m in range(minimum_products, end + 1) if m != 1)


def class_budget(d, h, s, minimum_products=0):
    result = budget_without_one(d, s, minimum_products)
    if minimum_products <= 1:
        result = min(result, kernel_complement(d, h) + s * s)
    return result


def intervals(d, k, b):
    low, high = (k - 1) ** 2 * b + 1, min(k * k * b, d * d)
    if b == 1:
        assert high == k * k
        return ((low, high - 1, k - 1, "below_square"),
                (high, high, k, "at_square"))
    return ((low, high, k, "general"),)


def class_parameters(d, k, b, low, high, h):
    n = d * d
    assert (k - 1) ** 2 * b < low <= high <= min(k * k * b, n)
    assert h in (k - 1, k)
    f = kernel_complement(d, h) if b == 1 else geometric_complement(d, b)
    return {"selected_lower": low, "selected_upper": high, "factor_rank_bound": h,
            "complement": f, "row_rest": max(n - high, f),
            "column_rest": n - b, "row_surplus": max(0, low + f - n)}


def class_bound(d, k, b, low, high, h):
    p = class_parameters(d, k, b, low, high, h)
    r, c = p["row_rest"], p["column_rest"]
    result = Fraction(2 * d * d + p["row_surplus"], 2 * k - 1)
    for s, weight in enumerate(weights(k), 1):
        mu, nu = chord_coefficients(k, s)
        left = max(Fraction(class_budget(d, h, s, b)), s * s * b + mu * r + nu * c)
        right = max(Fraction(class_budget(d, h, s)), b + nu * r + mu * c)
        result += weight * (left + right)
    return result


def dimension_certificate(d):
    best, selected, count, below, square = Fraction(d * d), {"k": 1}, 1, 0, 0
    assert exact_above_sqrt2_power(best, d)
    minima = []
    for k in range(2, d + 1):
        local, local_parameters = None, None
        for b in range(1, (d * d - 1) // ((k - 1) ** 2) + 1):
            parts = intervals(d, k, b)
            assert parts[0][0] == (k - 1) ** 2 * b + 1
            assert parts[-1][1] == min(k * k * b, d * d)
            for i in range(1, len(parts)):
                assert parts[i - 1][1] + 1 == parts[i][0]
            for low, high, h, label in parts:
                value = class_bound(d, k, b, low, high, h)
                assert exact_above_sqrt2_power(value, d), (d, k, b, label)
                count += 1
                below += label == "below_square"
                square += label == "at_square"
                parameters = {"k": k, "b": b, "subcase": label,
                              **class_parameters(d, k, b, low, high, h)}
                if local is None or value < local:
                    local, local_parameters = value, parameters
                if value < best:
                    best, selected = value, parameters
        minima.append({"k": k, "bound": encoded_fraction(local), "parameters": local_parameters})
    assert below == square == d - 1
    return {"d": d, "bound": encoded_fraction(best), "above_target": True,
            "scalar_subclasses": count, "rank_one_below_square": below,
            "rank_one_at_square": square, "minimizing_class": selected,
            "minimum_by_k": minima}


def independent_controls():
    maxima = closures = factor_bounds = 0
    for d in range(2, 81):
        for h in range(1, d + 1):
            assert kernel_complement(d, h) == max(
                s * (d - max(0, s - d + h)) for s in range(d + 1))
            maxima += 1
            for p in range(1, h + 1):
                for q in (p, d):
                    assert min(p, q) ** 2 <= p * q
                    assert kernel_complement(d, p) >= kernel_complement(d, h)
                    factor_bounds += 1
        for h in sorted({1, d // 2, d - 1, d}):
            for s in range(d + 1):
                for minimum in (0, 1, 2, d, 2 * d - 1, 2 * d + 1):
                    direct = min(class_complement(d, h, m) + s * s * m
                                 for m in range(minimum, max(minimum, 4 * d) + 1))
                    assert class_budget(d, h, s, minimum) == direct
                    closures += 1
    moving = 0
    for d in range(4, 7):
        for p in range(2, d + 1):
            s = d - p // 2
            ell = s - d + p
            q = d - ell
            Ablock = [[int(i == j) for j in range(s)] for i in range(ell)]
            Dblock = [[int(i == j) for j in range(d - s)] for i in range(q)]
            assert ell + d - s == p
            assert s * q == kernel_complement(d, p)
            for prime in (2, 3, 5):
                moving_case(d, s, q, [Ablock], [Dblock], 3, prime)
                moving += 1
    return {"kernel_integer_maxima": maxima, "factor_rank_controls": factor_bounds,
            "class_minimum_closures": closures, "moving_rank_deficient_cases": moving}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 577 ** 2 - 2 * 408 ** 2 == 1
    assert 1393 ** 2 - 2 * 985 ** 2 == -1
    controls = independent_controls()
    for d in range(2, 65):
        assert exact_above_sqrt2_power(harmonic_bound(d), d)
    dimensions = [dimension_certificate(d) for d in range(65, 72)]
    control = class_bound(72, 72, 1, 72 * 72, 72 * 72, 72)
    assert exact_below_sqrt2_power(control, 72)
    names = ("check_geometric_subspace_compression", "check_geometric_border_cutoff",
             "check_weighted_cut_rank_budget", "check_weighted_cut_residual_budget")
    modules = [Path(__file__).resolve()] + [Path(__import__(name).__file__) for name in names]
    result = {"status": "PASS", "scope": "Supplied formal crossing limits; generic Laurent cut ranks",
              "previous_cutoff": 64, "last_certified_dimension": 71,
              "controls": controls, "dimensions": dimensions,
              "scalar_subclasses": sum(row["scalar_subclasses"] for row in dimensions),
              "dimension_72_control": {"k": 72, "b": 1, "subcase": "at_square",
                  "bound": encoded_fraction(control), "below_target": True,
                  "meaning": "This certificate is insufficient; no realizable seed is asserted."},
              "exact_fitting_or_half_imported": False,
              "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in modules}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'dimensions'}, indent=2))
    for row in dimensions:
        print(json.dumps({key: value for key, value in row.items() if key != 'minimum_by_k'}))


if __name__ == '__main__':
    main()
