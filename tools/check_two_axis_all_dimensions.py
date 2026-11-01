#!/usr/bin/env python3
"""Exact checks for the dimension-independent two-axis cube regrouping.

The checker verifies endpoint regrouping, cut incidence, width arithmetic,
and the inherited three-axis pressure certificate. It does not construct a
four-or-higher-axis family and does not re-prove the two-axis transfer.
"""
from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import math

import certify_cube_optimized_type as optimized
from cube_homogeneous_gaussian import augmentation_crossblock, exact_monomial_smith, kernel


def endpoint(vals: tuple[int, ...], d: int) -> int:
    n = 1 << d
    flip = n - 1
    return int(all(vals[v] == vals[v ^ flip]
                   for v in range(n) if v < (v ^ flip)))


def regrouped_endpoint(vals: tuple[int, ...], d: int,
                       parents: tuple[int, int]) -> int:
    hidden = [a for a in range(d) if a not in parents]
    groups: dict[tuple[int, ...], dict[tuple[int, ...], tuple[int, ...]]] = {}
    for v in range(1 << d):
        key = tuple((v >> a) & 1 for a in parents)
        hkey = tuple((v >> a) & 1 for a in hidden)
        groups.setdefault(key, {})[hkey] = tuple(vals[v] for a in hidden)
    keys = sorted(groups)
    labels = [tuple(groups[key][h] for h in sorted(groups[key]))
              for key in keys]
    # The antipode reverses the hidden coordinates. Read the two antipodal
    # parent quadrants in the opposite internal order.
    labels[2] = labels[2][::-1]
    labels[3] = labels[3][::-1]
    x, y, u, v = labels
    return int(x == v and y == u)


def check_endpoints(d: int) -> dict:
    cases = 0
    pairs = list(combinations(range(d), 2))
    for mask in range(1 << (1 << d)):
        vals = tuple((mask >> v) & 1 for v in range(1 << d))
        target = endpoint(vals, d)
        for parents in pairs:
            assert regrouped_endpoint(vals, d, parents) == target
            cases += 1
    return {"dimension": d, "pairings": len(pairs),
            "assignments": 1 << (1 << d), "endpoint_cases": cases}


def check_incidence(d: int) -> dict:
    rows = []
    for parents in combinations(range(d), 2):
        hidden = [a for a in range(d) if a not in parents]
        for cut in range(d):
            split = 0
            for key in range(4):
                sides = set()
                for hkey in range(1 << (d - 2)):
                    v = 0
                    for j, a in enumerate(parents):
                        if (key >> (1 - j)) & 1:
                            v |= 1 << a
                    for j, a in enumerate(hidden):
                        if (hkey >> j) & 1:
                            v |= 1 << a
                    sides.add((v >> cut) & 1)
                if len(sides) == 2:
                    split += 1
            if cut in parents:
                assert split == 0
            else:
                assert split == 4
            rows.append({"parents": list(parents), "cut": cut,
                         "superlegs_split": split})
    return {"dimension": d, "rows": rows}


def check_widths(rho: F) -> list[dict]:
    rows = []
    for d in range(3, 10):
        total = d * (1 << (d - 1)) * rho
        parent_pair = F(2, d) * total
        lam_prime = 1 << (d - 2)
        assert parent_pair == 4 * rho * lam_prime
        threshold = 4 * rho * lam_prime
        rows.append({
            "d": d,
            "lambda_prime": lam_prime,
            "per_cube_total": [total.numerator, total.denominator],
            "parent_pair": [parent_pair.numerator, parent_pair.denominator],
            "four_rho_lambda_prime": [threshold.numerator, threshold.denominator],
            "separate_width_cap": 2 * lam_prime,
        })
    return rows


def profile(p: float, r: float) -> float:
    Z = (1 + r) * (1 + r * r) * (1 + r ** 4) * (1 + r ** 5)
    H = -(p * math.log(p) + (1 - p) * math.log(1 - p))
    return (math.log(Z) - 4 * p * math.log(r)) / (4 * H)


def profile_grid() -> dict:
    best = (float("inf"), None)
    for i in range(1, 400):
        p = i / 400
        for j in range(1, 400):
            r = j / 400
            value = profile(p, r)
            if value < best[0]:
                best = (value, (p, r))
    return {"grid": 400, "minimum": best[0],
            "at": [best[1][0], best[1][1]],
            "scope": "bounded numerical diagnostic, not a global optimum proof"}


def check_certificate() -> dict:
    target = F(312359221, 500000000)
    coarse = F(5, 8)
    assert 2 * coarse * coarse < 1
    assert 2 * target * target < 1
    cert = optimized.certificate(target)
    assert cert["passed"]
    negative = optimized.certificate(F(624718441, 10 ** 9))
    assert negative["rejected_at_this_point"]
    profiles = [exact_monomial_smith(augmentation_crossblock(kernel(), axis), axis)["orders"]
                for axis in range(3)]
    assert profiles == [[2, 4, 8, 10]] * 3
    return {"coarse_rho": [coarse.numerator, coarse.denominator],
            "optimized_rho": [target.numerator, target.denominator],
            "certificate": cert,
            "negative_control_rejected": True,
            "one_particle_profiles": profiles,
            "alpha_lower_bound_diagnostic": -math.log2(float(target))}


def main() -> None:
    endpoint_rows = [check_endpoints(3), check_endpoints(4)]
    incidence_rows = [check_incidence(3), check_incidence(4)]
    width_rows = check_widths(F(5, 8))
    certificate = check_certificate()
    diagnostic = profile_grid()
    sources = ["check_two_axis_all_dimensions.py",
               "certify_cube_optimized_type.py",
               "cube_homogeneous_gaussian.py"]
    report = {
        "status": "passed",
        "endpoint_regrouping": endpoint_rows,
        "cut_incidence": incidence_rows,
        "width_arithmetic": width_rows,
        "three_axis_certificate": certificate,
        "profile_diagnostic": diagnostic,
        "source_sha256": {name: sha256(Path("tools", name).read_bytes()).hexdigest()
                          for name in sources},
        "scope": ("Exact endpoint regrouping and width arithmetic for a global "
                  "two-axis use of a d-axis cube family, d>=3. No d>=4 family, "
                  "no new factoring exponent, and no re-proof of the two-axis "
                  "transfer are supplied."),
        "factoring_exponent_improved": False,
    }
    Path("evidence/two-axis-all-dimensions.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed",
                      "endpoint_dimensions": [3, 4],
                      "width_dimensions": list(range(3, 10)),
                      "profile_grid_minimum": diagnostic["minimum"],
                      "factoring_exponent_improved": False}, sort_keys=True))


if __name__ == "__main__":
    main()
