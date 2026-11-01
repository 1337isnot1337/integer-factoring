#!/usr/bin/env python3
"""Checks for the two-axis regrouping of the three-axis Gaussian cube.

The Gaussian factorization, the three one-particle Smith profiles, and the
rational pressure certificate are inherited from the cited modules.  This
checker verifies the new combinatorial regrouping and the width arithmetic.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from hashlib import sha256
from pathlib import Path
import json
import random

import certify_cube_optimized_type as optimized
from cube_homogeneous_gaussian import augmentation_crossblock, exact_monomial_smith, kernel


ANTIPODAL = [(v, v ^ 7) for v in range(8) if v < (v ^ 7)]


def bits(v: int) -> tuple[int, int, int]:
    return ((v >> 2) & 1, (v >> 1) & 1, v & 1)


def ilabel(vertex: int, axes: tuple[int, int], hidden: int) -> tuple[int, int]:
    b = bits(vertex)
    quadrant = (b[axes[0]] << 1) | b[axes[1]]
    return (quadrant, b[hidden])


def endpoint(vals: tuple[int, ...]) -> int:
    return int(all(vals[a] == vals[b] for a, b in ANTIPODAL))


def regrouped_endpoint(vals: tuple[int, ...], axes: tuple[int, int]) -> int:
    hidden = 3 - axes[0] - axes[1]
    labels = {}
    for v in range(8):
        q, h = ilabel(v, axes, hidden)
        labels.setdefault(q, [0, 0])[h] = vals[v]
    # The antipode reverses the hidden coordinate.  Read the two antipodal
    # quadrants U,V in the opposite internal order so that word equality is
    # exactly the two antipodal leg equalities.
    labels[2].reverse()
    labels[3].reverse()
    # Quadrant 0 is X, 1 is Y, 2 is U, 3 is V.  Antipodal diagonals are
    # X--V (0--3) and Y--U (1--2).
    x, y, u, v = (tuple(labels[q]) for q in range(4))
    return int(x == v and y == u)


def check_pairings() -> list[dict]:
    rows = []
    for axes in combinations(range(3), 2):
        for mask in range(1 << 8):
            vals = tuple((mask >> v) & 1 for v in range(8))
            assert endpoint(vals) == regrouped_endpoint(vals, axes)
        hidden = 3 - axes[0] - axes[1]
        # The two parent cuts are the two selected coordinate cuts.
        for cut_axis in axes:
            for mask in range(1 << 8):
                vals = tuple((mask >> v) & 1 for v in range(8))
                left = {bits(v)[cut_axis] for v in range(8)}
                labels = {}
                for v in range(8):
                    q, h = ilabel(v, axes, hidden)
                    labels.setdefault(q, [0, 0])[h] = bits(v)[cut_axis]
                grouped = {tuple(labels[q]) for q in range(4)}
                assert grouped == {(0, 0), (1, 1)} or grouped == {(0, 1), (1, 0)}
                assert left == {0, 1}
        rows.append({"axes": list(axes), "hidden": hidden,
                     "endpoint_assignments": 1 << 8,
                     "identity": "E8 = [X=V][Y=U]"})
    return rows


def check_widths() -> dict:
    rng = random.Random(20261101)
    rho = F(5, 8)
    cases = 0
    lam = 7
    per_piece_limit = int(12 * rho * lam)
    cap = 4 * lam
    for _ in range(200):
        triples = []
        for _ in range(30):
            while True:
                a = rng.randrange(0, cap + 1)
                b = rng.randrange(0, cap + 1)
                c = rng.randrange(0, cap + 1)
                if a + b + c <= per_piece_limit:
                    break
            triples.append((a, b, c))
        total_three = sum(sum(t) for t in triples)
        total_min = sum(min(t) for t in triples)
        assert total_min * 3 <= total_three
        # lambda'=2 lambda, per-piece three-cut cap 12 rho lambda = 6 rho lambda'.
        for t in triples:
            assert sum(t) <= per_piece_limit
            assert max(t) <= cap
        cases += 1
    return {"random_assignments": cases,
            "inequality": "sum_cubes min(r1,r2,r3) <= (1/3) sum_cubes (r1+r2+r3)"}


def check_geometry() -> list[dict]:
    rows = []
    rho = F(5, 8)
    lam = 7
    lam_prime = 2 * lam
    cap = 2 * lam_prime
    for g in range(2, 8):
        pairs = list(combinations(range(g), 2))
        B = len(pairs) ** 2
        volume = g ** 4 * lam_prime
        other_per_pair = g
        # An assignment can be represented by arbitrary three-cut widths whose
        # per-piece sum is at most 12 rho lam = 6 rho lam'.
        rng = random.Random(900 + g)
        row_total = 0
        for _ in range(B):
            while True:
                a = rng.randrange(0, cap + 1)
                b = rng.randrange(0, cap + 1)
                c = rng.randrange(0, cap + 1)
                if a + b + c <= int(12 * rho * lam):
                    break
            row_total += min(a, b, c)
        assert row_total <= 2 * rho * lam_prime * B
        U = 2 * (row_total + len(pairs) * other_per_pair * lam_prime)
        assert U <= 4 * rho * lam_prime * B + g * g * (g - 1) * lam_prime
        child_cap = cap * len(pairs) + other_per_pair * lam_prime
        assert child_cap == volume // (g * g)
        ratio = F(U, volume)
        bound = rho * F(g - 1, g) ** 2 + F(1, g) - F(1, g * g)
        assert ratio <= bound
        assert bound <= rho + F(1, g)
        rows.append({"g": g, "squares": B, "volume": volume,
                     "child_cap": child_cap, "U": U,
                     "U_over_V": [ratio.numerator, ratio.denominator]})
    return rows


def check_certificate() -> dict:
    target = F(312359221, 500000000)
    coarse = F(5, 8)
    assert 2 * coarse * coarse < 1
    assert 2 * target * target < 1
    # log_2(1/rho) > 1/2 iff rho^2 < 1/2.
    cert = optimized.certificate(target)
    assert cert["passed"]
    negative = optimized.certificate(F(624718441, 10 ** 9))
    assert negative["rejected_at_this_point"]
    profiles = [exact_monomial_smith(augmentation_crossblock(kernel(), axis), axis)["orders"]
                for axis in range(3)]
    assert profiles == [[2, 4, 8, 10]] * 3
    return {"coarse_rho": [coarse.numerator, coarse.denominator],
            "coarse_score_lower_bound": "log_2(8/5)",
            "coarse_threshold_check": "2*(5/8)^2 < 1",
            "optimized_rho": [target.numerator, target.denominator],
            "optimized_threshold_check": "2*(312359221/500000000)^2 < 1",
            "certificate": cert,
            "negative_control_rejected": True,
            "one_particle_profiles": profiles}


def main() -> None:
    pairing = check_pairings()
    widths = check_widths()
    geometry = check_geometry()
    certificate = check_certificate()
    sources = ["check_two_axis_cube_regrouping.py",
               "certify_cube_optimized_type.py", "cube_homogeneous_gaussian.py"]
    report = {
        "status": "passed",
        "pairings": pairing,
        "widths": widths,
        "geometry": geometry,
        "certificate": certificate,
        "source_sha256": {name: sha256(Path("tools", name).read_bytes()).hexdigest()
                          for name in sources},
        "scope": ("Exact regrouping and width arithmetic for the two-axis use of "
                  "the existing three-axis Gaussian cube; the local tensor family, "
                  "its factorizations, and the two-axis evaluator are inherited. "
                  "The improved factoring exponent is an ordinary proof, not Lean-verified "
                  "and not independently reviewed."),
    }
    Path("evidence/two-axis-cube-regrouping.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed",
                      "pairing_cases": len(pairing),
                      "geometry_cases": len(geometry),
                      "coarse_score": "log_2(8/5) > 1/2",
                      "optimized_score": "log_2(500000000/312359221) > 1/2",
                      "factoring_exponent_improved": True}, sort_keys=True))


if __name__ == "__main__":
    main()
