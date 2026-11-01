#!/usr/bin/env python3
"""Exact checks for the global-pairing replacement of the rejected two-axis
regrouping.

The audit in check_two_axis_cube_regrouping.py rejects the original per-cube
minimum step.  This checker verifies a different interface: one common
internal coordinate for the whole piece assignment, an exact identity on the
fixed parent labels, and the total parent-width bound used by the multiaxis
transfer.  The inherited Gaussian pieces, profiles and certificate are
replayed from the audit checker.
"""
from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import random

import check_two_axis_cube_regrouping as audit


AXES = [(0, 1), (0, 2), (1, 2)]


def old_vertex(q: int, h: int, axes: tuple[int, int]) -> int:
    """Physical vertex represented by parent label (q,h) for one pairing."""
    a, b = axes
    hidden = 3 - a - b
    h_old = h if q < 2 else 1 - h
    coords = [0, 0, 0]
    coords[a] = (q >> 1) & 1
    coords[b] = q & 1
    coords[hidden] = h_old
    return (coords[0] << 2) | (coords[1] << 1) | coords[2]


def pullback(newvals: tuple[int, ...], axes: tuple[int, int]) -> tuple[int, ...]:
    """Physical eight-bit assignment induced by a four-label assignment."""
    vals = [0] * 8
    for q in range(4):
        for h in range(2):
            vals[old_vertex(q, h, axes)] = newvals[(q << 1) | h]
    return tuple(vals)


def common_relabel(vals: tuple[int, ...], axes: tuple[int, int]) -> tuple[int, ...]:
    """Apply the common internal-coordinate map to an eight-bit tensor."""
    out = [0] * 8
    for q in range(4):
        for h in range(2):
            out[(q << 1) | h] = vals[old_vertex(q, h, axes)]
    return tuple(out)


def parent_endpoint(newvals: tuple[int, ...]) -> int:
    """E4(x,y,u,v)=[x=v][y=u] on four two-bit labels."""
    return int(
        all(newvals[h] == newvals[6 + h] for h in range(2))
        and all(newvals[2 + h] == newvals[4 + h] for h in range(2))
    )


def relabel_scalar(fn, axes: tuple[int, int]) -> tuple[int, ...]:
    return tuple(fn(pullback(tuple((mask >> v) & 1 for v in range(8)), axes))
                 for mask in range(256))


def check_common_map() -> list[dict]:
    rows = []
    for axes in AXES:
        for mask in range(256):
            vals = tuple((mask >> v) & 1 for v in range(8))
            relabeled = common_relabel(vals, axes)
            assert audit.regrouped_endpoint(vals, axes) == parent_endpoint(relabeled)
        # Common-map positive control: the audit's T0=E8+A, T1=A pair.
        def e8(vals):
            return audit.endpoint(vals)

        def t0(vals):
            return audit.endpoint(vals) ^ vals[1]

        def t1(vals):
            return vals[1]

        t0r = relabel_scalar(t0, axes)
        t1r = relabel_scalar(t1, axes)
        e4r = relabel_scalar(e8, axes)
        assert all((t0r[m] ^ t1r[m]) == e4r[m] for m in range(256))
        rows.append({"axes": list(axes), "assignments_checked": 256,
                     "common_map_identity": "P_c(E8)=E8 and E8=sum_i P_c(T_i)"})
    return rows


def check_piecewise_failure() -> dict:
    control = audit.check_independent_relabeling()
    assert control["assignments"] == 256
    assert control["mismatches"] == 128
    return control


def check_global_widths() -> dict:
    rho = F(5, 8)
    lam = 7
    lam_prime = 2 * lam
    cap = 2 * lam_prime
    per_piece_limit = int(6 * rho * lam_prime)
    # The audit's two-width counterexample.  The global choice makes the
    # hidden coordinate the largest total and uses the other two as parents.
    triples = [(0, 3, 4), (3, 0, 4)]
    totals = [sum(t[j] for t in triples) for j in range(3)]
    internal = max(range(3), key=lambda j: totals[j])
    parents = [totals[j] for j in range(3) if j != internal]
    selected = min(parents)
    global_min = min(totals)
    independent_min = sum(min(t) for t in triples)
    assert global_min == 3 and independent_min == 0
    assert selected == global_min
    assert sum(parents) * 3 <= 2 * sum(totals)
    example_selected = selected

    rng = random.Random(20261101)
    cases = 0
    for _ in range(200):
        rows = []
        for _ in range(30):
            while True:
                a = rng.randrange(0, cap + 1)
                b = rng.randrange(0, cap + 1)
                c = rng.randrange(0, cap + 1)
                if a + b + c <= per_piece_limit:
                    break
            rows.append((a, b, c))
        total_three = sum(sum(t) for t in rows)
        totals = [sum(t[j] for t in rows) for j in range(3)]
        internal = max(range(3), key=lambda j: totals[j])
        parents = [totals[j] for j in range(3) if j != internal]
        pair_total = sum(parents)
        selected = min(parents)
        assert pair_total * 3 <= 2 * total_three
        assert pair_total <= 4 * rho * lam_prime * len(rows)
        assert selected <= 2 * rho * lam_prime * len(rows)
        for t in rows:
            assert sum(t) <= per_piece_limit
            assert max(t) <= cap
        cases += 1
    return {"audit_counterexample": [list(t) for t in triples],
            "independent_minimum_sum": independent_min,
            "selected_global_parent_width": example_selected,
            "established_global_minimum": global_min,
            "random_assignments": cases,
            "interface": "common internal coordinate; pair_total <= "
                         "(2/3) sum_i S_i <= 4 rho lambda' B"}


def check_hidden_incidence() -> dict:
    control = audit.check_fixed_bundle_cuts()
    cuts = control["cuts"]
    assert [row["cut_on_whole_parent_labels"] for row in cuts] == [True, True, False]
    assert cuts[2]["split_superlegs"] == [0, 1, 2, 3]
    return control


def check_total_width_rank() -> dict:
    lam_prime = 32
    pieces = 4
    rho = F(5, 8)
    B = (lam_prime + 1) ** pieces
    endpoint_rank = 1 << (2 * lam_prime)
    uniform_width = int(2 * rho * lam_prime)
    uniform_upper = B * (1 << uniform_width)
    assert uniform_upper < endpoint_rank
    cap = 2 * lam_prime
    total_width_budget = int(2 * rho * lam_prime * B)
    assert cap <= total_width_budget
    rank_sum = (1 << cap) + (B - 1)
    assert rank_sum >= endpoint_rank
    return {"lambda_prime": lam_prime,
            "piece_exponent": pieces,
            "pieces": B,
            "endpoint_rank": endpoint_rank,
            "uniform_piece_width": uniform_width,
            "uniform_rank_upper": uniform_upper,
            "uniform_bound_below_endpoint": uniform_upper < endpoint_rank,
            "cap_width": cap,
            "total_width_budget": total_width_budget,
            "rank_with_one_cap_piece": rank_sum,
            "conclusion": "rank excludes the uniform per-piece fixed-cut "
                          "bound, not the total-width interface"}


def check_geometry() -> list[dict]:
    rows = []
    rho = F(5, 8)
    lam = 7
    lam_prime = 2 * lam
    cap = 2 * lam_prime
    per_piece_limit = int(6 * rho * lam_prime)
    for g in range(2, 8):
        pairs = list(range(g * (g - 1) // 2))
        B = len(pairs) ** 2
        volume = g ** 4 * lam_prime
        rng = random.Random(900 + g)
        triples = []
        for _ in range(B):
            while True:
                a = rng.randrange(0, cap + 1)
                b = rng.randrange(0, cap + 1)
                c = rng.randrange(0, cap + 1)
                if a + b + c <= per_piece_limit:
                    break
            triples.append((a, b, c))
        totals = [sum(t[j] for t in triples) for j in range(3)]
        internal = max(range(3), key=lambda j: totals[j])
        parents = [totals[j] for j in range(3) if j != internal]
        selected = min(parents)
        assert sum(parents) <= 4 * rho * lam_prime * B
        assert selected <= 2 * rho * lam_prime * B
        U = 2 * (selected + len(pairs) * g * lam_prime)
        assert U <= 4 * rho * lam_prime * B + g * g * (g - 1) * lam_prime
        child_cap = cap * len(pairs) + g * lam_prime
        assert child_cap == volume // (g * g)
        ratio = F(U, volume)
        bound = rho * F(g - 1, g) ** 2 + F(1, g) - F(1, g * g)
        assert ratio <= bound
        assert bound <= rho + F(1, g)
        rows.append({"g": g, "squares": B, "volume": volume,
                     "internal": internal, "selected_parent_width": selected,
                     "child_cap": child_cap, "U": U,
                     "U_over_V": [ratio.numerator, ratio.denominator]})
    return rows


def main() -> None:
    common = check_common_map()
    piecewise = check_piecewise_failure()
    widths = check_global_widths()
    hidden = check_hidden_incidence()
    rank = check_total_width_rank()
    geometry = check_geometry()
    certificate = audit.check_certificate()
    sources = ["check_two_axis_global_pairing.py",
               "check_two_axis_cube_regrouping.py",
               "certify_cube_optimized_type.py",
               "cube_homogeneous_gaussian.py"]
    report = {
        "status": "replacement_interface_checked",
        "scoreboard": "remains 1/2 until independent acceptance",
        "common_map": common,
        "piecewise_counterexample_preserved": piecewise,
        "global_widths": widths,
        "hidden_incidence": hidden,
        "total_width_rank": rank,
        "geometry": geometry,
        "certificate": certificate,
        "source_sha256": {name: sha256(Path("tools", name).read_bytes()).hexdigest()
                          for name in sources},
        "scope": ("Global-pairing interface replacement for the rejected per-cube "
                  "two-axis regrouping.  The common internal-coordinate map fixes "
                  "the parent labels for a whole assignment; the width hypothesis "
                  "is the total parent-axis bound, not a uniform per-piece fixed-cut "
                  "bound.  The inherited Gaussian pieces, profiles, certificate and "
                  "two-axis evaluator are not re-proved.  The claim is an ordinary "
                  "proof, not Lean-verified and not independently accepted."),
    }
    Path("evidence/two-axis-global-pairing-repair.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed",
                      "common_map_pairings": len(common),
                      "piecewise_mismatches": piecewise["mismatches"],
                      "geometry_cases": len(geometry),
                      "selected_global_parent_width": widths["selected_global_parent_width"],
                      "rank_feasible": rank["rank_with_one_cap_piece"] >= rank["endpoint_rank"],
                      "coarse_score": "log_2(8/5) > 1/2",
                      "optimized_score": "log_2(500000000/312359221) > 1/2",
                      "scoreboard": "1/2 until independent acceptance"}, sort_keys=True))


if __name__ == "__main__":
    main()
