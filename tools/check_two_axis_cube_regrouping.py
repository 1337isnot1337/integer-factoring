#!/usr/bin/env python3
"""Audit the rejected two-axis regrouping of the Gaussian cube.

The Gaussian factorization, the three one-particle Smith profiles, and the
rational pressure certificate are inherited from the cited modules.  This
checker retains the valid endpoint/arithmetic controls and adds exact
counterexamples to independent summand relabeling and the fixed-bundle cut.
None of these controls establishes the rejected one-cut transfer.
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


def check_independent_relabeling() -> dict:
    """A common endpoint symmetry is not an independent summand symmetry."""
    def transpose(v):
        return ((v & 1) << 2) | (v & 2) | ((v >> 2) & 1)
    assert all(transpose(v ^ 7) == (transpose(v) ^ 7) for v in range(8))
    failures=[]
    for mask in range(256):
        vals=tuple((mask >> v) & 1 for v in range(8))
        permuted=tuple(vals[transpose(v)] for v in range(8))
        original=endpoint(vals)
        assert endpoint(permuted)==original
        # T0=E8+A, T1=A, A=x_001. Their original sum is exactly E8.
        t0=original ^ vals[1]
        t1=vals[1]
        assert (t0 ^ t1)==original
        rotated_t1=permuted[1]
        changed=t0 ^ rotated_t1
        if changed!=original:
            failures.append(dict(mask=mask,endpoint=original,relabelled_sum=changed))
    assert len(failures)==128
    return dict(assignments=256,mismatches=len(failures),first=failures[0],
                decomposition='T0=E8+x_001, T1=x_001 over F2',
                relabeling='Swap cube coordinates 0 and 2 in T1 only',
                conclusion='Independent summand relabeling does not preserve the endpoint')


def check_fixed_bundle_cuts() -> dict:
    axes=(0,1);hidden=2
    labels={q:[] for q in range(4)}
    for v in range(8):
        q,h=ilabel(v,axes,hidden)
        labels[q].append(v)
    rows=[]
    for axis in range(3):
        split=[q for q,vertices in labels.items()
               if len({bits(v)[axis] for v in vertices})>1]
        rows.append(dict(axis=axis,split_superlegs=split,
                         cut_on_whole_parent_labels=not split))
    assert rows[0]['split_superlegs']==rows[1]['split_superlegs']==[]
    assert rows[2]['split_superlegs']==[0,1,2,3]
    triples=[(0,3,4),(3,0,4)]
    local=sum(min(x) for x in triples)
    global_best=min(sum(x[a] for x in triples) for a in range(3))
    assert local==0 and global_best==3
    return dict(pairing_axes=axes,hidden_axis=hidden,cuts=rows,
                width_triples=triples,sum_independent_minima=local,
                minimum_global_axis_sum=global_best,
                conclusion='Hidden-axis cut splits all labels; independent minima are not one global cut')


def check_fixed_cut_rank() -> dict:
    # E4 is a permutation matrix at the fixed row flattening, for every
    # alphabet size. Here Q={0,1}^2 gives rank |Q|^2=16.
    size=4
    rows=[]
    for x in range(size):
        for y in range(size):
            row=sum(int(x==v and y==u) << (size*u+v)
                    for u in range(size) for v in range(size))
            rows.append(row)
    assert len(set(rows))==16 and all(r and r&(r-1)==0 for r in rows)
    lam_prime=32;exponent=4;rho=F(5,8)
    piece_bound=(lam_prime+1)**exponent
    channel_width=int(2*rho*lam_prime)
    upper=piece_bound*(1 << channel_width)
    endpoint_rank=1 << (2*lam_prime)
    assert upper<endpoint_rank
    return dict(finite_endpoint_rank=16,lambda_prime=lam_prime,
                illustrative_piece_exponent=exponent,piece_bound=piece_bound,
                claimed_channel_width=channel_width,summed_rank_upper=upper,
                endpoint_rank=endpoint_rank,
                conclusion='A polynomial-piece fixed-cut decomposition cannot have uniform width 2*rho*lambda+O(log lambda) for rho<1')


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
    relabeling = check_independent_relabeling()
    bundle_cuts = check_fixed_bundle_cuts()
    fixed_cut_rank = check_fixed_cut_rank()
    sources = ["check_two_axis_cube_regrouping.py",
               "certify_cube_optimized_type.py", "cube_homogeneous_gaussian.py"]
    report = {
        "status": "audit_passed_proposed_transfer_rejected",
        "pairings": pairing,
        "widths": widths,
        "geometry": geometry,
        "certificate": certificate,
        "independent_relabeling": relabeling,
        "fixed_bundle_cuts": bundle_cuts,
        "fixed_cut_rank": fixed_cut_rank,
        "factoring_exponent_improved": False,
        "source_sha256": {name: sha256(Path("tools", name).read_bytes()).hexdigest()
                          for name in sources},
        "scope": ("Endpoint pairings and formal width inequalities survive, but "
                  "the one-cut transfer is rejected: independently relabeling pieces "
                  "does not preserve their sum, and the hidden cut splits parent labels. "
                  "Uniform submaximal fixed-cut widths also contradict endpoint rank. "
                  "No improved complete factoring exponent follows."),
    }
    Path("evidence/two-axis-cube-regrouping.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report['status'],
                      "pairing_cases": len(pairing),
                      "geometry_cases": len(geometry),
                      "independent_relabeling_mismatches": relabeling['mismatches'],
                      "factoring_exponent_improved": False}, sort_keys=True))


if __name__ == "__main__":
    main()
