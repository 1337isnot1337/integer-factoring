#!/usr/bin/env python3
"""Exact single-block test for recombinations before cube orientation.

The adaptive global-pairing step is rejected.  The favorable-reweighting
diagnostic excludes literal favorable orientations.  This checker asks a
different question in the same m=5, weight-1 Gaussian coefficient family:
can a recombination of the 22 pieces, followed by one fixed parent
pairing, restore E4 inside the Gaussian width regime?

There are only three fixed pairings.  The preimage of E4 under each is the
same recombination R of 11 pieces.  It has exact cut ranks (16,16,16),
widths (4,4,4), and width sum 12, above the 5/8 budget 6*rho*lambda'=15/2.
The endpoint can be restored only at effective rho=1.  This is a finite
single-block diagnostic, not an asymptotic or multi-block theorem.
"""
from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from check_cube_favorable_reweighting import span
from check_two_axis_global_pairing import (
    AXES, check_assignment_dependent_pairing, pullback,
)


def binary_rank(rows: list[int]) -> int:
    pivots = {}
    for value in rows:
        while value:
            lead = value.bit_length() - 1
            if lead in pivots:
                value ^= pivots[lead]
            else:
                pivots[lead] = value
                break
    return len(pivots)


def rank_at_axis(tensor: int, axis: int) -> int:
    left = [v for v in range(8) if not ((v >> (2 - axis)) & 1)]
    right = [v for v in range(8) if (v >> (2 - axis)) & 1]
    rows = [0] * 16
    for mask in range(256):
        if (tensor >> mask) & 1:
            row = sum(((mask >> v) & 1) << j for j, v in enumerate(left))
            column = sum(((mask >> v) & 1) << j for j, v in enumerate(right))
            rows[row] ^= 1 << column
    return binary_rank(rows)


def reorient(tensor: int, axes: tuple[int, int]) -> int:
    out = 0
    for newmask in range(256):
        oldmask = sum(bit << i for i, bit in enumerate(
            pullback(tuple((newmask >> v) & 1 for v in range(8)), axes)))
        if (tensor >> oldmask) & 1:
            out |= 1 << newmask
    return out


def inverse_reorient(tensor: int, axes: tuple[int, int]) -> int:
    out = 0
    for newmask in range(256):
        if (tensor >> newmask) & 1:
            oldmask = sum(bit << i for i, bit in enumerate(
                pullback(tuple((newmask >> v) & 1 for v in range(8)), axes)))
            out |= 1 << oldmask
    return out


def report() -> dict:
    source = check_assignment_dependent_pairing()
    endpoint = int(source["expected_truth_bits"], 16)
    pieces = sorted((tuple(piece["exponent"]), int(piece["original_truth_bits"], 16))
                    for piece in source["pieces"])
    tensors = [tensor for _, tensor in pieces]
    assert len(tensors) == 22

    fixed_pairings = {}
    for axes in AXES:
        total = 0
        for tensor in tensors:
            total ^= reorient(tensor, axes)
        assert total == endpoint
        fixed_pairings[str(axes)] = True

    preimages = []
    for hidden in range(3):
        axes = tuple(axis for axis in range(3) if axis != hidden)
        candidate = inverse_reorient(endpoint, axes)
        assert reorient(candidate, axes) == endpoint
        preimages.append(candidate)
    assert len(set(preimages)) == 1
    recombination = preimages[0]

    representation = span(tensors, recombination)
    assert representation["endpoint_in_span"]
    subset = representation["solution"]
    assert subset is not None
    total = 0
    for index in subset:
        total ^= tensors[index]
    assert total == recombination

    ranks = [rank_at_axis(recombination, axis) for axis in range(3)]
    widths = [(rank - 1).bit_length() if rank else 0 for rank in ranks]
    assert ranks == [16, 16, 16]
    assert widths == [4, 4, 4]
    width_sum = sum(widths)
    rho = F(5, 8)
    lam_prime = 2
    cap = 2 * lam_prime
    budget = 6 * rho * lam_prime
    assert width_sum > budget
    effective_rho = F(width_sum, 6 * lam_prime)
    assert effective_rho == 1
    assert all(width <= cap for width in widths)

    dual = (1 << 240) | (1 << 255)
    assert (endpoint & dual).bit_count() % 2 == 1
    assert (recombination & dual).bit_count() % 2 == 1

    return {
        "family": "m=5, weight 1, Q={00001,00010}",
        "coefficient_pieces": len(tensors),
        "fixed_pairing_positive_controls": fixed_pairings,
        "preimages_coincide": True,
        "recombination_truth_bits": hex(recombination),
        "recombination_piece_indices": subset,
        "recombination_piece_count": len(subset),
        "recombination_exponents": [list(pieces[index][0]) for index in subset],
        "exact_cut_ranks": ranks,
        "exact_cut_widths": widths,
        "width_sum": width_sum,
        "per_cut_cap": cap,
        "gaussian_width_budget": str(budget),
        "budget_exceeded": True,
        "effective_rho": str(effective_rho),
        "favorable_dual_on_recombination": (recombination & dual).bit_count() % 2,
        "favorable_dual_on_endpoint": (endpoint & dual).bit_count() % 2,
        "scope": ("Fixed 22-piece single-block Gaussian family only. Every "
                  "recombination restored to E4 by one fixed parent pairing "
                  "is this full-rank recombination, which leaves the 5/8 "
                  "width regime. Larger alphabets, other families, coupled "
                  "blocks and a new adaptive identity remain open."),
        "factoring_exponent_improved": False,
        "source_sha256": {
            name: sha256(Path("tools", name).read_bytes()).hexdigest()
            for name in (
                "check_recombination_before_orientation.py",
                "check_two_axis_global_pairing.py",
                "check_cube_favorable_reweighting.py",
                "cube_homogeneous_gaussian.py",
            )
        },
    }


def main() -> None:
    data = report()
    Path("evidence/recombination-before-orientation.json").write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "passed",
        "recombination_piece_count": data["recombination_piece_count"],
        "exact_cut_ranks": data["exact_cut_ranks"],
        "exact_cut_widths": data["exact_cut_widths"],
        "width_sum": data["width_sum"],
        "gaussian_width_budget": data["gaussian_width_budget"],
        "effective_rho": data["effective_rho"],
        "factoring_exponent_improved": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
