#!/usr/bin/env python3
"""Exact graded-projector ranks, boundary profiles and unfiltered span controls."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from check_quantum_braid_fusion_seam import apply_e, cut_ranks, rank, routes
from paired_path_crossing import power
from quantum_braid_fusion_seam import FIELD, INVERSE_DELTA, Q
from quantum_braid_graded_seam import (
    coefficient_tensor,
    free_counts,
    fusion_rank,
    projector_sites,
    rank_bounds,
)
from quantum_braid_path import EMBED, local_action

ROOT = Path(__file__).resolve().parents[1]


def reference_fusion(length, left, right):
    # Fundamental fusion follows A4; fusion by label3 reverses its vertices.
    objects = [1] * left + [3] * ((length - left) // 3)
    objects += [1] * right + [3] * ((length - right) // 3)
    state = [1, 0, 0, 0]
    for label in objects:
        if label == 3:
            state = state[::-1]
        else:
            state = [
                sum(state[j] for j in (h - 1, h + 1) if 0 <= j <= 3) for h in range(4)
            ]
    return state[0]


def right_e(matrix, route):
    diagonal, source, target = route
    result = FIELD[matrix, diagonal[None, :]]
    result[:, source] ^= FIELD[matrix[:, target], diagonal[source][None, :]]
    return result


def projector_matrix(length, left, right, rs):
    matrix = np.eye(len(paths(2 * length)), dtype=np.uint8)
    for site in projector_sites(length, left, right):
        a, b = apply_e(matrix, rs[site]), apply_e(matrix, rs[site + 1])
        matrix = (
            matrix
            ^ a
            ^ b
            ^ FIELD[
                INVERSE_DELTA,
                apply_e(b, rs[site]) ^ apply_e(a, rs[site + 1]),
            ]
        )
    return matrix


def direct_entry(four_paths, left, right, value):
    field = Field(0b10011)
    length = len(four_paths[0])
    states = {four_paths[0] + four_paths[1]: 1}

    def e(vector, site):
        result = {}
        for p, v in vector.items():
            for q, c in local_action(p, site):
                result[q] = result.get(q, 0) ^ field.mul(v, EMBED[c])
        return {p: v for p, v in result.items() if v}

    def combine(terms):
        result = {}
        for vector, scalar in terms:
            for p, v in vector.items():
                result[p] = result.get(p, 0) ^ field.mul(v, scalar)
        return {p: v for p, v in result.items() if v}

    for diagonal in range(2 * length - 1):
        if diagonal == length - 1:
            for site in projector_sites(length, left, right):
                a, b = e(states, site), e(states, site + 1)
                states = combine(
                    (
                        (states, 1),
                        (a, 1),
                        (b, 1),
                        (e(b, site), INVERSE_DELTA),
                        (e(a, site + 1), INVERSE_DELTA),
                    )
                )
            continue
        for row in range(length):
            column = diagonal - row
            if 0 <= column < length:
                e_weight = value if row == column else 1
                states = combine(
                    ((states, Q), (e(states, length + row - 1 - column), e_weight))
                )
    return field.mul(
        power(Q, length - length * length // 2),
        states.get(four_paths[2] + four_paths[3], 0),
    )


def check():
    reference_controls = 0
    for length in range(2, 34, 2):
        for left, right in product(free_counts(length), repeat=2):
            assert reference_fusion(length, left, right) == fusion_rank(left + right)
            reference_controls += 1
    selected = {
        (4, 1, 1),
        (6, 0, 0),
        (6, 3, 3),
        (8, 2, 2),
        (8, 2, 5),
        (8, 5, 5),
        (8, 8, 8),
    }
    projector_controls, profiles, spans = [], [], []
    coefficient_checks, scalar_checks, polynomial_checks = 0, 0, 0
    field = Field(0b10011)
    for length in (4, 6, 8):
        rs = routes(length)
        ps = paths(length)
        d = len(ps)
        pieces = []
        for left, right in product(free_counts(length), repeat=2):
            sites = projector_sites(length, left, right)
            projector = projector_matrix(length, left, right, rs)
            for site in sites:
                for generator in (site, site + 1):
                    assert not np.any(apply_e(projector, rs[generator]))
                    assert not np.any(right_e(projector, rs[generator]))
            if sites:
                actual = rank(projector)
            else:
                assert np.array_equal(projector, np.eye(len(projector), dtype=np.uint8))
                actual = len(projector)
            expected = fusion_rank(left + right)
            assert actual == expected == reference_fusion(length, left, right)
            projector_controls.append(
                {
                    "length": length,
                    "free_left": left,
                    "free_right": right,
                    "exact_rank": actual,
                    "both_internal_annihilations": True,
                }
            )
            tensor = coefficient_tensor(length, left, right)
            pairs = [cut_ranks(v) for v in tensor]
            for b, (r, c) in enumerate(pairs):
                rb, cb = rank_bounds(length, left, right, b)
                assert r <= min(d * d, rb) and c <= min(d * d, cb)
                coefficient_checks += 1
                if np.any(tensor[b]):
                    pieces.extend((tensor[b], tensor[b].transpose(0, 2, 1, 3)))
            if (length, left, right) in selected:
                support = np.argwhere(np.any(tensor, axis=0))
                assert len(support)
                samples = list(dict.fromkeys((tuple(support[0]), tuple(support[-1]))))
                for ids in samples:
                    coefficients = tensor[(slice(None),) + ids]
                    for value in range(16):
                        evaluated = 0
                        for c in coefficients[::-1]:
                            evaluated = field.mul(evaluated, value) ^ int(c)
                        assert evaluated == direct_entry(
                            tuple(ps[i] for i in ids), left, right, value
                        )
                        scalar_checks += 1
                    polynomial_checks += 1
            profiles.append(
                {
                    "length": length,
                    "labels": d,
                    "free_left": left,
                    "free_right": right,
                    "central_rank": expected,
                    "rank_pairs_by_b": pairs,
                    "all_joint_bounds_pass": True,
                }
            )
        identity = np.zeros((d, d, d, d), dtype=np.uint8)
        swap = identity.copy()
        for p, q in product(range(d), repeat=2):
            identity[p, q, p, q] = swap[p, q, q, p] = 1
        pieces.extend((identity, identity.transpose(0, 2, 1, 3)))
        matrix = np.array(pieces).reshape(len(pieces), -1)
        dimension = rank(matrix)
        membership = rank(np.concatenate((matrix, swap.reshape(1, -1)))) == dimension
        assert membership == (length == 4)
        assert dimension == {4: 13, 6: 80, 8: 146}[length]
        spans.append(
            {
                "length": length,
                "labels": d,
                "projector_pairs": len(free_counts(length)) ** 2,
                "listed_nonzero_pieces": len(pieces),
                "span_dimension": dimension,
                "crossing_in_unfiltered_span": membership,
                "scope": "All prefix choices and degrees, both cut orientations and the two pairings; no asymptotic extrapolation.",
            }
        )
    sources = (
        "research/invertible-fusion-graded-seams.md",
        "tools/quantum_braid_graded_seam.py",
        "tools/check_quantum_braid_graded_seam.py",
        "tools/quantum_braid_fusion_seam.py",
        "tools/check_quantum_braid_fusion_seam.py",
        "tools/quantum_braid_nilpotent_seam.py",
        "tools/quantum_braid_path.py",
        "tools/check_fixed_field_descent.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_block_braid.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
    )
    return {
        "status": "pass",
        "field": "F16; order-five height-three vacuum quotient",
        "independent_fusion_multiplicity_controls": reference_controls,
        "global_projector_controls": projector_controls,
        "complete_coefficient_rank_pairs_checked": coefficient_checks,
        "independent_scalar_evaluations": scalar_checks,
        "independent_complete_sampled_polynomials": polynomial_checks,
        "profiles": profiles,
        "unfiltered_span_controls": spans,
        "scope": "Uniform programmable central rank factors and exact small controls; current prefix spans miss5/13labels; no stronger complete factoring rate.",
        "source_sha256": {
            f: sha256((ROOT / f).read_bytes()).hexdigest() for f in sources
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
