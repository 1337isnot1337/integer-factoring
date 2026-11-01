#!/usr/bin/env python3
"""Vacuum-interface factors and complete triplet-seam coefficient profiles."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from paired_path_crossing import power
from quantum_braid_fusion_seam import (
    FIELD,
    INVERSE_DELTA,
    Q,
    coefficient_tensor,
)
from quantum_braid_path import EMBED, local_action
from quantum_braid_nilpotent_seam import coefficient_tensor as nilpotent_tensor

ROOT = Path(__file__).resolve().parents[1]
INVERSE = np.array([0] + [power(a, 14) for a in range(1, 16)], dtype=np.uint8)


def rank(matrix):
    a = matrix.copy()
    r = 0
    for column in range(a.shape[1]):
        nonzero = np.flatnonzero(a[r:, column])
        if not len(nonzero):
            continue
        p = r + nonzero[0]
        a[[r, p]] = a[[p, r]]
        a[r] = FIELD[INVERSE[a[r, column]], a[r]]
        factors = a[:, column].copy()
        factors[r] = 0
        a ^= FIELD[factors[:, None], a[r, None, :]]
        r += 1
        if r == len(a):
            break
    return r


def cut_ranks(tensor):
    d = tensor.shape[0]
    return [
        rank(tensor.reshape(d * d, d * d)),
        rank(tensor.transpose(0, 2, 1, 3).reshape(d * d, d * d)),
    ]


def routes(length):
    ps = paths(2 * length)
    indices = {p: i for i, p in enumerate(ps)}
    result = []
    for site in range(2 * length - 1):
        diagonal = np.zeros(len(ps), dtype=np.uint8)
        source, target = [], []
        for i, p in enumerate(ps):
            for q, c in local_action(p, site):
                if q == p:
                    diagonal[i] = EMBED[c]
                else:
                    source.append(i)
                    target.append(indices[q])
        result.append(
            (diagonal, np.array(source, dtype=int), np.array(target, dtype=int))
        )
    return result


def apply_e(matrix, route):
    diagonal, source, target = route
    result = FIELD[diagonal[:, None], matrix]
    result[target] ^= FIELD[diagonal[source, None], matrix[source]]
    return result


def projector_controls():
    records = []
    fusion_vacuum, fusion_two = 1, 0
    for length in range(1, 9):
        fusion_vacuum, fusion_two = fusion_two, fusion_vacuum + fusion_two
        rs = routes(length)
        projector = np.eye(len(paths(2 * length)), dtype=np.uint8)
        for site in range(0, 2 * length, 2):
            projector ^= FIELD[INVERSE_DELTA, apply_e(projector, rs[site])]
        actual = rank(projector)
        assert actual == fusion_vacuum
        if length % 2 == 0:
            assert actual == len(paths(length))
        records.append(
            {
                "central_pairs": length,
                "global_vacuum_dimension": len(projector),
                "exact_projector_rank": actual,
            }
        )
    return records


def fixed_subset_tensor(length, selected, rs):
    ps, qs = paths(length), paths(2 * length)
    indices = {p: i for i, p in enumerate(qs)}
    d = len(ps)
    boundary = np.array([indices[p + q] for p in ps for q in ps])
    states = np.zeros((len(qs), d * d), dtype=np.uint8)
    states[boundary, np.arange(d * d)] = 1
    for _, j, k in sorted((j + k, j, k) for j in range(length) for k in range(length)):
        action = apply_e(states, rs[length + j - 1 - k])
        if j + k == length - 1:
            states ^= FIELD[INVERSE_DELTA, action]
        elif j == k:
            states = action if j in selected else FIELD[Q, states]
        else:
            states = FIELD[Q, states] ^ action
    result = states[boundary].T.reshape(d, d, d, d)
    return FIELD[power(Q, length - length * length // 2), result]


def direct_entry(four_paths, value):
    field = Field(0b10011)
    length = len(four_paths[0])
    states = {four_paths[0] + four_paths[1]: 1}
    for diagonal in range(2 * length - 1):
        for row in range(length):
            column = diagonal - row
            if not 0 <= column < length:
                continue
            site = length + row - 1 - column
            central = diagonal == length - 1
            own = 1 if central else Q
            e_weight = INVERSE_DELTA if central else value if row == column else 1
            following = {}
            for p, v in states.items():
                following[p] = following.get(p, 0) ^ field.mul(own, v)
                for q, c in local_action(p, site):
                    term = field.mul(field.mul(v, e_weight), EMBED[c])
                    following[q] = following.get(q, 0) ^ term
            states = {p: v for p, v in following.items() if v}
    return field.mul(
        power(Q, length - length * length // 2),
        states.get(four_paths[2] + four_paths[3], 0),
    )


def mixed_span_profile(length, triplet):
    nilpotent = np.array(EMBED, dtype=np.uint8)[nilpotent_tensor(length)]
    top, d = nilpotent.shape[0] - 1, nilpotent.shape[2]
    pieces = []
    for gap in range(3):
        for b in range(length + 1):
            tensor = nilpotent[top - gap, b]
            if np.any(tensor):
                pieces.extend((tensor, tensor.transpose(0, 2, 1, 3)))
    for b in range((length + 4) // 5 + 1):
        tensor = triplet[b]
        if np.any(tensor):
            pieces.extend((tensor, tensor.transpose(0, 2, 1, 3)))
    identity = np.zeros((d, d, d, d), dtype=np.uint8)
    swap = identity.copy()
    for p, q in product(range(d), repeat=2):
        identity[p, q, p, q] = swap[p, q, q, p] = 1
    pieces.extend((identity, identity.transpose(0, 2, 1, 3)))
    matrix = np.array(pieces).reshape(len(pieces), -1)
    dimension = rank(matrix)
    membership = rank(np.concatenate((matrix, swap.reshape(1, -1)))) == dimension
    assert membership == (length == 4)
    return {
        "nilpotent_central_gaps": [0, 1, 2],
        "triplet_boundary_limit": (length + 4) // 5,
        "includes_cut_reflections_and_two_pairings": True,
        "listed_nonzero_pieces": len(pieces),
        "span_dimension": dimension,
        "crossing_in_span": membership,
        "scope": "Complete finite full-alphabet span test; no extrapolation to larger widths.",
    }


def check():
    random, field = Random(771932), Field(0b10011)
    for _ in range(24):
        a = np.array(
            [[random.randrange(16) for _ in range(11)] for _ in range(9)],
            dtype=np.uint8,
        )
        assert rank(a) == field.rank(a.tolist())
    projector = projector_controls()
    assert 21**2 * 5**25 < 2**67
    profiles, subsets_checked, scalar_checks = [], 0, 0
    for length in (4, 6, 8):
        ps = paths(length)
        d = len(ps)
        tensor = coefficient_tensor(length)
        rs = routes(length)
        sums = np.zeros_like(tensor)
        subset_ranks = []
        for b in range(length + 1):
            ranks = []
            for selected in combinations(range(length), b):
                fixed = fixed_subset_tensor(length, selected, rs)
                c = cut_ranks(fixed)[1]
                assert c <= len(paths(2 * b))
                ranks.append(c)
                sums[b] ^= fixed
                subsets_checked += 1
            subset_ranks.append(
                {
                    "boundary_gates": b,
                    "subsets": len(ranks),
                    "column_ranks": sorted(set(ranks)),
                    "vacuum_bound": len(paths(2 * b)),
                }
            )
        assert np.array_equal(sums, tensor)
        pairs = [cut_ranks(v) for v in tensor]
        for b, (r, c) in enumerate(pairs):
            assert r <= min(b + 1, length - b + 1) * d
            assert c <= comb(length, b) * len(paths(2 * b))
        for ids in ((d - 1, 0, 0, d - 1), (0, 1, 1, 0)):
            coefficients = tensor[(slice(None),) + ids]
            for value in range(16):
                evaluated = 0
                for c in coefficients[::-1]:
                    evaluated = field.mul(evaluated, value) ^ int(c)
                assert evaluated == direct_entry(tuple(ps[i] for i in ids), value)
                scalar_checks += 1
        profiles.append(
            {
                "length": length,
                "labels": d,
                "boundary_degree": length,
                "rank_pairs_by_b": pairs,
                "complete_subset_interfaces": subset_ranks,
                "candidate_middle_boundary_limit": (length + 4) // 5,
                "all_joint_bounds_pass": True,
                "mixed_candidate_span": mixed_span_profile(length, tensor),
            }
        )
    sources = (
        "research/fusion-projector-seams.md",
        "tools/quantum_braid_fusion_seam.py",
        "tools/check_quantum_braid_fusion_seam.py",
        "tools/quantum_braid_path.py",
        "tools/quantum_braid_nilpotent_seam.py",
        "tools/check_fixed_field_descent.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_block_braid.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
    )
    return {
        "status": "pass",
        "field": "F16 with the height-three order-five path quotient",
        "independent_field_rank_controls": 24,
        "global_triplet_projector_controls": projector,
        "complete_selected_subset_interfaces": subsets_checked,
        "independent_scalar_evaluations": scalar_checks,
        "independent_complete_sampled_polynomials": 6,
        "profiles": profiles,
        "middle_window_integer_certificate": {"lhs": 21**2 * 5**25, "rhs": 2**67},
        "scope": "Uniform vacuum-interface and fusion-projector rank factors; complete finite controls; mixed crossing span and stronger factoring rate remain unproved.",
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
