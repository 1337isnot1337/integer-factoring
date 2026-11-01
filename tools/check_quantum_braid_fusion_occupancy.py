#!/usr/bin/env python3
"""Replay the complete fusion-occupancy crossing and both actual cut ranks."""

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
from quantum_braid_fusion_occupancy import (
    central_projectors,
    central_rank,
    coefficient_tensor,
    rank_bounds,
)
from quantum_braid_fusion_seam import FIELD, INVERSE_DELTA, Q
from quantum_braid_path import EMBED, local_action
from quantum_braid_spectral_fusion import matrix_product

ROOT = Path(__file__).resolve().parents[1]


def reference_multiplicity(triplets):
    vacuum, charged = 1, 0
    for _ in range(triplets):
        vacuum, charged = charged, vacuum + charged
    return vacuum


def literal_entry(four_paths, u, y, multiply, inverses):
    """Independent scalar dictionary walk; never constructs a polynomial array."""
    n = len(four_paths[0])
    states = {four_paths[0] + four_paths[1]: 1}
    for _, row, column in sorted((j + k, j, k) for j in range(n) for k in range(n)):
        if row + column == n - 1:
            constant, coefficient = u, multiply[1 ^ u][inverses[7]]
        elif row == column:
            constant, coefficient = Q, y
        else:
            constant, coefficient = Q, 1
        following = {p: multiply[constant][v] for p, v in states.items() if constant}
        site = n + row - 1 - column
        for p, v in states.items():
            for target, scalar in local_action(p, site):
                value = multiply[coefficient][multiply[EMBED[scalar]][v]]
                following[target] = following.get(target, 0) ^ value
        states = {p: v for p, v in following.items() if v}
    scalar = 1
    for _ in range(n + n * n // 2):
        scalar = multiply[scalar][inverses[Q]]
    return multiply[scalar][states.get(four_paths[2] + four_paths[3], 0)]


def evaluate_entry(coefficients, u, y, multiply):
    answer = 0
    for row in coefficients[::-1]:
        value = 0
        for c in row[::-1]:
            value = multiply[value][y] ^ int(c)
        answer = multiply[answer][u] ^ value
    return answer


def check():
    field = Field(0b10011)
    multiply = tuple(tuple(field.mul(a, b) for b in range(16)) for a in range(16))
    inverses = tuple(field.power(a, 14) if a else 0 for a in range(16))
    assert np.array_equal(FIELD, multiply)
    projector_records, profiles = [], []
    placement_checks = rank_pairs = scalar_checks = 0
    for n in (4, 6, 8):
        ps, qs = paths(n), paths(2 * n)
        d, size = len(ps), len(qs)
        rs = routes(n)
        projectors = central_projectors(n)
        literal_sums = np.zeros_like(projectors)
        ranks_by_occupancy = [set() for _ in range(n + 1)]
        for mask in range(1 << n):
            image = np.eye(size, dtype=np.uint8)
            for pair in range(n):
                cup = FIELD[INVERSE_DELTA, apply_e(image, rs[2 * pair])]
                image = image ^ cup if mask >> pair & 1 else cup
            k = mask.bit_count()
            actual = rank(image)
            assert actual == reference_multiplicity(k)
            literal_sums[k] ^= image
            ranks_by_occupancy[k].add(actual)
            placement_checks += 1
        assert np.array_equal(literal_sums, projectors)
        assert np.array_equal(np.bitwise_xor.reduce(projectors, axis=0), np.eye(size))
        totals = [rank(p) for p in projectors]
        assert totals == [central_rank(n, k) for k in range(n + 1)]
        assert sum(totals) == size
        if n == 4:
            for i, j in product(range(n + 1), repeat=2):
                assert np.array_equal(
                    matrix_product(projectors[i], projectors[j]),
                    projectors[i] if i == j else np.zeros_like(projectors[0]),
                )
        if n <= 6:
            assert totals == [field.rank(p.tolist()) for p in projectors]
        projector_records.append(
            {
                "length": n,
                "vacuum_states": size,
                "literal_placements": 1 << n,
                "placement_ranks_by_occupancy": [sorted(s) for s in ranks_by_occupancy],
                "complete_sector_ranks": totals,
                "sum_is_identity": True,
            }
        )

        tensor = coefficient_tensor(n)
        weights = np.array([power(Q, 2 * k) for k in range(n + 1)], dtype=np.uint8)
        endpoint = np.bitwise_xor.reduce(
            FIELD[weights[:, None, None, None, None, None], tensor], axis=(0, 1)
        )
        swap = np.fromfunction(
            lambda p, q, r, s: (p == s) & (q == r), (d, d, d, d), dtype=int
        ).astype(np.uint8)
        assert np.array_equal(endpoint, swap)
        pairs = [[cut_ranks(a) for a in row] for row in tensor]
        for k, row in enumerate(pairs):
            for b, (r, c) in enumerate(row):
                rb, cb = rank_bounds(n, k, b)
                assert r <= rb and c <= cb, (n, k, b, r, c, rb, cb)
                rank_pairs += 1
        assert not np.any(tensor[1])
        samples = ((0, 0, 0, 0), (d - 1, 0, 0, d - 1), (0, 1, 0, 1))
        for ids in samples:
            coefficients = tensor[(slice(None), slice(None)) + ids]
            for u, y in product(range(16), repeat=2):
                assert evaluate_entry(coefficients, u, y, multiply) == literal_entry(
                    tuple(ps[i] for i in ids), u, y, multiply, inverses
                ), (n, ids, u, y)
                scalar_checks += 1
        profiles.append(
            {
                "length": n,
                "labels": d,
                "degree_by_variable": [n, n],
                "all_coefficient_rank_pairs": pairs,
                "all_joint_bounds_pass": True,
                "complete_scale_one_crossing": True,
                "identically_zero_single_triplet_sector": True,
                "complete_polynomial_entries_independently_certified": len(samples),
            }
        )
    sources = (
        "research/fusion-occupancy-resolution.md",
        "tools/quantum_braid_fusion_occupancy.py",
        "tools/check_quantum_braid_fusion_occupancy.py",
        "tools/quantum_braid_graded_seam.py",
        "tools/quantum_braid_fusion_seam.py",
        "tools/check_quantum_braid_fusion_seam.py",
        "tools/quantum_braid_spectral_fusion.py",
        "tools/quantum_braid_path.py",
        "tools/paired_path_crossing.py",
        "tools/check_fixed_field_descent.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
    )
    return {
        "status": "pass",
        "literal_placement_projectors_checked": placement_checks,
        "complete_coefficient_rank_pairs_checked": rank_pairs,
        "independent_bivariate_scalar_evaluations": scalar_checks,
        "complete_sampled_polynomials_certified": 9,
        "projector_controls": projector_records,
        "complete_profiles": profiles,
        "scope": "Uniform complete occupancy crossing resolution and exact central ranks; no favorable full rank sum or stronger factoring score.",
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
