#!/usr/bin/env python3
"""Exact finite controls for complementary asymmetric nilpotent-seam edges."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from check_quantum_braid_path_profiles import cut_ranks, rank
from quantum_braid_nilpotent_seam import (
    FIELD,
    IDENTITY,
    N,
    Q,
    coefficient_tensor,
    matrix_product,
    normalization_row,
)
from quantum_braid_path import EMBED, local_action
from quantum_braid_rank_one import add_row, pack

ROOT = Path(__file__).resolve().parents[1]


def auxiliary_controls(random):
    assert np.array_equal(matrix_product(Q, Q), IDENTITY ^ FIELD[2, Q])
    qp, npower = IDENTITY.copy(), IDENTITY.copy()
    for _ in range(5):
        qp = matrix_product(qp, Q)
    for _ in range(4):
        npower = matrix_product(npower, N)
    assert np.array_equal(qp, IDENTITY) and not np.any(npower)
    assert matrix_product(matrix_product(matrix_product(N, N), N), Q)[0, 0] == 1
    count = 0
    # Each central e_i has eigenvalues0,alpha. Check the complete central
    # auxiliary polynomial on independently chosen simultaneous eigensectors.
    for length in (2, 4, 6, 8, 10, 12, 16, 20, 24, 32):
        for _ in range(8):
            polynomial = [IDENTITY.copy()]
            for epsilon in [random.choice((0, 2)) for _ in range(length)]:
                local = (
                    IDENTITY,
                    Q ^ IDENTITY ^ FIELD[epsilon, IDENTITY ^ N],
                    FIELD[epsilon, N],
                )
                following = [
                    np.zeros((4, 4), dtype=np.uint8) for _ in range(len(polynomial) + 2)
                ]
                for i, a in enumerate(polynomial):
                    for j, b in enumerate(local):
                        following[i + j] ^= matrix_product(b, a)
                polynomial = following
            bound = 2 * length - length // 4
            assert not np.any(polynomial[bound + 1 :])
            count += 1
    field = Field(0b111)
    for _ in range(24):
        matrix = np.array(
            [[random.randrange(4) for _ in range(13)] for _ in range(9)], dtype=np.uint8
        )
        assert rank(matrix) == field.rank(matrix.tolist())
    return count


def cup_image_controls():
    """Full global path matrices, not only the middle-zero boundary block."""
    records = []
    for length in (4, 6, 8):
        ps = paths(2 * length)
        indices = {p: i for i, p in enumerate(ps)}
        free = length // 4
        cases = 0
        for selected in combinations(range(length), length - free):
            matrix = np.eye(len(ps), dtype=np.uint8)
            for pair in selected:
                diagonal = np.zeros(len(ps), dtype=np.uint8)
                source, target = [], []
                for i, p in enumerate(ps):
                    for q, c in local_action(p, 2 * pair):
                        if q == p:
                            diagonal[i] = c
                        else:
                            source.append(i)
                            target.append(indices[q])
                source, target = np.array(source, dtype=int), np.array(
                    target, dtype=int
                )
                following = FIELD[diagonal[:, None], matrix]
                following[target] ^= FIELD[diagonal[source, None], matrix[source]]
                matrix = following
            assert rank(matrix) == len(paths(2 * free))
            cases += 1
        records.append(
            {
                "length": length,
                "free_pairs": free,
                "placements_checked": cases,
                "exact_projector_rank": len(paths(2 * free)),
            }
        )
    return records


def direct_entry(four_paths, x, y):
    """Literal scalar propagation over F16; no batched indexed XOR operations."""
    field = Field(0b10011)
    length = len(four_paths[0])
    q = [[EMBED[int(v)] for v in row] for row in Q]
    nil = N.tolist()
    states = {four_paths[0] + four_paths[1]: (1, 0, 0, 0)}

    def apply(a, vector):
        result = []
        for row in a:
            value = 0
            for c, v in zip(row, vector):
                value ^= field.mul(c, v)
            result.append(value)
        return tuple(result)

    for diagonal in range(2 * length - 1):
        for row in range(length):
            column = diagonal - row
            if not 0 <= column < length:
                continue
            site = length + row - 1 - column
            if diagonal == length - 1:
                a = [
                    [
                        int(i == j) ^ field.mul(x, q[i][j] ^ int(i == j))
                        for j in range(4)
                    ]
                    for i in range(4)
                ]
                b = [
                    [
                        field.mul(x, int(i == j))
                        ^ field.mul(x ^ field.mul(x, x), nil[i][j])
                        for j in range(4)
                    ]
                    for i in range(4)
                ]
            else:
                a = q
                scalar = y if row == column else 1
                b = [[scalar * int(i == j) for j in range(4)] for i in range(4)]
            following = {}
            for p, vector in states.items():
                targets = [(p, apply(a, vector))]
                e_vector = apply(b, vector)
                for s, c in local_action(p, site):
                    targets.append((s, tuple(field.mul(EMBED[c], v) for v in e_vector)))
                for s, v in targets:
                    old = following.get(s, (0, 0, 0, 0))
                    following[s] = tuple(a ^ b for a, b in zip(old, v))
            states = {p: v for p, v in following.items() if any(v)}
    value = 0
    for c, v in zip(
        normalization_row(length),
        states.get(four_paths[2] + four_paths[3], (0, 0, 0, 0)),
    ):
        value ^= field.mul(EMBED[int(c)], v)
    return value


def evaluate(coefficients, x, y):
    field = Field(0b10011)
    result = 0
    for row in coefficients[::-1]:
        value = 0
        for c in row[::-1]:
            value = field.mul(value, y) ^ EMBED[int(c)]
        result = field.mul(result, x) ^ value
    return result


def strip_profiles(tensor, identity, swap):
    a, b = tensor.shape[:2]
    space = {}
    add_row(space, *pack(identity.reshape(-1)))
    target = pack(swap.reshape(-1))
    records = []
    for radius in range(b):
        for i, j in product(range(a), range(b)):
            if min(a - 1 - i, j) == radius:
                add_row(space, *pack(tensor[i, j].reshape(-1)))
        probe = space.copy()
        add_row(probe, *target)
        member = len(probe) == len(space)
        records.append(
            {
                "width": radius,
                "span_with_identity": len(space),
                "crossing_in_span": member,
            }
        )
        if member:
            break
    assert records[-1]["crossing_in_span"]
    return records


def check():
    random = Random(912438)
    degree_controls = auxiliary_controls(random)
    cup_controls = cup_image_controls()
    profiles, scalar_controls = [], 0
    for length, expected_span, first_width in ((4, 10, 2), (6, 69, 6), (8, 118, 8)):
        ps = paths(length)
        d = len(ps)
        tensor = coefficient_tensor(length)
        a, b = tensor.shape[:2]
        identity = np.zeros(tensor.shape[2:], dtype=np.uint8)
        swap = identity.copy()
        for p, q in product(range(d), repeat=2):
            identity[p, q, p, q] = swap[p, q, q, p] = 1
        assert np.array_equal(
            np.bitwise_xor.reduce(tensor.reshape(a * b, d, d, d, d), axis=0), swap
        )
        assert rank(tensor.reshape(a * b, -1)) == expected_span
        top_pairs = [cut_ranks(v) for v in tensor[-1]]
        low_pairs = [cut_ranks(v) for v in tensor[:, 0]]
        for j, pair in enumerate(top_pairs):
            assert pair[0] <= 8 * (length + 1) ** 3 * (j + 1) * len(
                paths(2 * (length // 4))
            )
        for i, pair in enumerate(low_pairs):
            assert pair[1] <= 16 * (i + 1)
        strips = strip_profiles(tensor, identity, swap)
        assert strips[-1]["width"] == first_width
        actual_degree = max(i for i in range(a) if np.any(tensor[i]))
        if length <= 6:
            for ids in ((d - 1, 0, 0, d - 1), (0, 1, 1, 0)):
                coefficients = tensor[(slice(None), slice(None)) + ids]
                # The grid has more points in each axis than the respective
                # degree bounds, hence checks the entire sampled bivariate entry.
                for x, y in product(range(16), range(8)):
                    assert evaluate(coefficients, x, y) == direct_entry(
                        tuple(ps[i] for i in ids), x, y
                    )
                    scalar_controls += 1
        profiles.append(
            {
                "length": length,
                "labels": d,
                "central_degree_bound": a - 1,
                "actual_central_degree": actual_degree,
                "boundary_degree_bound": b - 1,
                "coefficient_span_dimension": expected_span,
                "highest_allowed_central_edge_rank_pairs": top_pairs,
                "zero_boundary_edge_rank_pairs": low_pairs,
                "strip_membership": strips,
                "crossing_endpoint": True,
                "strict_exposure": "(1+t)I+t*G(t,t^(2n+1))",
            }
        )
    sources = [
        "research/nilpotent-seam-crossings.md",
        "tools/quantum_braid_nilpotent_seam.py",
        "tools/check_quantum_braid_nilpotent_seam.py",
        "tools/quantum_braid_rank_one.py",
        "tools/quantum_braid_path.py",
        "tools/check_quantum_braid_path_profiles.py",
        "tools/check_fixed_field_descent.py",
        "tools/search_auxiliary_extension_span.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_block_braid.py",
    ]
    return {
        "status": "pass",
        "field": "F4; independent sampled entry controls over F16",
        "central_eigensector_polynomial_degree_controls": degree_controls,
        "mandatory_cup_projector_controls": cup_controls,
        "packed_rank_reference_controls": 24,
        "independent_scalar_evaluations": scalar_controls,
        "independent_complete_bivariate_entries": 4,
        "profiles": profiles,
        "scope": "Uniform complementary asymmetric edge bounds, not a whole-crossing decomposition; finite edge spans miss the crossing; no stronger factoring score.",
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
