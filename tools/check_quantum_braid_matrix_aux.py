#!/usr/bin/env python3
"""Exact endpoints, leading factors and diagram controls for a matrix auxiliary."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from check_quantum_braid_path_profiles import cut_ranks, rank
from quantum_block_braid import braid_word
from quantum_braid_matrix_aux import (
    IDENTITY,
    Q,
    coefficient_tensor,
    leading_entry,
    matrix_product,
    normalization_row,
)
from quantum_braid_path import local_action

ROOT = Path(__file__).resolve().parents[1]


def diagram(length):
    """Literal TL connectivity; inputs first, outputs second; count removed loops."""
    m = 2 * length
    partners = list(range(m, 2 * m)) + list(range(m))
    loops = 0
    for site in braid_word(length, length)[1::2]:
        a, b = m + site, m + site + 1
        x, y = partners[a], partners[b]
        if x == b:
            loops += 1
        else:
            partners[x], partners[y] = y, x
            partners[a], partners[b] = b, a
    return partners, loops


def expected_pairs(n):
    return (
        [(i, i + 1) for i in range(0, n, 2)]
        + [(n + i, 2 * n + i) for i in range(n - 1)]
        + [(2 * n - 1, 4 * n - 1), (3 * n - 1, 3 * n)]
        + [(3 * n + i, 3 * n + i + 1) for i in range(1, n - 1, 2)]
    )


def zigzag_closure_loops(partners):
    remaining = set(range(len(partners)))
    loops = 0
    while remaining:
        stack = [next(iter(remaining))]
        loops += 1
        while stack:
            i = stack.pop()
            if i in remaining:
                remaining.remove(i)
                stack.extend((partners[i], i ^ 1))
    return loops


def direct_entry(four_paths, value, field, embedding):
    """Scalar dictionary propagation over F256, independent of the batched XORs."""
    n = len(four_paths[0])
    states = {four_paths[0] + four_paths[1]: (1, 0)}
    a = ((1 ^ value, value), (value, 1 ^ field.mul(embedding[3], value)))
    b = ((value, value ^ field.mul(value, value)), (0, value))

    def apply(matrix, vector):
        return tuple(
            field.mul(row[0], vector[0]) ^ field.mul(row[1], vector[1])
            for row in matrix
        )

    for site in braid_word(n, n):
        following = {}
        for p, vector in states.items():
            targets = [(p, apply(a, vector))]
            e_vector = apply(b, vector)
            for q, c in local_action(p, site):
                targets.append((q, tuple(field.mul(embedding[c], x) for x in e_vector)))
            for q, v in targets:
                old = following.get(q, (0, 0))
                following[q] = old[0] ^ v[0], old[1] ^ v[1]
        states = {p: v for p, v in following.items() if any(v)}
    vector = states.get(four_paths[2] + four_paths[3], (0, 0))
    left = normalization_row(n)
    return field.mul(embedding[int(left[0])], vector[0]) ^ field.mul(
        embedding[int(left[1])], vector[1]
    )


def check():
    q_power = IDENTITY.copy()
    for _ in range(5):
        q_power = matrix_product(q_power, Q)
    assert np.array_equal(q_power, IDENTITY)
    diagram_cases = 0
    for n in range(2, 130, 2):
        partners, loops = diagram(n)
        assert loops == 0
        assert {(i, p) for i, p in enumerate(partners) if i < p} == set(
            expected_pairs(n)
        )
        assert zigzag_closure_loops(partners) == n
        diagram_cases += 1
    field = Field(0x11B)
    alpha = next(a for a in range(2, 256) if field.mul(a, a) ^ a ^ 1 == 0)
    embedding = (0, 1, alpha, alpha ^ 1)
    random = Random(891204)
    profiles, leading_controls, entry_controls = [], 0, 0
    for n, span in ((4, 16), (6, 55), (8, 97)):
        ps = paths(n)
        d = len(ps)
        tensor = coefficient_tensor(n)
        assert len(tensor) - 1 == 3 * n * n // 2
        identity, swap = np.zeros(tensor.shape[1:], dtype=np.uint8), np.zeros(
            tensor.shape[1:], dtype=np.uint8
        )
        for p, q in product(range(d), repeat=2):
            identity[p, q, p, q] = swap[p, q, q, p] = 1
        assert np.array_equal(tensor[0], identity)
        assert np.array_equal(np.bitwise_xor.reduce(tensor, axis=0), swap)
        for ids in product(range(d), repeat=4):
            assert tensor[(-1,) + ids] == leading_entry(tuple(ps[i] for i in ids))
            leading_controls += 1
        assert cut_ranks(tensor[-1]) == [d, d]
        actual_span = rank(tensor.reshape(len(tensor), -1))
        assert actual_span == span
        if n <= 6:
            samples = [(d - 1, 0, 0, d - 1)] + [
                tuple(random.randrange(d) for _ in range(4)) for _ in range(2)
            ]
            # 64 distinct values exceed the degree (at most54), so these are
            # independent COMPLETE polynomial controls for each sampled entry.
            for ids in samples:
                coefficients = tensor[(slice(None),) + ids]
                for value in range(64):
                    expected = 0
                    for c in coefficients[::-1]:
                        expected = field.mul(expected, value) ^ embedding[int(c)]
                    assert expected == direct_entry(
                        tuple(ps[i] for i in ids), value, field, embedding
                    )
                    entry_controls += 1
        profiles.append(
            {
                "length": n,
                "labels": d,
                "degree": len(tensor) - 1,
                "coefficient_span_dimension": actual_span,
                "leading_cut_ranks": cut_ranks(tensor[-1]),
                "leading_nonzero_entries": int(np.count_nonzero(tensor[-1])),
                "strict_endpoints": True,
            }
        )
    sources = (
        "research/quantum-braid-matrix-auxiliary.md",
        "tools/check_quantum_braid_matrix_aux.py",
        "tools/quantum_braid_matrix_aux.py",
        "tools/quantum_braid_path.py",
        "tools/quantum_block_braid.py",
        "tools/paired_path_crossing.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
        "tools/check_quantum_braid_path_profiles.py",
        "tools/check_fixed_field_descent.py",
        "tools/search_auxiliary_extension_span.py",
    )
    return {
        "status": "pass",
        "field": "F4; independent entry controls over F256",
        "diagram_lengths_checked": diagram_cases,
        "leading_tensor_entries_checked": leading_controls,
        "independent_scalar_evaluations": entry_controls,
        "independent_complete_polynomial_entries": 6,
        "profiles": profiles,
        "scope": "Uniform auxiliary construction and leading factorization; finite spans only; no favorable whole-family rank sum or factoring improvement.",
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
