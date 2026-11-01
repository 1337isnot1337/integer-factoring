#!/usr/bin/env python3
"""Complete F4 path-braid profiles and ambient reference controls; needs NumPy."""

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
from paired_path_crossing import MULTIPLY as M16, power
from quantum_block_braid import polynomial_entry as ambient_entry
from quantum_braid_path import (
    EMBED,
    MULTIPLY as M4,
    coefficient_tensor,
    crossing_weights,
    polynomial_entry,
)
from search_auxiliary_extension_span import extension_rank

ROOT = Path(__file__).resolve().parents[1]
MULTIPLY = np.array(M4, dtype=np.uint8)


def rank(matrix):
    # The existing packed eliminator uses basis (alpha,alpha^2), with 3=1.
    rows = []
    for row in matrix:
        low = int.from_bytes(
            np.packbits(row & 1, bitorder="little").tobytes(), "little"
        )
        high = int.from_bytes(
            np.packbits(row >> 1, bitorder="little").tobytes(), "little"
        )
        rows.append((low ^ high, low))
    return extension_rank(rows)


def cut_ranks(tensor):
    d = tensor.shape[0]
    return [
        rank(tensor.reshape(d * d, d * d)),
        rank(tensor.transpose(0, 2, 1, 3).reshape(d * d, d * d)),
    ]


def evaluate(coefficients, value):
    result = np.zeros(coefficients.shape[1:], dtype=np.uint8)
    for coefficient in coefficients[::-1]:
        result = MULTIPLY[value, result] ^ coefficient
    return result


def check():
    random = Random(771)
    field = Field(0b111)
    for _ in range(32):
        n, m = random.randrange(1, 15), random.randrange(1, 20)
        a = np.array(
            [[random.randrange(4) for _ in range(m)] for _ in range(n)], dtype=np.uint8
        )
        assert rank(a) == field.rank(a.tolist())
    profiles, ambient_entries, coefficient_controls, scalar_controls = [], 0, 0, 0
    for n, expected_a, expected_b in ((4, 5, 5), (6, 33, 29), (8, 62, 51)):
        ps = paths(n)
        d = len(ps)
        a = coefficient_tensor(n)
        weights = np.array(crossing_weights(n), dtype=np.uint8)
        b = MULTIPLY[weights[:, None, None, None, None], a]
        identity, crossing = np.zeros(a.shape[1:], dtype=np.uint8), np.zeros(
            a.shape[1:], dtype=np.uint8
        )
        for p, q in product(range(d), repeat=2):
            identity[p, q, p, q] = 1
            crossing[p, q, q, p] = 1
        assert np.array_equal(b[0], identity)
        assert np.array_equal(np.bitwise_xor.reduce(b, axis=0), crossing)
        assert rank(a.reshape(len(a), -1)) == expected_a
        assert rank(b.reshape(len(b), -1)) == expected_b
        if n <= 6:
            phase = power(8, -n * n // 2)
            for ids in product(range(d), repeat=4):
                expected = ambient_entry(tuple(ps[i] for i in ids))
                old = [
                    M16[power(8, n * n // 2 - j)][EMBED[int(x)]]
                    for j, x in enumerate(a[(slice(None),) + ids])
                ]
                actual = [0] * (n * n + 2)
                for j, value in enumerate(old):
                    actual[j] ^= M16[phase][value]
                    actual[j + 1] ^= M16[1 ^ phase][value]
                assert actual == expected
                ambient_entries += 1
                coefficient_controls += len(actual)
        for _ in range(8):
            ids = tuple(random.randrange(d) for _ in range(4))
            assert b[(slice(None),) + ids].tolist() == polynomial_entry(
                tuple(ps[i] for i in ids)
            )
            scalar_controls += 1
        raw_pairs = [cut_ranks(t) for t in a]
        crossing_pairs = [cut_ranks(t) for t in b]
        assert all(raw_pairs[j][0] == raw_pairs[n * n - j][1] for j in range(n * n + 1))
        profiles.append(
            {
                "length": n,
                "labels": d,
                "path_states": len(paths(2 * n)),
                "full_boundary_entries": d**4,
                "unprojected_span_dimension": expected_a,
                "crossing_span_dimension": expected_b,
                "retained_coefficient_indices": [j for j, w in enumerate(weights) if w],
                "unprojected_rank_pairs": raw_pairs,
                "crossing_rank_pairs": crossing_pairs,
                "central_rank_pair": raw_pairs[n * n // 2],
                "full_alphabet_rank_reversal": True,
                "F4_evaluation_profiles": [
                    {
                        "value": s,
                        "unprojected": cut_ranks(evaluate(a, s)),
                        "crossing": cut_ranks(evaluate(b, s)),
                    }
                    for s in range(4)
                ],
            }
        )
    files = (
        "tools/check_quantum_braid_path_profiles.py",
        "tools/quantum_braid_path.py",
        "tools/quantum_block_braid.py",
        "tools/paired_path_crossing.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
        "tools/check_fixed_field_descent.py",
        "tools/search_auxiliary_extension_span.py",
    )
    return {
        "status": "pass",
        "field": "F2[alpha]/(alpha^2+alpha+1)",
        "packed_rank_reference_controls": 32,
        "full_ambient_polynomial_entries": ambient_entries,
        "full_ambient_coefficient_controls": coefficient_controls,
        "batch_scalar_entry_controls": scalar_controls,
        "profiles": profiles,
        "scope": "Complete finite path profiles and ambient comparisons; no asymptotic rank theorem or finite recombination optimum.",
        "source_sha256": {
            f: sha256((ROOT / f).read_bytes()).hexdigest() for f in files
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
