#!/usr/bin/env python3
"""Exact reversal alphabets, leading factors and finite nilpotent dressing controls."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from bounded_spin_crossing import paths, suffix_counts
from check_fixed_field_descent import Field
from paired_path_crossing import (
    MULTIPLY as M16,
    path_entry,
    power,
    splitting_blocks,
    step,
)
from quantum_braid_path import MULTIPLY as M4, coefficient_tensor, crossing_weights
from quantum_braid_reversal import (
    binary_layout,
    binary_path,
    reverse_label,
    reverse_path,
    diagram_norm,
    leading_entry,
    nilpotent_dressing,
)
from check_quantum_braid_path_profiles import cut_ranks, evaluate

ROOT = Path(__file__).resolve().parents[1]
FIELD = np.array(M4, dtype=np.uint8)


def identity_tensor(d):
    result = np.zeros((d,) * 4, dtype=np.uint8)
    for p, q in product(range(d), repeat=2):
        result[p, q, p, q] = 1
    return result


def crossing_tensor(d):
    result = np.zeros((d,) * 4, dtype=np.uint8)
    for p, q in product(range(d), repeat=2):
        result[p, q, q, p] = 1
    return result


def check():
    local, words, labels = 0, 0, 0
    for q in range(2, 16):
        try:
            blocks = splitting_blocks(q)
        except ValueError:
            continue

        def integer(k):
            return M16[power(q, k) ^ power(q, -k)][power(q ^ power(q, -1), -1)]

        def potential(h, r):
            value = power(q, (r + 1) * (h - r))
            for k in range(1, r + 1):
                value = M16[value][integer(k)]
            for k in range(1, h - r + 1):
                value = M16[value][integer(k)]
            return value

        for h in range(4):
            for direction in (1, -1):
                if not 0 <= h + direction <= 3:
                    continue
                for r in range(h + 1):
                    for bit in (0, 1):
                        result = step(h, direction, r, bit, blocks)
                        if result is None:
                            continue
                        nr, c = result
                        following = step(h + direction, -direction, nr, 1 - bit, blocks)
                        assert following is not None and following[0] == r
                        expected = M16[c][
                            M16[potential(h + direction, nr)][
                                power(potential(h, r), -1)
                            ]
                        ]
                        assert following[1] == expected
                        local += 1
    for n in range(2, 11, 2):
        for p in paths(n):
            reverse = reverse_path(p)
            assert diagram_norm(p) == diagram_norm(reverse)
            for word in range(1 << n):
                reflected = sum((1 - (word >> i & 1)) << (n - 1 - i) for i in range(n))
                assert path_entry(p, reflected) == path_entry(reverse, word)
                words += 1
    for width in range(1, 13):
        ps = [binary_path(width, i) for i in range(1 << width)]
        assert len(set(ps)) == 1 << width
        n, _ = binary_layout(width)
        assert suffix_counts(n)[n][0] >= 1 << width
        assert n == 4 or suffix_counts(n - 2)[n - 2][0] < 1 << width
        for j, p in enumerate(ps):
            assert ps[reverse_label(width, j)] == reverse_path(p)
            assert len(p) == n and sum(p) == 0
            height = 0
            for s in p:
                height += s
                assert 0 <= height <= 3
            labels += 1
    full_leading_entries, polynomial_controls, dressings = 0, 0, []
    field = Field(0b111)
    for n, width, scalar in ((4, 1, 1), (6, 2, 3), (8, 3, 3)):
        ps = paths(n)
        a = coefficient_tensor(n)
        expected = np.array(
            [leading_entry(four) for four in product(ps, repeat=4)], dtype=np.uint8
        ).reshape((len(ps),) * 4)
        assert np.array_equal(a[-1], expected)
        full_leading_entries += len(ps) ** 4
        b = FIELD[
            np.array(crossing_weights(n), dtype=np.uint8)[:, None, None, None, None], a
        ]
        selected = [ps.index(binary_path(width, i)) for i in range(1 << width)]
        d = len(selected)
        selected_tensor = b[:, selected][:, :, selected][:, :, :, selected][
            :, :, :, :, selected
        ]
        row = selected_tensor[-1].transpose(0, 1, 3, 2).reshape(d * d, d * d)
        assert not any(x for r in field.matmul(row.tolist(), row.tolist()) for x in r)
        assert cut_ranks(selected_tensor[-1]) == [1, d * d]
        dressed = nilpotent_dressing(selected_tensor, scalar)
        assert np.array_equal(dressed[0], identity_tensor(d))
        assert np.array_equal(
            np.bitwise_xor.reduce(dressed, axis=0), crossing_tensor(d)
        )
        if n <= 6:
            original_rows = selected_tensor.transpose(0, 1, 2, 4, 3).reshape(
                len(selected_tensor), d * d, d * d
            )
            u = [
                np.eye(d * d, dtype=np.uint8).tolist(),
                FIELD[scalar, row].tolist(),
                FIELD[scalar, row].tolist(),
            ]
            left = np.zeros((len(original_rows) + 2, d * d, d * d), dtype=np.uint8)
            for i, coefficient in enumerate(u):
                for j, a in enumerate(original_rows):
                    left[i + j] ^= np.array(
                        field.matmul(coefficient, a.tolist()), dtype=np.uint8
                    )
            reference = np.zeros((len(original_rows) + 4, d * d, d * d), dtype=np.uint8)
            for i, a in enumerate(left):
                for j, coefficient in enumerate(u):
                    reference[i + j] ^= np.array(
                        field.matmul(a.tolist(), coefficient), dtype=np.uint8
                    )
            actual = np.zeros_like(reference)
            actual[: len(dressed)] = dressed.transpose(0, 1, 2, 4, 3).reshape(
                len(dressed), d * d, d * d
            )
            assert np.array_equal(actual, reference)
            polynomial_controls += reference.size
        for z in range(4):
            original = (
                evaluate(selected_tensor, z).transpose(0, 1, 3, 2).reshape(d * d, d * d)
            )
            f = M4[scalar][M4[z][1 ^ z]]
            u = np.eye(d * d, dtype=np.uint8) ^ FIELD[f, row]
            expected_row = field.matmul(
                field.matmul(u.tolist(), original.tolist()), u.tolist()
            )
            actual = evaluate(dressed, z).transpose(0, 1, 3, 2).reshape(d * d, d * d)
            assert np.array_equal(actual, np.array(expected_row, dtype=np.uint8))
        dressings.append(
            {
                "length": n,
                "binary_width": width,
                "labels": d,
                "selected_indices": selected,
                "scalar": scalar,
                "original_degree": len(selected_tensor) - 1,
                "dressed_degree": len(dressed) - 1,
                "new_leading_rank_pair": cut_ranks(dressed[-1]),
            }
        )
        if n == 8:
            lex = b[:, :8, :8, :8, :8][-1].transpose(0, 1, 3, 2).reshape(64, 64)
            assert any(
                x for row in field.matmul(lex.tolist(), lex.tolist()) for x in row
            )
    assert dressings[0]["dressed_degree"] == 14 and dressings[1]["dressed_degree"] == 35
    files = (
        "tools/quantum_braid_reversal.py",
        "tools/check_quantum_braid_reversal.py",
        "tools/quantum_braid_path.py",
        "tools/quantum_block_braid.py",
        "tools/paired_path_crossing.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fixed_field_descent.py",
        "tools/check_quantum_braid_path_profiles.py",
        "tools/search_auxiliary_extension_span.py",
    )
    return {
        "status": "pass",
        "local_reflection_identities": local,
        "word_reflection_controls": words,
        "reversal_stable_binary_labels": labels,
        "full_leading_tensor_entries": full_leading_entries,
        "independent_polynomial_product_entries": polynomial_controls,
        "lexicographic_eight_label_square_zero_negative_control": True,
        "finite_dressings": dressings,
        "scope": "Uniform reversal alphabets and leading nilpotence; finite pole cancellation examples, without a favorable rank-pressure theorem.",
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
