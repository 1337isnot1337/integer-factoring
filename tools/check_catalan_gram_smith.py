#!/usr/bin/env python3
"""Verify finite Smith profiles of full Catalan cup Gram matrices."""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from catalan_cup_crossing import cup_basis, cup_frame, triangular_dual
from catalan_gram_smith import loop_count, smith_valuation_profile
from paired_path_crossing import power
from quantum_braid_spectral_fusion import FIELD, matrix_product

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    2: {1: 1},
    4: {1: 2},
    6: {1: 4, 3: 1},
    8: {1: 8, 3: 6},
    10: {1: 16, 3: 26},
    12: {1: 32, 3: 100},
    14: {1: 64, 3: 364, 7: 1},
    16: {1: 128, 3: 1288, 7: 14},
}


@lru_cache(maxsize=None)
def quantum_integer(q, index):
    if index == 0:
        return 0
    delta = q ^ power(q, -1)
    assert delta
    return FIELD[power(q, index) ^ power(q, -index), power(delta, -1)]


@lru_cache(maxsize=None)
def split_block(q, height, magnetic):
    return (
        (1, quantum_integer(q, magnetic)),
        (
            power(q, -height + magnetic - 1),
            FIELD[
                power(q, magnetic),
                quantum_integer(q, height - magnetic + 1),
            ],
        ),
    )


def dyck_paths(length):
    paths = []

    def visit(height, path):
        if len(path) == length:
            if height == 0:
                paths.append(tuple(path))
            return
        if height < length // 2:
            visit(height + 1, path + [1])
        if height:
            visit(height - 1, path + [-1])

    visit(0, [])
    return tuple(paths)


def path_entry(path, word, q):
    height = magnetic = 0
    value = 1
    for position, direction in enumerate(path):
        bit = word >> position & 1
        sector = magnetic + bit
        following = sector - int(direction == -1)
        if not 0 <= following <= height + direction:
            return 0
        scalar = 1
        if 0 < sector < height + 1:
            scalar = split_block(q, height, sector)[bit][int(direction == -1)]
        value = FIELD[value, scalar]
        magnetic = following
        height += direction
    return value if magnetic == 0 else 0


def path_norm(path, q):
    height = 0
    value = 1
    for direction in path:
        if direction == -1:
            value = FIELD[
                value,
                FIELD[
                    q,
                    FIELD[
                        quantum_integer(q, height),
                        quantum_integer(q, height + 1),
                    ],
                ],
            ]
        height += direction
    assert height == 0 and value
    return value


def physical_gram_controls():
    checked = 0
    q = 8
    delta = q ^ power(q, -1)
    for length in (2, 4, 6, 8):
        matchings, _ = cup_basis(length)
        frame = cup_frame(matchings, range(1 << length), q)
        actual = matrix_product(frame.T, frame)
        cups = length // 2
        expected = [
            [
                FIELD[
                    power(q, cups),
                    power(delta, loop_count(first, second)),
                ]
                for second in matchings
            ]
            for first in matchings
        ]
        assert actual.tolist() == expected
        checked += actual.size
    return checked


def unrestricted_path_controls():
    q = 3
    delta = q ^ power(q, -1)
    assert delta
    controls = []
    for length in (2, 4, 6, 8):
        mates, dyck_words = cup_basis(length)
        paths = dyck_paths(length)
        assert len(paths) == len(mates)
        vectors = np.array(
            [
                [path_entry(path, word, q) for path in paths]
                for word in range(1 << length)
            ],
            dtype=np.uint8,
        )
        norms = [path_norm(path, q) for path in paths]
        gram = matrix_product(vectors.T, vectors)
        expected = np.zeros_like(gram)
        np.fill_diagonal(expected, norms)
        assert np.array_equal(gram, expected)

        physical_frame = cup_frame(mates, range(1 << length), q)
        coordinate_dual = triangular_dual(cup_frame(mates, dyck_words, q))
        coordinates = matrix_product(coordinate_dual, vectors[list(dyck_words)])
        assert np.array_equal(matrix_product(physical_frame, coordinates), vectors)
        controls.append(
            {
                "length": length,
                "labels": len(paths),
                "path_word_entries": (1 << length) * len(paths),
                "orthogonal_nonzero_norms": len(norms),
                "regular_cup_coordinates": True,
            }
        )
    return controls


def check():
    physical_entries = physical_gram_controls()
    path_controls = unrestricted_path_controls()
    profiles = []
    for length, expected in EXPECTED.items():
        profile = smith_valuation_profile(length, precision=8)
        assert profile["gram_exponent_counts"] == expected
        assert profile["unresolved_at_least"] == 0
        profiles.append(profile)
    sources = (
        "research/catalan-gram-dual-smith.md",
        "tools/catalan_gram_smith.py",
        "tools/check_catalan_gram_smith.py",
        "tools/catalan_cup_crossing.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_braid_spectral_fusion.py",
    )
    return {
        "status": "pass",
        "physical_loop_gram_entries_checked_over_F16": physical_entries,
        "unrestricted_q_path_basis_controls_over_F16": path_controls,
        "truncated_toeplitz_profiles": profiles,
        "scope": (
            "Exact Smith valuations through 16 strands and unrestricted path "
            "controls through 8 strands; no favorable rank sum or factoring gain."
        ),
        "source_sha256": {
            path: sha256((ROOT / path).read_bytes()).hexdigest()
            for path in sources
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
