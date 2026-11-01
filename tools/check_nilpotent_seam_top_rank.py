#!/usr/bin/env python3
"""Certify the closed leading fusion automaton and its all-length rank series."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_fixed_field_descent import Field
from check_quantum_braid_fusion_seam import cut_ranks
from nilpotent_seam_rank_bounds import central_weights
from nilpotent_seam_top_rank import four_block_rank, top_auxiliary_rank
from quantum_braid_graded_seam import fusion_rank
from quantum_braid_nilpotent_seam import (
    FIELD,
    IDENTITY,
    N,
    Q,
    coefficient_tensor,
    matrix_product,
)
from quantum_braid_path import EMBED

ROOT = Path(__file__).resolve().parents[1]
F4 = Field(0b111)
X2 = Q ^ IDENTITY
X0 = X2 ^ FIELD[2, IDENTITY ^ N]
Y0 = FIELD[2, N]
GATES = ((IDENTITY, X0, Y0), (IDENTITY, X2, np.zeros_like(N)))


def literal_product(word):
    values = [IDENTITY.copy()]
    for bit in word:
        following = [np.zeros_like(N) for _ in range(len(values) + 2)]
        for i, a in enumerate(GATES[bit]):
            for j, b in enumerate(values):
                following[i + j] ^= matrix_product(a, b)
        values = following
    return values


def normalize(a):
    if not np.any(a):
        return None
    first = int(a.flat[np.flatnonzero(a)[0]])
    return FIELD[F4.power(first, 2), a].tobytes()


def decode(key):
    return np.frombuffer(key, dtype=np.uint8).reshape(4, 4)


def automaton():
    generators = []
    local_controls = []
    for word in product((0, 1), repeat=4):
        values = literal_product(word)
        assert not np.any(values[8])
        key = normalize(values[7])
        local_controls.append({"word": list(word), "top_matrix": values[7].tolist()})
        if key is not None:
            assert sum(word) <= 1
            generators.append((key, sum(word)))
    assert len(generators) == 5
    keys, indices, transitions = [normalize(IDENTITY)], {normalize(IDENTITY): 0}, []
    for key in keys:
        row = []
        for g, _ in generators:
            target = normalize(matrix_product(decode(g), decode(key)))
            if target is not None and target not in indices:
                indices[target] = len(keys)
                keys.append(target)
            row.append(None if target is None else indices[target])
        transitions.append(row)
    assert len(keys) == len(transitions) == 22
    # Recheck every transition against the complete closed state set.
    for i, key in enumerate(keys):
        for j, (g, _) in enumerate(generators):
            target = normalize(matrix_product(decode(g), decode(key)))
            assert transitions[i][j] == (None if target is None else indices[target])
    return generators, keys, transitions, local_controls


def transfer(vector, generators, transitions):
    result = [0] * len(vector)
    for i, row in enumerate(transitions):
        v, c = vector[2 * i : 2 * i + 2]
        for (_, charge), target in zip(generators, row):
            if target is not None:
                result[2 * target] += c if charge else v
                result[2 * target + 1] += v + c if charge else c
    return result


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def check():
    generators, keys, transitions, controls = automaton()
    n2 = matrix_product(N, N)
    outputs = [[0] * 44 for _ in range(2)]
    for i, key in enumerate(keys):
        outputs[0][2 * i] = F4.rank(decode(key).tolist())
        outputs[1][2 * i] = F4.rank(matrix_product(n2, decode(key)).tolist())
    initial = [1, 0] + [0] * 42
    vectors = [initial]
    for _ in range(64):
        vectors.append(transfer(vectors[-1], generators, transitions))
    series = [[dot(o, v) for v in vectors] for o in outputs]
    assert series[0] == [four_block_rank(m) for m in range(65)]
    assert series[1] == [top_auxiliary_rank(4 * m + 2) for m in range(65)]
    denominators = ([1, -7, 17, -17, 7, -1], [1, -3, 1])
    numerators = ([4, -24, 55, -52, 22], [2, -4])
    for o, series_values, den, num in zip(outputs, series, denominators, numerators):
        degree = len(den) - 1
        residual = [
            sum(c * vectors[degree - j][i] for j, c in enumerate(den))
            for i in range(44)
        ]
        for _ in range(44):
            assert dot(o, residual) == 0
            residual = transfer(residual, generators, transitions)
        convolution = [
            sum(den[j] * series_values[i - j] for j in range(i + 1))
            for i in range(degree)
        ]
        assert convolution == list(num)
    # Independent full word polynomials, before any four-letter grouping.
    word_controls = []
    for n in (4, 6, 8):
        a = 2 * n - n // 4
        total = 0
        for word in product((0, 1), repeat=n):
            values = literal_product(word)
            assert all(not np.any(x) for x in values[a + 1 :])
            total += fusion_rank(sum(word)) * F4.rank(values[a].tolist())
        assert total == top_auxiliary_rank(n)
        word_controls.append({"length": n, "literal_words": 1 << n, "rank": total})
    # Refined whole-coefficient upper bounds; equality is deliberately not assumed.
    profiles = []
    for n in (4, 6, 8):
        top = np.array(EMBED, dtype=np.uint8)[coefficient_tensor(n)[-1]]
        pairs = [cut_ranks(a) for a in top]
        bound = top_auxiliary_rank(n)
        assert bound <= 4 * central_weights(n)[2 * n - n // 4]
        for b, (r, _) in enumerate(pairs):
            assert r <= min(b + 1, n - b + 1) * bound
        profiles.append({"length": n, "central_rank": bound, "rank_pairs": pairs})
    assert profiles[1]["central_rank"] == 2
    assert all(r == c == 0 for r, c in profiles[1]["rank_pairs"])
    sources = (
        "research/nilpotent-leading-rank-automaton.md",
        "tools/nilpotent_seam_top_rank.py",
        "tools/check_nilpotent_seam_top_rank.py",
        "tools/nilpotent_seam_rank_bounds.py",
        "tools/quantum_braid_nilpotent_seam.py",
        "tools/quantum_braid_graded_seam.py",
        "tools/quantum_braid_fusion_seam.py",
        "tools/check_quantum_braid_fusion_seam.py",
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
        "projective_states": len(keys),
        "integer_counters": 44,
        "complete_closure_transitions": len(keys) * len(generators),
        "observable_krylov_controls": 88,
        "local_word_controls": controls,
        "generators": [
            {"matrix": decode(k).tolist(), "triplets": c} for k, c in generators
        ],
        "states": [decode(k).tolist() for k in keys],
        "transitions": transitions,
        "rank_outputs": outputs,
        "series_denominators": denominators,
        "series_numerators": numerators,
        "first_rank_values": [s[:17] for s in series],
        "literal_word_controls": word_controls,
        "complete_leading_coefficient_profiles": profiles,
        "scope": "Exact uniform leading central auxiliary ranks; contracted ranks may vanish; no favorable complete crossing rank sum or stronger factoring score.",
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
