#!/usr/bin/env python3
"""Local radical, path-action, metric and F4 descent controls for the braid compiler."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from paired_path_crossing import MULTIPLY as M16, path_entry as ambient_frame, power
from quantum_block_braid import polynomial_entry as ambient_entry
from quantum_braid_path import (
    EMBED,
    MULTIPLY as M4,
    coefficient_entry,
    crossing_weights,
    local_action,
    polynomial_entry,
)

ROOT = Path(__file__).resolve().parents[1]
F16, F4 = Field(0b10011), Field(0b111)
Q = 8


def qi(k):
    return M16[power(Q, k) ^ power(Q, -k)][power(Q ^ power(Q, -1), -1)]


def transpose(a):
    return list(map(list, zip(*a)))


def add(a, b):
    return [[x ^ y for x, y in zip(u, v)] for u, v in zip(a, b)]


def scale(a, c, multiply):
    return [[multiply[x][c] for x in row] for row in a]


def symmetrizer(h):
    i = [[0] * (h + 1) for _ in range(1 << h)]
    for w in range(1 << h):
        inversions = sum(
            j - (w & ((1 << j) - 1)).bit_count() for j in range(h) if w >> j & 1
        )
        i[w][w.bit_count()] = power(Q, -inversions)
    gram = F16.matmul(transpose(i), i)
    inverse_gram = [
        [power(gram[r][r], -1) if r == s else 0 for s in range(h + 1)]
        for r in range(h + 1)
    ]
    assert all(gram[r][r] for r in range(h + 1))
    return i, gram, F16.matmul(F16.matmul(i, inverse_gram), transpose(i))


def primal(path, word, final_magnetic):
    h, r, value = 0, 0, 1
    for j, direction in enumerate(path):
        bit = word >> j & 1
        nr = r + bit - int(direction == -1)
        if not 0 <= nr <= h + direction:
            return 0
        scalar = (
            (1 if bit == 0 else power(Q, -h + r))
            if direction == 1
            else (qi(r) if bit == 0 else M16[power(Q, r + 1)][qi(h - r)])
        )
        value = M16[value][scalar]
        h, r = h + direction, nr
    return value if r == final_magnetic else 0


def generator_vector(vector, site):
    result = [0] * len(vector)
    for w, value in enumerate(vector):
        bits = w >> site & 3
        if bits in (1, 2):
            result[w] ^= M16[value][power(Q, -1) if bits == 1 else Q]
            result[w ^ (3 << site)] ^= value
    return result


def local_path_matrix(start, length, site):
    candidates = []
    for steps in product((1, -1), repeat=length):
        h = start
        for step in steps:
            h += step
            if not 0 <= h <= 3:
                break
        else:
            candidates.append((1,) * start + steps)
    indices = {p: i for i, p in enumerate(candidates)}
    result = [[0] * len(candidates) for _ in candidates]
    for j, p in enumerate(candidates):
        for q, c in local_action(p, start + site):
            result[indices[q]][j] ^= c
    return candidates, result


def check():
    for a, b in product(range(4), repeat=2):
        assert M16[EMBED[a]][EMBED[b]] == EMBED[M4[a][b]]
    assert Q ^ power(Q, -1) == EMBED[2] and qi(5) == 0
    syms = [symmetrizer(h) for h in range(5)]
    metric_entries = 0
    for h in range(4):
        u = [[0] * (2 * h + 2) for _ in range(2 * h + 2)]
        for r in range(h + 1):
            u[r][r] = 1
            u[h + 1 + r][r + 1] = power(Q, -h + r)
        for m in range(1, h + 1):
            u[m][h + 1 + m] = qi(m)
            u[h + m][h + 1 + m] = M16[power(Q, m)][qi(h - m + 1)]
        g = syms[h][1]
        old = [
            [
                g[r % (h + 1)][s % (h + 1)] if r // (h + 1) == s // (h + 1) else 0
                for s in range(2 * h + 2)
            ]
            for r in range(2 * h + 2)
        ]
        values = [syms[h + 1][1][r][r] for r in range(h + 2)]
        if h:
            down_scale = M16[Q][M16[qi(h)][qi(h + 1)]]
            values += [M16[down_scale][syms[h - 1][1][r][r]] for r in range(h)]
        expected = [
            [values[r] if r == s else 0 for s in range(2 * h + 2)]
            for r in range(2 * h + 2)
        ]
        assert F16.matmul(F16.matmul(transpose(u), old), u) == expected
        metric_entries += (2 * h + 2) ** 2
    _, _, f4 = syms[4]
    partial_trace = [
        [M16[Q][f4[i][j]] ^ M16[power(Q, -1)][f4[i + 8][j + 8]] for j in range(8)]
        for i in range(8)
    ]
    assert not any(x for row in partial_trace for x in row)
    assert F16.rank(f4) == 5
    # Exact local recoupling, including the nonzero negligible height-four branch.
    recoupling_entries = 0
    for h in range(4):
        branches = (1,) if h == 0 else (1, -1)
        for r in range(h + 1):
            vectors = {
                b: [primal((1,) * h + (b, -b), w, r) for w in range(1 << (h + 2))]
                for b in branches
            }
            for branch, vector in vectors.items():
                c = M16[qi(h + 2 if branch == 1 else h)][power(qi(h + 1), -1)]
                total = [0] * len(vector)
                for v in vectors.values():
                    total = [x ^ M16[c][y] for x, y in zip(total, v)]
                assert generator_vector(vector, h) == total
                recoupling_entries += len(vector)
    local_relations = 0
    for h in range(4):
        _, a = local_path_matrix(h, 3, 0)
        _, b = local_path_matrix(h, 3, 1)
        assert F4.matmul(a, a) == scale(a, 2, M4)
        assert F4.matmul(b, b) == scale(b, 2, M4)
        assert F4.matmul(F4.matmul(a, b), a) == a
        assert F4.matmul(F4.matmul(b, a), b) == b
        local_relations += 4
        candidates, a = local_path_matrix(h, 4, 0)
        _, b = local_path_matrix(h, 4, 1)
        _, c = local_path_matrix(h, 4, 2)
        identity = [
            [int(i == j) for j in range(len(candidates))]
            for i in range(len(candidates))
        ]
        f2 = add(identity, scale(a, 3, M4))
        f3 = add(f2, F4.matmul(F4.matmul(f2, b), f2))
        f4_path = add(f3, scale(F4.matmul(F4.matmul(f3, c), f3), 2, M4))
        assert not any(x for row in f4_path for x in row)
        local_relations += 1
    norm_words = 0
    for n in range(2, 11, 2):
        for p in paths(n):
            h, nu = 0, 1
            for step in p:
                if step == -1:
                    nu = M16[nu][M16[Q][M16[qi(h)][qi(h + 1)]]]
                h += step
            for word in range(1 << n):
                assert (
                    ambient_frame(p, word, Q, True)
                    == M16[ambient_frame(p, word, Q)][power(nu, -1)]
                )
                norm_words += 1
    entries, reference_entries = 0, 0
    for n in (2, 4, 6):
        ps = paths(n)
        for ids in product(range(len(ps)), repeat=4):
            four = tuple(ps[i] for i in ids)
            values = polynomial_entry(four)
            total = 0
            for value in values:
                total ^= value
            assert values[0] == int(ids[0] == ids[2] and ids[1] == ids[3])
            assert total == int(ids[0] == ids[3] and ids[1] == ids[2])
            entries += 1
            if n <= 4 or ids in ((1, 3, 2, 4), (4, 3, 1, 2), (0, 4, 4, 0)):
                a = coefficient_entry(four)
                old = [M16[power(Q, n * n // 2 - j)][EMBED[x]] for j, x in enumerate(a)]
                delta = power(Q, -n * n // 2)
                strict = [0] * (n * n + 2)
                for j, value in enumerate(old):
                    strict[j] ^= M16[delta][value]
                    strict[j + 1] ^= M16[1 ^ delta][value]
                assert strict == ambient_entry(four)
                reference_entries += 1
    coordinates = {
        EMBED[a] ^ M16[Q][EMBED[b]]: (a, b) for a, b in product(range(4), repeat=2)
    }
    assert len(coordinates) == 16
    phase_controls = 0
    for n in range(2, 41, 2):
        weights = crossing_weights(n)
        phase = n * n // 2 % 5
        tau = (0, 1, 0, 2, 3)[phase]
        for j, w in enumerate(weights):
            a, b = coordinates[power(Q, n * n // 2 - j)]
            assert w == a ^ M4[tau][b]
            phase_controls += 1
        assert weights[0] == 1
    assert coordinates[Q][0] == 0 and coordinates[power(Q, 2)][0] == 1
    files = (
        "tools/quantum_braid_path.py",
        "tools/check_quantum_braid_path.py",
        "tools/quantum_block_braid.py",
        "tools/paired_path_crossing.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fixed_field_descent.py",
        "tools/check_fourier_multiplication_gates.py",
    )
    return {
        "status": "pass",
        "field": "F4",
        "metric_identity_entries": metric_entries,
        "height_four_projector_rank": 5,
        "height_four_partial_quantum_trace_zero": True,
        "exact_recoupling_entries": recoupling_entries,
        "local_TL_and_wall_relations": local_relations,
        "norm_dual_word_controls": norm_words,
        "complete_small_endpoint_entries": entries,
        "ambient_reference_polynomials": reference_entries,
        "independent_phase_projection_controls": phase_controls,
        "scope": "Local quotient and phase descent identities support an all-length F4 path compiler; no favorable uniform two-cut rank bound.",
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
