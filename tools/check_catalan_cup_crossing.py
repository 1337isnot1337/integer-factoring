#!/usr/bin/env python3
"""Independent physical checks for Catalan cup frames and complete crossings."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from catalan_cup_crossing import (
    binary_length,
    cup_action,
    cup_basis,
    cup_frame,
    crossing_layout,
    crossing_tensor,
    field_kronecker,
    label_count,
    multiply,
    triangular_dual,
)
from check_fixed_field_descent import Field
from check_quantum_braid_fusion_seam import rank as f16_rank
from paired_path_crossing import power
from quantum_block_braid import braid_word
from quantum_braid_spectral_fusion import FIELD, matrix_product, physical_generator

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = Field(0b10011)


def binary_rank(matrix):
    packed = np.packbits(matrix, axis=1, bitorder="little")
    basis = {}
    for row in packed:
        value = int.from_bytes(row.tobytes(), "little")
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def rank(matrix, q):
    return binary_rank(matrix) if q == 1 else f16_rank(matrix)


def cut_ranks(tensor, q):
    d = tensor.shape[0]
    return [
        rank(tensor.reshape(d * d, d * d), q),
        rank(tensor.transpose(0, 2, 1, 3).reshape(d * d, d * d), q),
    ]


def literal_cup_entry(mate, word, q):
    result = 1
    for left, right in enumerate(mate):
        if left < right:
            a, b = (word >> left) & 1, (word >> right) & 1
            if a == b:
                return 0
            if a == 0:
                result = REFERENCE.mul(result, q)
    return result


def reference_inverse(matrix):
    n = len(matrix)
    a = [
        list(map(int, row)) + [int(i == j) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        inverse = REFERENCE.power(a[j][j], 14)
        a[j] = [REFERENCE.mul(inverse, v) for v in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                c = a[i][j]
                a[i] = [v ^ REFERENCE.mul(c, w) for v, w in zip(a[i], a[j])]
    return np.array([row[n:] for row in a], dtype=np.uint8)


def physical_polynomial(length, first, second, q):
    """Direct 2^(2n) binary-word evolution; no link-state routes or output map."""
    mates, dyck_words = cup_basis(length)
    minor = np.array(
        [[literal_cup_entry(p, w, q) for p in mates] for w in dyck_words],
        dtype=np.uint8,
    )
    decoder = reference_inverse(minor)
    f = np.array(
        [literal_cup_entry(mates[first], w, q) for w in range(1 << length)],
        dtype=np.uint8,
    )
    g = np.array(
        [literal_cup_entry(mates[second], w, q) for w in range(1 << length)],
        dtype=np.uint8,
    )
    states = FIELD[g[:, None], f[None, :]].reshape(-1, 1)
    words = np.arange(1 << (2 * length))
    for site in braid_word(length, length):
        bits = (words >> site) & 3
        mixed = (bits == 1) | (bits == 2)
        diagonal = np.where(bits == 1, power(q, -1), np.where(bits == 2, q, 0))
        following = np.zeros((len(words), states.shape[1] + 1), dtype=np.uint8)
        following[:, :-1] = FIELD[q, states]
        following[:, 1:] ^= FIELD[diagonal[:, None], states]
        following[words[mixed] ^ (3 << site), 1:] ^= states[mixed]
        states = following
    selected = [p | (r << length) for p in dyck_words for r in dyck_words]
    decoded = multiply(field_kronecker(decoder, decoder), states[selected], q)
    gamma = power(q, -length * length // 2)
    decoded = FIELD[gamma, decoded]
    strict = np.zeros((len(decoded), decoded.shape[1] + 1), dtype=np.uint8)
    strict[:, :-1] = FIELD[gamma, decoded]
    strict[:, 1:] ^= FIELD[1 ^ gamma, decoded]
    return strict


def retained_wall():
    mates, words = cup_basis(8)
    indices = {p: i for i, p in enumerate(mates)}
    es = []
    for site in range(3):
        e = np.zeros((len(mates), len(mates)), dtype=np.uint8)
        for j, p in enumerate(mates):
            target, scalar = cup_action(p, site, 8)
            e[indices[target], j] = scalar
        es.append(e)
    identity = np.eye(len(mates), dtype=np.uint8)
    delta = 8 ^ power(8, -1)
    f3 = (
        identity
        ^ es[0]
        ^ es[1]
        ^ FIELD[
            power(delta, -1),
            matrix_product(es[0], es[1]) ^ matrix_product(es[1], es[0]),
        ]
    )
    f4 = f3 ^ FIELD[delta, matrix_product(f3, matrix_product(es[2], f3))]
    frame = cup_frame(mates, range(256), 8)
    gram = matrix_product(frame.T, frame)
    output = crossing_layout(4, 8)[3]
    assert np.array_equal(matrix_product(f4, f4), f4)
    assert f16_rank(f4) == 1 and not np.any(matrix_product(gram, f4))
    assert np.any(matrix_product(output, f4))
    return {
        "physical_strands": 8,
        "primal_rank": 14,
        "gram_rank": f16_rank(gram),
        "retained_f4_rank": 1,
        "gram_annihilates_f4": True,
        "coordinate_product_decoder_sees_f4": True,
    }


def check():
    frame_controls, profiles = [], []
    local_entries = polynomial_entries = rank_pairs = 0
    for q in (1, 8):
        for n in (2, 4, 6, 8, 10):
            mates, words = cup_basis(n)
            assert len(mates) == label_count(n)
            minor = cup_frame(mates, words, q)
            inverse = triangular_dual(minor)
            assert np.array_equal(inverse, reference_inverse(minor))
            assert np.array_equal(multiply(inverse, minor, q), np.eye(len(mates)))
            full = cup_frame(mates, range(1 << n), q)
            assert rank(full, q) == len(mates)
            gram = multiply(full.T, full, q)
            if q == 1:
                assert not np.any(gram)
            frame_controls.append(
                {
                    "q": q,
                    "length": n,
                    "labels": len(mates),
                    "primal_rank": len(mates),
                    "gram_rank": rank(gram, q),
                }
            )
            if n <= 8:
                index = {p: i for i, p in enumerate(mates)}
                for site in range(n - 1):
                    actual = multiply(physical_generator(n, site, q), full, q)
                    expected = np.zeros_like(full)
                    for j, p in enumerate(mates):
                        target, scalar = cup_action(p, site, q)
                        expected[:, j] = FIELD[scalar, full[:, index[target]]]
                    assert np.array_equal(actual, expected)
                    local_entries += actual.size
        for n in (4, 6, 8):
            tensor = crossing_tensor(n, q)
            d = tensor.shape[1]
            identity = np.fromfunction(
                lambda p, r, s, t: (p == s) & (r == t), (d, d, d, d), dtype=int
            ).astype(np.uint8)
            swap = identity.transpose(0, 1, 3, 2)
            assert np.array_equal(tensor[0], identity)
            assert np.array_equal(np.bitwise_xor.reduce(tensor, axis=0), swap)
            if q == 1:
                assert np.all(tensor <= 1)
            for first, second in ((0, d - 1), (d - 1, 0)):
                reference = physical_polynomial(n, first, second, q)
                compiled = tensor[:, first, second].reshape(len(tensor), d * d).T
                assert np.array_equal(compiled, reference[:, : len(tensor)])
                assert not np.any(reference[:, len(tensor) :])
                polynomial_entries += reference.size
            pairs = [cut_ranks(a, q) for a in tensor]
            if q == 1:
                for j in (0, len(tensor) // 2, len(tensor) - 1):
                    a = tensor[j].reshape(d * d, d * d)
                    assert binary_rank(a) == f16_rank(a)
            rank_pairs += len(pairs)
            profiles.append(
                {
                    "q": q,
                    "length": n,
                    "labels": d,
                    "polynomial_degree": len(tensor) - 1,
                    "strict_endpoints": True,
                    "coefficient_span_dimension": rank(
                        tensor.reshape(len(tensor), -1), q
                    ),
                    "complete_coefficient_rank_pairs": pairs,
                }
            )
    lengths = [binary_length(width) for width in range(25)]
    for width, n in enumerate(lengths):
        assert 1 << width <= label_count(n) < 4 * (1 << width)
        assert n == 2 or label_count(n - 2) < 1 << width
    wall = retained_wall()
    sources = (
        "research/catalan-cup-crossings.md",
        "tools/catalan_cup_crossing.py",
        "tools/check_catalan_cup_crossing.py",
        "tools/quantum_braid_spectral_fusion.py",
        "tools/quantum_block_braid.py",
        "tools/paired_path_crossing.py",
        "tools/check_quantum_braid_fusion_seam.py",
        "tools/quantum_braid_fusion_seam.py",
        "tools/quantum_braid_nilpotent_seam.py",
        "tools/quantum_braid_path.py",
        "tools/check_fixed_field_descent.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
    )
    return {
        "status": "pass",
        "frame_controls": frame_controls,
        "independent_physical_local_action_entries": local_entries,
        "independent_physical_polynomial_entries": polynomial_entries,
        "complete_coefficient_rank_pairs_checked": rank_pairs,
        "binary_selected_lengths": lengths,
        "retained_negligible_wall_control": wall,
        "complete_profiles": profiles,
        "scope": "Uniform Catalan primal frames and coordinate duals with strict crossings; no favorable complete rank sum or stronger factoring score.",
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
