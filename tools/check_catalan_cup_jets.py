#!/usr/bin/env python3
"""Check full cup Hasse jets against independent truncated-ring physical words."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from catalan_cup_crossing import binary_product, crossing_tensor, cup_basis
from catalan_cup_jets import crossing_jets, cup_series, inverse_series
from check_catalan_cup_crossing import binary_rank
from quantum_block_braid import braid_word

ROOT = Path(__file__).resolve().parents[1]


class ReferenceRing:
    """Small finite-control table only; the uniform compiler has no such table."""

    def __init__(self, order):
        self.order = order
        size = 1 << order
        self.table = np.zeros((size, size), dtype=np.uint8)
        for a in range(size):
            for b in range(size):
                value = 0
                for j in range(order):
                    if b >> j & 1:
                        value ^= a << j
                self.table[a, b] = value & (size - 1)

    def power_q(self, exponent):
        value = 3 if exponent >= 0 else (1 << self.order) - 1
        result = 1
        for _ in range(abs(exponent)):
            result = int(self.table[result, value])
        return result

    def cup(self, mate, word):
        value = 1
        for left, right in enumerate(mate):
            if left < right:
                a, b = word >> left & 1, word >> right & 1
                if a == b:
                    return 0
                if a == 0:
                    value = int(self.table[value, 3])
        return value

    def inverse(self, minor):
        n = len(minor)
        augmented = [
            list(map(int, row)) + [int(i == j) for j in range(n)]
            for i, row in enumerate(minor)
        ]
        for j in range(n):
            assert augmented[j][j] == 1
            for i in range(n):
                if i != j and augmented[i][j]:
                    c = augmented[i][j]
                    augmented[i] = [
                        int(v ^ self.table[c, w])
                        for v, w in zip(augmented[i], augmented[j])
                    ]
        assert [row[:n] for row in augmented] == np.eye(n, dtype=int).tolist()
        return np.array([row[n:] for row in augmented], dtype=np.uint8)

    def matrix_product(self, a, b):
        return np.bitwise_xor.reduce(self.table[a[:, :, None], b[None, :, :]], axis=1)

    def physical_polynomial(self, n, first, second):
        mates, dyck = cup_basis(n)
        minor = [[self.cup(p, w) for p in mates] for w in dyck]
        decoder = self.inverse(minor)
        f = np.array([self.cup(mates[first], w) for w in range(1 << n)], dtype=np.uint8)
        g = np.array(
            [self.cup(mates[second], w) for w in range(1 << n)], dtype=np.uint8
        )
        states = self.table[g[:, None], f[None, :]].reshape(-1, 1)
        words = np.arange(1 << (2 * n))
        inverse_q = (1 << self.order) - 1
        for site in braid_word(n, n):
            bits = words >> site & 3
            mixed = (bits == 1) | (bits == 2)
            diagonal = np.where(bits == 1, inverse_q, np.where(bits == 2, 3, 0))
            following = np.zeros((len(words), states.shape[1] + 1), dtype=np.uint8)
            following[:, :-1] = self.table[3, states]
            following[:, 1:] ^= self.table[diagonal[:, None], states]
            following[words[mixed] ^ (3 << site), 1:] ^= states[mixed]
            states = following
        kronecker = self.table[decoder[:, None, :, None], decoder[None, :, None, :]]
        kronecker = kronecker.reshape(len(mates) ** 2, len(mates) ** 2)
        selected = [p | (r << n) for p in dyck for r in dyck]
        decoded = self.matrix_product(kronecker, states[selected])
        gamma = self.power_q(-n * n // 2)
        decoded = self.table[gamma, decoded]
        strict = np.zeros((len(decoded), decoded.shape[1] + 1), dtype=np.uint8)
        strict[:, :-1] = self.table[gamma, decoded]
        strict[:, 1:] ^= self.table[1 ^ gamma, decoded]
        return strict


def profile(tensor):
    d = tensor.shape[2]
    rows, all_pieces, good_pieces = [], [], []
    for layer in tensor:
        ranks = []
        for a in layer:
            r = binary_rank(a.reshape(d * d, d * d))
            c = binary_rank(a.transpose(0, 2, 1, 3).reshape(d * d, d * d))
            ranks.append([r, c])
            if np.any(a):
                pieces = (a.reshape(-1), a.transpose(0, 2, 1, 3).reshape(-1))
                all_pieces.extend(pieces)
                if (r * c) ** 5 < d**14:
                    good_pieces.extend(pieces)
        rows.append(ranks)
    all_pieces, good_pieces = np.array(all_pieces), np.array(good_pieces)
    dimension = binary_rank(good_pieces)
    swap = np.bitwise_xor.reduce(tensor[0], axis=0).reshape(1, -1)
    member = binary_rank(np.concatenate((good_pieces, swap))) == dimension
    return {
        "all_coefficient_rank_pairs": rows,
        "coefficient_span_dimension": binary_rank(tensor.reshape(-1, d**4)),
        "span_with_cut_reflections": binary_rank(all_pieces),
        "exact_filter": "(r*c)^5 < L^14; gamma=7/5",
        "filtered_pieces": len(good_pieces),
        "filtered_span_dimension": dimension,
        "crossing_in_filtered_span": member,
    }


def check():
    profiles = []
    physical_entries = coefficient_pairs = inverse_controls = 0
    for n, order in ((4, 2), (4, 4), (4, 8), (6, 2), (6, 4), (6, 8), (8, 4)):
        tensor = crossing_jets(n, order)
        d = tensor.shape[2]
        baseline = crossing_tensor(n, 1)
        assert np.array_equal(tensor[0, : len(baseline)], baseline)
        assert not np.any(tensor[0, len(baseline) :])
        identity = np.fromfunction(
            lambda p, q, r, s: (p == r) & (q == s), (d, d, d, d), dtype=int
        ).astype(np.uint8)
        swap = identity.transpose(0, 1, 3, 2)
        assert np.array_equal(tensor[0, 0], identity) and not np.any(tensor[1:, 0])
        endpoints = np.bitwise_xor.reduce(tensor, axis=1)
        assert np.array_equal(endpoints[0], swap) and not np.any(endpoints[1:])
        mates, words = cup_basis(n)
        minor = cup_series(mates, words, order)
        inverse = inverse_series(minor)
        for degree in range(order):
            product = np.zeros((d, d), dtype=np.uint8)
            for j in range(degree + 1):
                product ^= binary_product(minor[j], inverse[degree - j])
            assert np.array_equal(
                product, np.eye(d) if degree == 0 else np.zeros((d, d))
            )
            inverse_controls += 1
        if (n, order) in ((4, 8), (6, 8), (8, 4)):
            reference = ReferenceRing(order).physical_polynomial(n, 0, d - 1)
            compiled = np.bitwise_or.reduce(
                np.array(
                    [
                        tensor[j, :, 0, d - 1].reshape(tensor.shape[1], d * d).T << j
                        for j in range(order)
                    ]
                ),
                axis=0,
            )
            assert np.array_equal(compiled, reference[:, : tensor.shape[1]])
            assert not np.any(reference[:, tensor.shape[1] :])
            physical_entries += order * reference.size
        record = {
            "length": n,
            "labels": d,
            "jet_order": order,
            "t_degree": tensor.shape[1] - 1,
            "strict_jet_endpoints": True,
            **profile(tensor),
        }
        if n == 6:
            assert not record["crossing_in_filtered_span"]
        profiles.append(record)
        coefficient_pairs += order * tensor.shape[1]
    sources = (
        "research/catalan-cup-jets.md",
        "tools/catalan_cup_jets.py",
        "tools/check_catalan_cup_jets.py",
        "tools/catalan_cup_crossing.py",
        "tools/check_catalan_cup_crossing.py",
        "tools/quantum_block_braid.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_braid_spectral_fusion.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
        "tools/check_fixed_field_descent.py",
    )
    return {
        "status": "pass",
        "independent_physical_binary_jet_entries": physical_entries,
        "inverse_series_coefficients_checked": inverse_controls,
        "complete_coefficient_rank_pairs_checked": coefficient_pairs,
        "complete_profiles": profiles,
        "scope": "Uniform effective Hasse-jet enlargement of a strict F2 crossing; no favorable complete decomposition or stronger factoring bound.",
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
