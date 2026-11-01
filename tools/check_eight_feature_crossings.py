#!/usr/bin/env python3
"""Verify full local parity matrix products using independent four-leg tensors."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import random
import subprocess
from tempfile import TemporaryDirectory

from check_markov_pairing_transfer import matrix_product, ranks_of_coefficients


def entry(t, labels, matrices, boundary, off_shift=0):
    x, y, u, v = labels
    if x ^ y ^ u ^ v:
        return 0
    value = [1, 0, 0, 1]
    for i in range(t):
        xi, gi, hi = (x >> i) & 1, ((x ^ y) >> i) & 1, ((x ^ u) >> i) & 1
        b = matrices[4 * xi + 2 * gi + hi]
        factor = [int(gi == hi and r == c) ^ ((3 << (0 if r == c else off_shift)) if b >> (2 * r + c) & 1 else 0)
                  for r in range(2) for c in range(2)]
        value = matrix_product(value, factor)
    return value[0] if boundary == 0 else value[1] ^ value[2] ^ value[3]


def direct(t, matrices, boundary, off_shift=0):
    dimension = 1 << t
    tensors = {}
    for labels in product(range(dimension), repeat=4):
        x, y, u, v = labels
        polynomial = entry(t, labels, matrices, boundary, off_shift)
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        index = ((x * dimension + y) * dimension + u) * dimension + v
        for degree in range(t * (1 + off_shift) + 1):
            if polynomial >> degree & 1:
                tensors[degree] = tensors.get(degree, 0) | (1 << index)
    return tensors


def verify_short_witnesses(native):
    matrices, boundary, shift = native["matrices"], native["boundary"], native["off_shift"]
    def value(word):
        labels = [0] * 4
        for i, q in enumerate(word):
            x, g, h = q >> 2, (q >> 1) & 1, q & 1
            for j, bit in enumerate((x, x ^ g, x ^ h, x ^ g ^ h)):
                labels[j] |= bit << i
        return entry(len(word), labels, matrices, boundary, shift)
    if word := native["bit_flip_witness"]:
        assert value(word) != value([q ^ 4 for q in word])
    if word := native["order_witness"]:
        swapped = word.copy()
        i = native["order_swap"]
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        assert value(word) != value(swapped)
    if word := native["off_crossing_witness"]:
        assert any(((q >> 1) & 1) != (q & 1) for q in word)
        assert value(word)
    for q in native["identity_letters"]:
        assert ((q >> 1) & 1) == (q & 1) and matrices[q] == 0
    for q in native["paired_identity_letters"]:
        assert q in native["identity_letters"] and matrices[q ^ 4] == 0


def profile(binary, directory, t, matrices, boundary, off_shift=None):
    output = directory / "current.json"
    command = [str(binary), "profile", str(t), str(output), str(boundary), *map(str, matrices)]
    if off_shift is not None:
        command.append(str(off_shift))
    subprocess.run(command,
                   check=True, stdout=subprocess.DEVNULL)
    native = json.loads(output.read_text())
    if "off_shift" in native:
        verify_short_witnesses(native)
    return native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    rng = random.Random(20261030)
    witness = [0, 2, 4, 0, 1, 2, 4, 0]
    cases = [[0] * 8, [9, 0, 9, 0] * 2, witness]
    cases += [[rng.randrange(16) for _ in range(8)] for _ in range(9)]
    with TemporaryDirectory(prefix="eight-feature-ranks-") as temporary:
        directory = Path(temporary)
        full, circulant = directory / "full", directory / "circulant"
        for source, binary in (("eight_feature_crossings.cpp", full),
                               ("parity_matrix_crossings.cpp", circulant)):
            subprocess.run(["g++", "-O2", "-std=c++17", str(root / source), "-o", str(binary)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        count = 0
        for boundary in (0, 1):
            for matrices in cases:
                native = profile(full, directory, 3, matrices, boundary)
                tensors = direct(3, matrices, boundary)
                rows, columns = ranks_of_coefficients(3, tensors, 4)
                assert rows == native["row_ranks"] and columns == native["column_ranks"]
                count += 1
        shifted_count = 0
        for boundary in (0, 1):
            for shift in (1, 2):
                for matrices in (witness, cases[-1]):
                    native = profile(full, directory, 3, matrices, boundary, shift)
                    tensors = direct(3, matrices, boundary, shift)
                    rows, columns = ranks_of_coefficients(3, tensors, len(native["row_ranks"]))
                    assert rows == native["row_ranks"] and columns == native["column_ranks"]
                    shifted_count += 1
        for boundary in (0, 1):
            for matrices in ([0, 2, 4, 0], [10, 12, 10, 3], [rng.randrange(16) for _ in range(4)]):
                a = profile(full, directory, 6, matrices * 2, boundary)
                b = profile(circulant, directory, 6, matrices, boundary)
                assert a["row_ranks"] == b["row_ranks"] and a["column_ranks"] == b["column_ranks"]
        wide = profile(full, directory, 8, [9, 0, 9, 0] * 2, 0)
        evidence_path = root.parent / "evidence" / "finite-state-pairing-screen.json"
        if not evidence_path.exists():
            evidence_path = root / "finite-state-pairing-screen.json"
        baseline = json.loads(evidence_path.read_text())["selected_larger_profiles"]["markov-two-pair-baseline-t8.json"]
        assert wide["row_ranks"] == baseline["row_ranks"] and wide["column_ranks"] == baseline["column_ranks"]
    for t in range(2, 7):
        for shift in (0, 1, 2):
            assert entry(t, (0, 2, 1, 3), witness, 0, shift) == 5 << (2 * shift)
            assert entry(t, (0, 1, 2, 3), witness, 0, shift) == 0
            assert entry(t, (0, 0, 0, 0), witness, 0, shift) == 1
            assert entry(t, (1, 1, 1, 1), witness, 0, shift) == 2
    result = {
        "status": "passed", "field": "F_2",
        "independent_full_tensor_cases": count,
        "independent_delayed_transition_cases": shifted_count,
        "independent_endpoint_entries": (count + shifted_count) * 16 ** 3,
        "independent_six_position_circulant_profiles": 6,
        "packed_eight_position_baseline": 1,
        "independently_checked_short_witness_profiles": count + shifted_count + 7,
        "persistent_order_and_bit_flip_witness_lengths": list(range(2, 7)),
        "persistent_witness_offdiagonal_shifts": [0, 1, 2],
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((root / "eight_feature_crossings.cpp").read_bytes()).hexdigest(),
        "scope": "Exact endpoint and finite ranks only; no asymptotic contraction claim.",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
