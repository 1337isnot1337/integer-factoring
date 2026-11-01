#!/usr/bin/env python3
"""Independently verify the four-feature parity matrix-product family over F_2."""
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


def entry_polynomial(t, labels, matrices, boundary):
    x, y, u, v = labels
    if x ^ y ^ u ^ v:
        return 0
    value = [1, 0, 0, 1]
    for i in range(t):
        g = ((x ^ y) >> i) & 1
        h = ((x ^ u) >> i) & 1
        b = matrices[2 * g + h]
        factor = [int(g == h and r == c) ^ (3 if b >> (2 * r + c) & 1 else 0)
                  for r in range(2) for c in range(2)]
        value = matrix_product(value, factor)
    return value[0] if boundary == 0 else value[1] ^ value[2] ^ value[3]


def direct_tensors(t, matrices, boundary):
    dimension = 1 << t
    tensors = {}
    for x, y, u, v in product(range(dimension), repeat=4):
        polynomial = entry_polynomial(t, (x, y, u, v), matrices, boundary)
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        index = ((x * dimension + y) * dimension + u) * dimension + v
        for degree in range(t + 1):
            if polynomial >> degree & 1:
                tensors[degree] = tensors.get(degree, 0) | (1 << index)
    return tensors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    rng = random.Random(20261030)
    matrices = [[0, 0, 0, 0], [9, 0, 9, 0], [1, 2, 4, 8]]
    matrices += [[rng.randrange(16) for _ in range(4)] for _ in range(9)]
    with TemporaryDirectory(prefix="parity-crossing-") as temporary:
        directory = Path(temporary)
        binary = directory / "parity_crossings"
        subprocess.run(["g++", "-O2", "-std=c++17", str(root / "parity_matrix_crossings.cpp"),
                        "-o", str(binary)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        count = 0
        for boundary in (0, 1):
            for choice in matrices:
                output = directory / "profile.json"
                subprocess.run([str(binary), "profile", "3", str(output), str(boundary),
                                *map(str, choice)], check=True, stdout=subprocess.DEVNULL)
                native = json.loads(output.read_text())
                tensors = direct_tensors(3, choice, boundary)
                rows, columns = ranks_of_coefficients(3, tensors, 4)
                assert rows == native["row_ranks"]
                assert columns == native["column_ranks"]
                count += 1
    # The E01/E10 order is detected by the (0,0) boundary entry.
    choice = [1, 2, 4, 8]
    assert entry_polynomial(2, (0, 2, 1, 3), choice, 0) == 5
    assert entry_polynomial(2, (0, 1, 2, 3), choice, 0) == 0
    # Vanishing mixed difference is precisely the old affine span in g,h.
    assert choice[0] ^ choice[1] ^ choice[2] ^ choice[3] != 0
    # This variant has the identity letter K_(0,0)=I, so the order witness
    # persists at every longer word length by padding with that letter.
    padded_choice = [0, 2, 4, 0]
    for t in range(2, 7):
        assert entry_polynomial(t, (0, 2, 1, 3), padded_choice, 0) == 5
        assert entry_polynomial(t, (0, 1, 2, 3), padded_choice, 0) == 0
    evidence = {
        "status": "passed", "field": "F_2",
        "independent_full_tensor_rank_cases": count,
        "independent_endpoint_entries": count * 16 ** 3,
        "position_order_witness_polynomial_masks": [5, 0],
        "nonzero_mixed_feature_matrix_mask": 15,
        "identity_padding_witness_matrix_masks": padded_choice,
        "identity_padding_verified_lengths": list(range(2, 7)),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((root / "parity_matrix_crossings.cpp").read_bytes()).hexdigest(),
        "scope": "Exact endpoint and finite ranks; no asymptotic rank or factoring-score claim.",
    }
    encoded = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
