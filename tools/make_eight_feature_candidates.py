#!/usr/bin/env python3
"""Generate a reproducible structural screen for full local parity matrices."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path

from audit_parity_matrix_screen import PRODUCT

DIAGONALS = [0, 9, 15, 16, 51, 85, 170, 204]
OFF_DIAGONALS = [255, 240, 204, 170, 192, 160, 136, 128, 34, 68, 60, 90]


def scalar_ranks(truth):
    rows, columns = [0] * 4, [0] * 4
    for x, g, h in product(range(2), repeat=3):
        if truth >> (4 * x + 2 * g + h) & 1:
            y, u, v = x ^ g, x ^ h, x ^ g ^ h
            rows[2 * x + y] |= 1 << (2 * u + v)
            columns[2 * x + u] |= 1 << (2 * y + v)
    def rank(values):
        basis = {}
        for value in values:
            while value:
                pivot = value.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = value
                    break
                value ^= basis[pivot]
        return len(basis)
    return rank(rows), rank(columns)


def efficient_diagonals():
    values = []
    for truth in range(256):
        a, b = scalar_ranks(truth), scalar_ranks(truth ^ 153)
        products = [a[0] * a[1], b[0] * b[1]]
        roots = list(map(isqrt, products))
        if all(r * r == p for r, p in zip(roots, products)) and sum(roots) == 4:
            values.append(truth)
    assert len(values) == 52
    return values


def matrices(parameters):
    a, b, p, q = parameters
    return [((a >> i) & 1) | (((p >> i) & 1) << 1)
            | (((q >> i) & 1) << 2) | (((b >> i) & 1) << 3)
            for i in range(8)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--sparse", action="store_true",
                        help="Use singleton-supported off-diagonal truth tables.")
    parser.add_argument("--all-efficient-diagonals", action="store_true",
                        help="Use all 52 scalar diagonal pencils with rank sum four.")
    args = parser.parse_args()
    off_diagonals = [1 << i for i in range(8)] if args.sparse else OFF_DIAGONALS
    diagonals = efficient_diagonals() if args.all_efficient_diagonals else DIAGONALS
    candidates = []
    discarded_flip_invariant = discarded_commuting = 0
    for parameters in product(diagonals, diagonals, off_diagonals, off_diagonals):
        values = matrices(parameters)
        if values[:4] == values[4:]:
            discarded_flip_invariant += 1
            continue
        if all(PRODUCT[a][b] == PRODUCT[b][a] for a in values for b in values):
            discarded_commuting += 1
            continue
        candidates.append(parameters)
    assert len(set(candidates)) == len(candidates)
    args.output.write_text("".join(" ".join(map(str, values)) + "\n" for values in candidates))
    metadata = {
        "scope": "Structural subset, not an exhaustive eight-matrix search. Auxiliary noncommutativity alone is not claimed to imply visible tensor order dependence.",
        "diagonal_truth_tables": diagonals,
        "offdiagonal_truth_tables": off_diagonals,
        "unfiltered_choices": len(diagonals) ** 2 * len(off_diagonals) ** 2,
        "sparse_offdiagonals": args.sparse,
        "all_efficient_diagonals": args.all_efficient_diagonals,
        "discarded_x_independent_matrices": discarded_flip_invariant,
        "discarded_commuting_matrices": discarded_commuting,
        "retained_choices": len(candidates),
        "candidate_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    if args.metadata:
        args.metadata.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
