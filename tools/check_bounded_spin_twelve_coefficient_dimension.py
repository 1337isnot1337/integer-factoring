#!/usr/bin/env python3
"""Exact dimension of the twelve-site bounded-spin coefficient algebra."""
from __future__ import annotations

from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry
from native_gf5_flint import rank, rref

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-coefficient-dimension.json'
PRIME = 5


def original_form(four_paths):
    n = len(four_paths[0])
    result = [0] * (n + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return result


def sampled_span(length, sample_size, seed):
    selected = paths(length)
    rng = random.Random(seed)
    rows, quadruples = [], []
    for _ in range(sample_size):
        quadruple = tuple(rng.choice(selected) for _ in range(4))
        rows.append(original_form(quadruple))
        quadruples.append(quadruple)
    matrix = np.array(rows, dtype=np.uint8)
    reduced, r = rref(matrix.T)
    row_pivots = [int(np.flatnonzero(row)[0]) for row in reduced[:r]]
    reduced, r = rref(matrix)
    column_pivots = [int(np.flatnonzero(row)[0]) for row in reduced[:r]]
    return matrix, quadruples, row_pivots, column_pivots


def main() -> None:
    ten_matrix, _, _, ten_pivots = sampled_span(10, 4000, 2)
    assert rank(ten_matrix) == 8
    assert ten_pivots == [0, 2, 3, 4, 5, 6, 7, 10]

    matrix, quadruples, row_pivots, pivots = sampled_span(12, 50000, 3)
    assert rank(matrix) == 10
    assert pivots == [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert len(row_pivots) == 10
    witness = [quadruples[i] for i in row_pivots]
    witness_matrix = matrix[row_pivots]
    assert rank(witness_matrix) == 10
    assert rank(witness_matrix[:, pivots]) == 10

    receipt = {
        'passed': True,
        'prime': PRIME,
        'length': 12,
        'sample_size': 50000,
        'rank_lower_bound': 10,
        'rank_upper_bound': 10,
        'dimension': 10,
        'upper_bound_argument': ('bounded-spin-three-jets.md: the three uniform '
                                 'jet identities give coefficient span at most n-2, '
                                 'which is 10 at length twelve.'),
        'coefficient_indices': pivots,
        'witness_quadruples': [[list(path) for path in quadruple] for quadruple in witness],
        'witness_minor_sha256': sha256(witness_matrix[:, pivots].tobytes()).hexdigest(),
        'ten_site_control': {
            'sample_size': 4000,
            'rank': 8,
            'coefficient_indices': ten_pivots,
            'published_basis_indices': [0, 2, 3, 4, 5, 6, 7, 10],
            'matches': ten_pivots == [0, 2, 3, 4, 5, 6, 7, 10],
        },
        'constructor_sha256': sha256(
            (ROOT / 'tools/bounded_spin_six_vertex.py').read_bytes()).hexdigest(),
        'scope': ('Exact dimension of the twelve-site coefficient span: the explicit '
                  'rank-10 minor gives lower bound 10 and the three jet identities '
                  'give upper bound 10. No twelve-site primary decomposition, cut-rank '
                  'bound, or score change is claimed.'),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'dimension': 10,
                      'coefficient_indices': pivots,
                      'ten_site_control': ten_pivots}))


if __name__ == '__main__':
    main()
