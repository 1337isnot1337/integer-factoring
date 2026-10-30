#!/usr/bin/env python3
"""Exact controls for the binary three-pairing observable rank detector."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import time

from check_binary_two_piece_f4_boundary import Field, QUADS, cuts, rank, kron


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    e = {q: int(q[0] == q[3] and q[1] == q[2]) for q in QUADS}
    r = {q: int(q[0] == q[1] and q[2] == q[3]) for q in QUADS}
    c = {q: int(q[0] == q[2] and q[1] == q[3]) for q in QUADS}
    records = {}
    for q in (2, 4):
        field = Field(q)
        cases = []
        for aa, bb, cc in product(range(q), repeat=3):
            tensor = {quad: field.add(field.add(field.mul(aa, e[quad]), field.mul(bb, r[quad])),
                                      field.mul(cc, c[quad])) for quad in QUADS}
            rr, cr = cuts(tensor, field)
            expected_r = 4 if field.add(aa, cc) else int(aa != 0)+int(bb != 0)
            expected_c = 4 if field.add(aa, bb) else int(aa != 0)+int(cc != 0)
            assert (rr, cr) == (expected_r, expected_c)
            detector = [[aa, field.add(aa, bb)], [field.add(aa, cc), 0]]
            dr = rank(detector, field)
            assert (dr == 0) == (aa == bb == cc == 0)
            assert rr*cr >= 4*dr*dr
            if q == 2:
                assert rr*cr == 4*dr*dr
            cases.append({'coefficients': [aa, bb, cc], 'cut_ranks': [rr, cr], 'detector_rank': dr})
        records[str(q)] = {'coefficient_cases': len(cases), 'all_cases': cases}

    field = Field(2)
    h = {quad: int(all(value == 0 for value in quad)) for quad in QUADS}
    assert cuts(h, field) == (1, 1)
    assert cuts({quad: e[quad] ^ h[quad] for quad in QUADS}, field) == (3, 3)
    crossing_image = [[1, 1], [1, 0]]
    assert rank(crossing_image, field) == 2
    tensorized = [[1]]
    tensorized_ranks = []
    for length in range(5):
        assert rank(tensorized, field) == 2**length
        tensorized_ranks.append(2**length)
        tensorized = kron(tensorized, crossing_image, field)

    # Strict observable example from observable-transition-realization.md.
    coefficients = [[(1, 1, 1), (1, 1, 1), (1, 1, 0)],
                    [(1, 1, 1), (1, 1, 1), (1, 1, 0)],
                    [(1, 0, 1), (1, 0, 1), (1, 0, 0)]]
    n = [[rank([[aa, aa ^ bb], [aa ^ cc, 0]], field) for aa, bb, cc in row]
         for row in coefficients]
    assert n == [[1, 1, 1], [1, 1, 1], [1, 1, 2]]
    vector = [1, 1, 1]
    path_bounds = []
    for length in range(8):
        bound = sum(vector)
        assert bound >= 2**length
        path_bounds.append(bound)
        vector = [sum(a*b for a, b in zip(row, vector)) for row in n]
    result = {
        'field_controls': records,
        'detector_crossing_rank': 2,
        'universal_cost_multiplier': 2,
        'tensorized_crossing_ranks': tensorized_ranks,
        'strict_observable_example_detector_rank_matrix': n,
        'strict_example_path_rank_upper_bounds': path_bounds,
        'coordinate_tensor_escape': {'H_cut_ranks': [1, 1], 'E_plus_H_cut_ranks': [3, 3]},
        'scope': 'Proof bounds the word-rank spectral transfer; it does not bound arbitrary grouped-piece true ranks.',
        'seconds': time.monotonic()-started,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value if key != 'field_controls' else
                      {q: record['coefficient_cases'] for q, record in value.items()}
                      for key, value in result.items()}, indent=2))


if __name__ == '__main__':
    main()
