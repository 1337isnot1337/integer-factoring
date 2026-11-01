#!/usr/bin/env python3
"""Exact controls for the four-piece F4 boundary pencil."""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import isqrt
from pathlib import Path
from random import Random
import time

from check_binary_two_piece_f4_boundary import Field, QUADS, cuts, rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    precision = 32
    unit = 1 << precision
    pairs = [(r, c) for r, c in product(range(1, 5), repeat=2)
             if (r, c) not in ((1, 3), (3, 1))]
    profile_count = 0
    equality_profiles = []
    for profiles in product(pairs, repeat=4):
        if any(sum(profiles[i][cut] for i in range(4) if i != omitted) < 4
               for omitted in range(4) for cut in range(2)):
            continue
        lower = sum(weight*isqrt(r*c << (2*precision))
                    for weight, (r, c) in zip((1, 1, 2, 2), profiles))
        assert lower >= 8*unit
        if lower == 8*unit:
            equality_profiles.append(profiles)
        profile_count += 1

    field = Field(4)
    rng = Random(20261031)
    omega, square = 2, 3
    h = [1, 1, 2, 2]
    endpoint_checks = transition_checks = triple_checks = lower_rows = 0
    initial = [rank([[int(x == v) if i < 2 else 0 for v in range(2)]
                     for x in range(2)], field) for i in range(4)]
    assert initial == [2, 2, 0, 0]
    terminals = [rank([[int(y == u == 0) if i == 0 else int(y == u == 1) if i == 1
                       else int(y == u == 0) ^ ((omega if i == 2 else square)
                                                if y == u == 1 else 0)
                       for u in range(2)] for y in range(2)], field)
                 for i in range(4)]
    assert terminals == h

    # These controls force F=0, G=0, A=0, B=0 and all-zero free tensors.
    families = []
    for mode in ('zero', 'F_zero', 'G_zero', 'A_zero', 'B_zero'):
        families.append((mode, [({quad: rng.randrange(4) for quad in QUADS},
                                 {quad: rng.randrange(4) for quad in QUADS})
                                for _ in range(4)]))
    for _ in range(64):
        families.append(('free', [({quad: rng.randrange(4) for quad in QUADS},
                                    {quad: rng.randrange(4) for quad in QUADS})
                                   for _ in range(4)]))
    zero_attainment = None
    for mode, frees in families:
        tensors = [[{} for _ in range(4)] for _ in range(4)]
        for x, y, u, v in QUADS:
            beta = [int(y == u == 0), int(y == u == 1),
                    int(y == u == 0) ^ (omega if y == u == 1 else 0),
                    int(y == u == 0) ^ (square if y == u == 1 else 0)]
            pair = [int(2*x+v == j) for j in range(4)]
            zero = [0]*4
            kets = [[int(x == v)]+pair+pair,
                    [int(x == v)]+[field.mul(omega, a) for a in pair]
                    +[field.mul(square, a) for a in pair],
                    [0]+pair+zero, [0]+zero+pair]
            endpoint = [0]*9
            for i in range(4):
                for coordinate in range(9):
                    endpoint[coordinate] ^= field.mul(kets[i][coordinate], beta[i])
            assert endpoint == [int(x == v and y == u)]+[0]*8
            endpoint_checks += 1
            for i, (fs, gs) in enumerate(frees):
                d = beta[i] if x == v else 0
                f, g = fs[x, y, u, v], gs[x, y, u, v]
                if mode == 'zero':
                    f = g = 0
                elif mode == 'F_zero':
                    f = 0
                elif mode == 'G_zero':
                    g = 0
                elif mode == 'A_zero':
                    g = f ^ d
                elif mode == 'B_zero':
                    g = field.mul(omega, field.mul(omega, f) ^ d)
                a = f ^ g ^ d
                b = field.mul(omega, f) ^ field.mul(square, g) ^ d
                values = [a, b, f, g]
                assert d == a ^ f ^ g
                assert d == b ^ field.mul(omega, f) ^ field.mul(square, g)
                assert d == field.mul(omega, a) ^ field.mul(square, b) ^ field.mul(square, f)
                assert d == field.mul(square, a) ^ field.mul(omega, b) ^ field.mul(omega, g)
                triple_checks += 4
                # A direct scalar-vector dot product checks the coefficients
                # of each transition without relying on its stated formula.
                for j, coefficients in enumerate(((1, 1, 1), (1, omega, square),
                                                    (0, 1, 0), (0, 0, 1))):
                    assert field.dot(coefficients, (d, f, g)) == values[j]
                    tensors[i][j][x, y, u, v] = values[j]
                    transition_checks += 1
        products = []
        for i in range(4):
            row = []
            for j in range(4):
                r, c = cuts(tensors[i][j], field)
                row.append(r*c)
            products.append(row)
            lower = sum(weight*isqrt(p << (2*precision)) for weight, p in zip(h, row))
            assert lower >= 4*h[i]*unit
            lower_rows += 1
        if mode == 'zero':
            assert products == [[4, 4, 0, 0], [4, 4, 0, 0],
                                [16, 16, 0, 0], [16, 16, 0, 0]]
            zero_attainment = {'rank_products': products, 'live_states': [0, 1], 'exact_base': 4}

    result = {
        'positive_rank_profiles_checked': profile_count,
        'equality_rank_profiles': equality_profiles,
        'F4_free_tensor_families_checked': len(families),
        'forced_zero_cases': ['F=G=0', 'F=0', 'G=0', 'A=0', 'B=0'],
        'vector_endpoint_checks': endpoint_checks,
        'transition_entry_checks': transition_checks,
        'triple_identity_checks': triple_checks,
        'exact_lower_vector_rows': lower_rows,
        'initial_boundary_ranks': initial,
        'terminal_boundary_ranks': terminals,
        'lower_vector': h,
        'attaining_construction': zero_attainment,
        'scope': 'Exact finite controls; the arbitrary-tensor extension-field statement is proved in the companion note.',
        'seconds': time.monotonic()-started,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
