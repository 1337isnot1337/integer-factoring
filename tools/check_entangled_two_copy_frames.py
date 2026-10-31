#!/usr/bin/env python3
"""Independent exact controls for physical vectors correlating two copies."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from entangled_two_copy_frames import EntangledTwoCopyFrames, FIELD
from f4_recentered_lattices import pack, unpack, multiply_matrix, f4_rank
from native_matrix_smith_f4 import MatrixOracle
from native_tensor_frame_f4 import sandwich
from search_f4_recentered_crossings import five_letter_source
from two_copy_crossing_couplings import CoupledCopies, endpoint_power


def verify():
    rng = random.Random(311026215); checks = Counter(); cases = []
    for n in (1, 2, 3, 5, 8):
        for precision in (128, 256):
            ms = [[[(rng.getrandbits(precision), rng.getrandbits(precision))
                    for _ in range(n)] for _ in range(n)] for _ in range(3)]
            p = [[[pack(a) for a in row] for row in m] for m in ms]
            expected = multiply_matrix(multiply_matrix(p[0], p[1], precision), p[2], precision)
            assert sandwich(*ms, precision=precision) == [[unpack(a) for a in row] for row in expected]
            checks['independent_full_precision_polynomial_sandwich_entries'] += n*n
    for q in (2, 3, 5):
        if q == 5:
            values, heights, _, _ = five_letter_source(); degree = 18
        else:
            values = []
            for index in range(q**4):
                x, y, u, v = ((index//q**j) % q for j in range(4))
                identity, crossing = int(x == u and y == v), int(x == v and y == u)
                values.append((identity | ((identity ^ crossing) << 1), 0))
            zero, charged = (0,)*q, tuple(range(-q, 0))
            heights = (zero, zero, charged, charged, charged, charged, zero, zero); degree = 1
        compiler = CoupledCopies(q, values, degree, heights, precision=128)
        frames = EntangledTwoCopyFrames(compiler)
        for case in range(4):
            if q == 5 and case == 0:
                vectors = [[int(i//q == i % q) for i in range(q*q)] for _ in range(4)]
            else:
                vectors = [[rng.randrange(4) for _ in range(q*q)] for _ in range(4)]
                for vector in vectors:
                    if not any(vector): vector[0] = 1
            polynomial = (endpoint_power((2, 4, 8, 16)[case]) << (15 if case == 0 else 4)) if q == 5 else (3, 5, 6, 3)[case]
            coefficient = (1, 2, 3, 2)[case]
            start = time.monotonic(); predicted, pivots = frames.profiles(vectors, polynomial, coefficient)
            structured_seconds = time.monotonic()-start
            start = time.monotonic(); originals = frames.original_matrices(vectors, polynomial, coefficient)
            actual = tuple(compiler.integer_orders(MatrixOracle(matrix, q**4, precision=128).profile()['orders'],
                                                   2*compiler.scalars[side]) for side, matrix in enumerate(originals))
            complete_seconds = time.monotonic()-start
            assert predicted == actual
            checks['complete_original_matrix_profile_comparisons'] += 2
            indices = list(range(q**8)) if q < 5 else [0, q**8-1]+[rng.randrange(q**8) for _ in range(2048)]
            for index in indices:
                first_index, second_index = divmod(index, q**4)
                first = compiler.coordinates[first_index]; second = compiler.coordinates[second_index]
                c = coefficient
                for leg in range(4):
                    c = FIELD[c][vectors[leg][q*first[leg]+second[leg]]]
                for side, matrix in enumerate(originals):
                    i1, j1 = compiler.location(first, side); i2, j2 = compiler.location(second, side)
                    i, j = q*q*i1+i2, q*q*j1+j2
                    expected = list(compiler.base[side][i][j])
                    weight = sum(compiler.heights[4*side+leg][first[leg]]+
                                 compiler.heights[4*side+leg][second[leg]] for leg in range(4))
                    for power in range(polynomial.bit_length()):
                        if polynomial >> power & 1:
                            exponent = 2*compiler.scalars[side]+weight+(-power if side == 0 else power)
                            for plane in (0, 1):
                                if c >> plane & 1: expected[plane] ^= 1 << exponent
                    assert tuple(expected) == matrix[i][j]
                    checks['independent_literal_weighted_entries'] += 1
            cases.append({'source_alphabet': q, 'vectors': vectors, 'polynomial': polynomial,
                          'F4_coefficient': coefficient,
                          'physical_vector_matrix_ranks': [f4_rank([v[q*i:q*(i+1)] for i in range(q)]) for v in vectors],
                          'profiles': [dict(sorted(Counter(side).items())) for side in actual],
                          'structured_pivot_counts': pivots,
                          'structured_seconds': structured_seconds, 'complete_matrix_seconds': complete_seconds})
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'entangled_two_copy_frames.py', 'native_tensor_frame_f4.py',
             'native_tensor_frame_f4.cpp', 'pure_two_copy_frames.py', 'native_rank_one_smith_f4.py',
             'native_rank_one_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    return {'scope': 'Independent finite polynomial frame identities, physical coefficient controls, and original-matrix profile replay.',
            'checks': dict(checks), 'cases': cases,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    text = json.dumps(verify(), indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')
