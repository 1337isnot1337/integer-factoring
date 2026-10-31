#!/usr/bin/env python3
"""Replay structured pure-coupling profiles against original complete matrices."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import time

from pure_two_copy_frames import PureTwoCopyFrames
from search_f4_recentered_crossings import five_letter_source
from two_copy_crossing_couplings import CoupledCopies, constant_slab_tensors, endpoint_power


def verify(large=False):
    counts = Counter(); results = []; rng = random.Random(311026211)
    def check(compiler, frames, first, second, polynomial, coefficient, label):
        start = time.monotonic()
        fast, pivots = frames.profiles(first, second, polynomial, coefficient)
        fast_seconds = time.monotonic()-start
        start = time.monotonic()
        complete = compiler.profiles(first, second, polynomial, coefficient)
        complete_seconds = time.monotonic()-start
        assert fast == complete, (label, fast, complete)
        counts['complete_original_matrix_profile_comparisons'] += 2
        counts['complete_regular_frame_identity_entries'] += sum(x['complete_frame_identity_entries'] for x in frames.controls)
        results.append({'label': label, 'source_alphabet': compiler.q,
                        'first_support': first, 'second_support': second,
                        'polynomial': polynomial, 'F4_coefficient': coefficient,
                        'profiles': [dict(sorted(Counter(row).items())) for row in complete],
                        'structured_pivot_counts': pivots,
                        'structured_seconds': fast_seconds, 'complete_matrix_seconds': complete_seconds})
    for q, charged in ((2, False), (2, True), (3, False)):
        values = []
        for index in range(q**4):
            x, y, u, v = ((index//q**j) % q for j in range(4))
            identity, crossing = int(x == u and y == v), int(x == v and y == u)
            values.append((identity | ((identity ^ crossing) << 1), 0))
        zero = (0,)*q
        height = tuple(range(-q, 0)) if charged else zero
        heights = (zero, zero, height, height, height, height, zero, zero)
        compiler = CoupledCopies(q, values, 1, heights, precision=128)
        frames = PureTwoCopyFrames(compiler)
        for case in range(8):
            supports = []
            for _ in range(2):
                legs = [tuple(i for i in range(q) if rng.randrange(2)) for _ in range(4)]
                legs = [leg if leg else (rng.randrange(q),) for leg in legs]
                supports.append(tuple(product(*legs)))
            check(compiler, frames, *supports, rng.choice((3, 5, 6)), rng.choice((1, 2, 3)),
                  'small-'+str(q)+'-'+str(charged)+'-'+str(case))
    if large:
        source_path = Path('evidence/source-slab-module-directions-h16.json')
        source = json.loads(source_path.read_text())['minimal_basis']['generators']
        generators = [[(entry >> 8) & 1 for entry in source[index]] for index in (5, 9)]
        generators.append([1]+[0]*124)
        supports = constant_slab_tensors(generators, 5)
        values, heights, _, _ = five_letter_source()
        compiler = CoupledCopies(5, values, 18, heights, precision=128)
        frames = PureTwoCopyFrames(compiler)
        for i, j, shift, power, coefficient in ((0, 0, 0, 36, 1), (0, 5, 4, 24, 2),
                                               (5, 9, 0, 32, 3), (10, 14, 8, 12, 2),
                                               (12, 12, 0, 2, 1), (1, 6, 6, 18, 1)):
            check(compiler, frames, supports[i]['support'], supports[j]['support'],
                  endpoint_power(power) << shift, coefficient,
                  'five-'+','.join(map(str, (i, j, shift, power, coefficient))))
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'pure_two_copy_frames.py', 'native_rank_one_smith_f4.cpp',
             'native_rank_one_smith_f4.py', 'two_copy_crossing_couplings.py',
             'f4_recentered_lattices.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    return {'scope': 'Finite independent original-matrix replay of pure-coupling structured profiles; no untested structured profile is certified.',
            'large_controls_included': large, 'checks': dict(counts), 'cases': results,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--large', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.dumps(verify(args.large), indent=2)+'\n'
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end='')
