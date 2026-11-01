#!/usr/bin/env python3
"""Replay multiple physical corrections against complete original cut matrices."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random
import sys
import time

from low_rank_two_copy_updates import LowRankTwoCopyFrames
from native_matrix_smith_f4 import MatrixOracle
from search_f4_recentered_crossings import five_letter_source
from two_copy_crossing_couplings import CoupledCopies, endpoint_power


def controls():
    rng = random.Random(311026219)
    updates = []
    for index in range(9):
        vectors = []
        for leg in range(4):
            vector = [0]*25
            for position in rng.sample(range(25), 3+(index+leg) % 2):
                vector[position] = rng.randrange(1, 4)
            vectors.append(vector)
        if index % 3 == 0:
            vectors[3], vectors[2] = vectors[0][:], vectors[1][:]
        updates.append({'vectors': vectors,
                        'polynomial': endpoint_power((2, 4, 8)[index % 3]) << (8+index % 4),
                        'coefficient': 1+index % 3})
    cancellation = [dict(update, coefficient=c) for update in updates[:4] for c in (1, 2, 3)]
    return [('two_terms', updates[:2]), ('six_terms', updates[:6]),
            ('nine_terms', updates), ('twelve_term_cancellation', cancellation)]


def verify():
    started = time.monotonic()
    source, heights, _, _ = five_letter_source()
    compiler = CoupledCopies(5, source, 18, heights, precision=128)
    frames = LowRankTwoCopyFrames(compiler)
    baseline = compiler.baseline_control()
    reports = []; comparisons = 0; endpoint_checks = 0; margins = []
    for name, updates in controls():
        matrices = frames.original_sum_matrices(updates)
        for matrix in matrices:
            for i, row in enumerate(matrix):
                first, second = divmod(i, 25)
                crossing_column = (5*(first % 5)+first//5)*25+5*(second % 5)+second//5
                for j, (low, high) in enumerate(row):
                    assert low.bit_count() % 2 == int(j == crossing_column) and high.bit_count() % 2 == 0
                    endpoint_checks += 1
        original = tuple(compiler.integer_orders(MatrixOracle(matrix, 625, precision=128).profile()['orders'],
                                                  2*compiler.scalars[side])
                         for side, matrix in enumerate(matrices))
        records = []
        for precision in (256, 512):
            profiles, certificates = frames.profiles(updates, work_precision=precision)
            assert profiles == original, (name, precision)
            comparisons += 2
            margins.extend(item['minimum_certified_pivot_error_margin'] for item in certificates)
            records.append({'work_precision': precision, 'certificates': certificates})
        if name == 'twelve_term_cancellation':
            assert original == baseline
        reports.append({'name': name, 'updates': updates,
                        'original_complete_profiles': [dict(sorted(Counter(side).items())) for side in original],
                        'certified_comparisons': records})
        print(json.dumps({'completed_case': name, 'elapsed_seconds': time.monotonic()-started}), file=sys.stderr, flush=True)
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'low_rank_two_copy_updates.py',
             'native_certified_low_rank_smith_f4.py', 'native_certified_low_rank_smith_f4.cpp',
             'native_certified_rank_one_smith_f4.cpp', 'native_rectangular_smith_f4.cpp',
             'entangled_two_copy_frames.py', 'pure_two_copy_frames.py', 'native_tensor_frame_f4.py',
             'native_tensor_frame_f4.cpp', 'two_copy_crossing_couplings.py')
    return {'scope': 'Complete original 625-by-625 matrices versus certified sequential pure-update profiles, including nine terms and exact 12-term cancellation.',
            'complete_profile_comparisons': comparisons, 'literal_endpoint_entries_checked': endpoint_checks,
            'minimum_certified_pivot_error_margin': min(margins), 'cases': reports,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names}}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
