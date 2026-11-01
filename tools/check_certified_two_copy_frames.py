#!/usr/bin/env python3
"""Certify retained 625-by-625 controls through tracked-error structured frames."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import entangled_two_copy_frames as entangled
from native_certified_rank_one_smith_f4 import certified_rank_one_profile
from pure_two_copy_frames import PureTwoCopyFrames
from search_f4_recentered_crossings import five_letter_source
from two_copy_crossing_couplings import CoupledCopies


def verify():
    root = Path(__file__).resolve().parent
    source, heights, _, _ = five_letter_source()
    compiler = CoupledCopies(5, source, 18, heights, precision=128)
    pure = PureTwoCopyFrames(compiler)
    correlated = entangled.EntangledTwoCopyFrames(compiler)
    reports = []; margins = []; counts = Counter()
    for path in (root.parent/'evidence/pure-two-copy-frame-controls.json',
                 root.parent/'evidence/entangled-two-copy-frame-controls.json'):
        record = json.loads(path.read_text())
        for name, digest in record['source_sha256'].items():
            assert sha256((root/name).read_bytes()).hexdigest() == digest, name
        for index, case in enumerate(record['cases']):
            if case['source_alphabet'] != 5:
                continue
            claimed = [{int(k): v for k, v in row.items()} for row in case['profiles']]
            for precision in (256, 512):
                if 'first_support' in case:
                    profiles = []
                    for side in (0, 1):
                        data = pure.structured_data(case['first_support'], case['second_support'],
                                                    case['polynomial'], case['F4_coefficient'], side)
                        result = certified_rank_one_profile(*data, cutoff=128, work_precision=precision)
                        margins.append(result['minimum_certified_pivot_error_margin'])
                        counts.update(result['pivot_counts'])
                        profiles.append([a-2*compiler.scalars[side] for a in result['orders']])
                else:
                    # A process-local adapter selects the independently implemented
                    # tracked-error routine; no shared source or live search changes.
                    def checked(diagonal, left, right, cutoff=128, work_precision=512):
                        result = certified_rank_one_profile(diagonal, left, right,
                                                            cutoff=cutoff, work_precision=precision)
                        margins.append(result['minimum_certified_pivot_error_margin'])
                        counts.update(result['pivot_counts'])
                        return result
                    entangled.rank_one_profile = checked
                    profiles, _ = correlated.profiles(case['vectors'], case['polynomial'], case['F4_coefficient'])
                assert [dict(Counter(row)) for row in profiles] == claimed
                reports.append({'source_record': path.name, 'case_index': index,
                                'work_precision': precision, 'complete_profile_count': 2})
    names = (Path(__file__).name, 'native_certified_rank_one_smith_f4.cpp',
             'native_certified_rank_one_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    return {'scope': 'Certified tracked-error elimination through exact finite-jet frames, compared with separately retained original complete-matrix controls.',
            'complete_profiles_checked': 2*len(reports), 'comparisons': reports,
            'minimum_certified_pivot_error_margin': min(margins),
            'pivot_counts': dict(counts),
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names}}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
