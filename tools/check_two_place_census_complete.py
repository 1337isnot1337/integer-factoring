#!/usr/bin/env python3
"""Verify full registered coverage, canonical frequencies, and every class representative."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import lzma
from pathlib import Path

from run_two_place_f4_census import TOTAL, q_planes, read_input, replay_class
from two_place_f4_census import build_input


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    archive = root / 'evidence/two-place-f4-census-complete.json.xz'
    packet = json.loads(lzma.decompress(archive.read_bytes()))
    manifest, summary = packet['manifest'], packet['summary']
    assert summary['status'] == 'complete' and summary['visited'] == summary['registered'] == TOTAL
    assert summary['unresolved'] == summary['endpoint_failures'] == 0
    for name, expected in manifest['source_sha256'].items():
        assert sha256(packet['as_run_sources'][name].encode()).hexdigest() == expected
    text, construction = build_input()
    assert sha256(text).hexdigest() == manifest['input_sha256']
    entries, _ = read_input(root / 'evidence/two-place-f4-census-input.txt')
    multiplier, offset = manifest['multiplier'], manifest['offset']
    assert multiplier % 2 == 1 and 0 <= offset < TOTAL
    inverse = pow(multiplier, -1, TOTAL)
    frequencies, position = Counter(), 0
    for chunk_record in packet['chunks']:
        chunk = chunk_record['data']
        canonical_bytes = (json.dumps(chunk, indent=2, sort_keys=True) + '\n').encode()
        assert sha256(canonical_bytes).hexdigest() == chunk_record['sha256']
        assert chunk['start'] == position and chunk['stop'] > position
        assert chunk['multiplier'] == multiplier and chunk['offset'] == offset
        assert chunk['visited'] == chunk['stop']-chunk['start']
        assert chunk['unresolved'] == chunk['endpoint_failures'] == 0
        assert sum(record['count'] for record in chunk['classes']) == chunk['visited']
        for record in chunk['classes']:
            code = record['representative']
            assert chunk['start'] <= inverse * (code-offset) % TOTAL < chunk['stop']
            assert record['q_planes'] == q_planes(code)
            assert len(record['orders']) == 50 and sum(record['orders'][25:]) == 444
            frequencies[tuple(record['orders'])] += record['count']
        position = chunk['stop']
    assert position == TOTAL and len(packet['chunks']) == 128
    assert sum(frequencies.values()) == TOTAL and len(frequencies) == 67
    canonical = sorted([{'orders': list(orders), 'count': count}
                        for orders, count in frequencies.items()], key=lambda row: row['orders'])
    assert canonical == packet['canonical_profile_frequencies']
    canonical_hash = sha256(json.dumps(canonical, separators=(',', ':'), sort_keys=True).encode()).hexdigest()
    assert canonical_hash == packet['canonical_profile_frequencies_sha256']
    assert frequencies == Counter({tuple(row['orders']): row['count'] for row in summary['classes']})
    independent, determinant_histogram = [], Counter()
    for record in summary['classes']:
        assert record['q_planes'] == q_planes(record['representative'])
        independent.append({'representative': record['representative'],
                            'q_planes': record['q_planes'], 'replay': replay_class(record, entries)})
        determinant_histogram[sum(record['orders'][:25])] += record['count']
    # Independent counts from the all-field exact determinant classification.
    expected_histogram = Counter({44: 12, 48: 3, 52: 1})
    expected_histogram.update({52-4*d: 3*4**(d-1) for d in range(3, 12)})
    assert determinant_histogram == expected_histogram
    report = {'status': 'pass', 'registered_coefficients_covered': TOTAL,
              'contiguous_chunks': len(packet['chunks']), 'complete_profile_classes': len(frequencies),
              'independent_fraction_free_complete_orders': 50*len(independent),
              'independent_class_replays': independent,
              'canonical_profile_frequencies_sha256': canonical_hash,
              'opposite_order_sum_histogram': dict(sorted(determinant_histogram.items())),
              'archive_sha256': sha256(archive.read_bytes()).hexdigest(),
              'construction': construction,
              'scope': 'Checks the entire registered coverage/count/hash ledger and independently '
              'replays every distinct complete profile representative. It does not repeat all '
              '4,194,304 native profile computations or certify numerical supporting-line minima. '
              'No all-field asymmetric obstruction or improved factoring exponent is claimed.',
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in
                                ['tools/check_two_place_census_complete.py', 'tools/two_place_f4_census.py',
                                 'tools/weighted_smith_oracle_f4.py', 'tools/check_joint_energy_conjugations.py']}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'registered_coefficients_covered', 'complete_profile_classes',
                       'independent_fraction_free_complete_orders')}))


if __name__ == '__main__':
    main()
