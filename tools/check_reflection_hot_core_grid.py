#!/usr/bin/env python3
"""Verify exact all-scale order domination in finite hot-core F4 grids.

The optional exhaustive pass replays all 65,536 coefficient rows per
record. Every distinct stored profile pair also gets a literal polynomial
construction, separate from the grid's packed linear-combination path.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import itertools
import json
from pathlib import Path
import time

from reflection_hot_core_profiles import profile_data
from search_reflection_hot_core_filters import compile_sources, exact_profiles, factor_from_row


def key(profiles):
    return tuple(tuple(row) for row in profiles)


def pair_orders(profiles):
    return tuple(sorted(a + b for a in profiles[0] for b in profiles[1]))


def verify(path, exhaustive=False):
    start = time.monotonic()
    record = json.loads(path.read_text())
    assert record['complete'] and not record['unresolved']
    assert record['completed_candidates'] == record['planned_candidates'] == 4 ** 8
    M, P = record['interval'], record['ramification']
    gaussian = profile_data(M, P, {0: {0: 1}})
    baseline = pair_orders(gaussian['profiles'])
    counts = Counter()
    margins = []
    for row in record['profile_classes']:
        factor = factor_from_row(row['representative'], P)
        actual = profile_data(M, P, factor)
        assert key(actual['profiles']) == key(row['profiles'])
        candidate = pair_orders(actual['profiles'])
        assert len(candidate) == len(baseline) == 64
        assert all(a <= b for a, b in zip(candidate, baseline))
        margins.append(min(b - a for a, b in zip(candidate, baseline)))
        profile = key(row['profiles'])
        assert profile not in counts and row['rows'] > 0
        counts[profile] = row['rows']
    assert sum(counts.values()) == 4 ** 8
    replayed = 0
    if exhaustive:
        compiled = [compile_sources(M, P, place) for place in ('zero', 'infinity')]
        fresh = Counter()
        for coefficients in itertools.product(range(4), repeat=8):
            fresh[exact_profiles(compiled, coefficients)] += 1
            replayed += 1
        assert fresh == counts
    return {'source': str(path), 'source_sha256': sha256(path.read_bytes()).hexdigest(),
            'interval': M, 'ramification': P, 'coefficient_rows': sum(counts.values()),
            'distinct_profile_pairs': len(counts), 'literal_profile_replays': 2 * len(counts),
            'exact_combined_order_inequalities': 64 * len(counts),
            'minimum_coordinate_margin': min(margins),
            'exhaustive_coefficient_rows_replayed': replayed,
            'all_scales_s_greater_than_one_no_gain': True,
            'seconds': time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records', nargs='*', type=Path)
    parser.add_argument('--exhaustive', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    paths = args.records or [Path(f'evidence/reflection-hot-grid-M{M}.json') for M in (4, 8, 16, 32, 64, 128)]
    records = []
    for path in paths:
        records.append(verify(path, args.exhaustive))
        print(json.dumps(records[-1]), flush=True)
    result = {'status': 'Exact finite-grid all-scale domination; natural gauges only.',
              'proof': 'Sorted pairwise Smith-order sums are coordinatewise at most the Gaussian sums.',
              'exhaustive': args.exhaustive, 'records': records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
