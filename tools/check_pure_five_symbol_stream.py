#!/usr/bin/env python3
"""Certify every profile in a frozen pure-correction search stream.

The original accelerator is not used. Explicitly checked finite frames
and the tracked-error structured engine certify each complete original
matrix profile. The regenerated parameter/profile stream must match the
historical rolling digest exactly. All-scale prefix dominance is then
checked using those certified integer orders.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import chain, islice
import json
from pathlib import Path
import random
import time

from check_f4_recentered_walks import certify_profile, literal_matrix
from convolution_prefix_dominance import convolution, dominated
from f4_recentered_lattices import f4_rank, multiply_matrix, pack
from native_certified_rank_one_smith_f4 import certified_rank_one_profile
from pure_five_symbol_frames import PureFiveSymbolFrames, literal_values
from replay_source_provenance import resolve_source_hashes
from search_f4_recentered_crossings import digest
from search_pure_five_symbol_corrections import random_case, structured_schedule


def heights(energy):
    half, zero = tuple(value//2 for value in energy), (0,)*5
    return zero, zero, half, half, half, half, zero, zero


def certify_frames(frame, energy):
    count = ranks = 0
    for side in (0, 1):
        matrix, scalar = literal_matrix(frame.base, 5, 18, heights(energy), side)
        first, second = frame.frames[side]
        left = [[pack(value) for value in row] for row in first]
        right = [[pack(value) for value in row] for row in second]
        actual = multiply_matrix(multiply_matrix(left, matrix, frame.jet_precision), right, frame.jet_precision)
        for i in range(25):
            for j in range(25):
                assert actual[i][j] == (pack(frame.diagonals[side][i]) if i == j else 0)
                count += 1
        for regular in (left, right):
            assert f4_rank([[value & 3 for value in row] for row in regular]) == 25
            ranks += 1
        certify_profile(matrix, scalar, frame.reference_profiles[side])
    return count, ranks


def certified_profiles(frame, vectors, polynomial, coefficient):
    profiles, margin, tail = [], None, None
    attempts = Counter()
    for side in (0, 1):
        d, u, v = frame.structured_data(vectors, polynomial, coefficient, side)
        result = None
        for precision in (512, 1024, 2048):
            attempts[precision] += 1
            try:
                result = certified_rank_one_profile(d, u, v, cutoff=frame.jet_precision,
                                                    work_precision=precision)
                break
            except ArithmeticError:
                if precision == 2048:
                    raise
        orders = result['orders']
        assert len(orders) == 25 and min(orders) >= 0 and max(orders) < frame.jet_precision
        gap = result['minimum_certified_pivot_error_margin']
        assert gap > 0
        margin = gap if margin is None else min(margin, gap)
        gap = frame.jet_precision-max(orders)
        tail = gap if tail is None else min(tail, gap)
        profiles.append(tuple(order-(18 if side == 0 else 0) for order in orders))
    return tuple(profiles), margin, tail, attempts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    record = json.loads(args.record.read_text())
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    base = tuple(map(tuple, record['reference_values']))
    assert digest(base) == record['reference_digest']
    gauges = record['gauges']
    frames = {name: PureFiveSymbolFrames(base, gauge['energies']) for name, gauge in gauges.items()}
    frame_entries = regular_ranks = 0
    reference_convolutions = {}
    for name, frame in frames.items():
        assert frame.reference_profiles == tuple(map(tuple, gauges[name]['reference_profiles']))
        a, b = certify_frames(frame, gauges[name]['energies'])
        frame_entries += a
        regular_ranks += b
        reference_convolutions[name] = convolution(frame.reference_profiles)
    counters = record['counters']
    assert counters['two_zero_subspace_queries'] <= 167238
    assert counters['structured_queries'] == counters['two_zero_subspace_queries']+counters['random_projective_queries']
    if record['completed_requested_schedule']:
        assert counters['two_zero_subspace_queries'] == 167238
        assert counters['random_projective_queries'] == record['arguments']['random_queries']
    rng = random.Random(record['arguments']['seed'])
    cases = chain(islice(structured_schedule(), counters['two_zero_subspace_queries']),
                  (random_case(rng) for _ in range(counters['random_projective_queries'])))
    rolling = sha256()
    profile_cache = {name: {} for name in gauges}
    domination = {name: Counter(dominated=0, nondominated=0) for name in gauges}
    attempts = Counter()
    queries = 0
    minimum_error_margin = minimum_tail_margin = None
    for vectors, polynomial, coefficient, stage in cases:
        profiles = {}
        for name, frame in frames.items():
            profiles[name], margin, tail, used = certified_profiles(frame, vectors, polynomial, coefficient)
            attempts.update(used)
            minimum_error_margin = margin if minimum_error_margin is None else min(minimum_error_margin, margin)
            minimum_tail_margin = tail if minimum_tail_margin is None else min(minimum_tail_margin, tail)
            if profiles[name] not in profile_cache[name]:
                profile_cache[name][profiles[name]] = dominated(profiles[name], reference_convolutions[name])
            domination[name]['dominated' if profile_cache[name][profiles[name]] else 'nondominated'] += 1
        rolling.update(json.dumps([vectors, polynomial, coefficient, profiles], separators=(',', ':'), sort_keys=True).encode()+b'\n')
        queries += 1
        if queries % 10000 == 0:
            print(json.dumps({'event': 'certified_prefix', 'queries': queries,
                              'minimum_pivot_error_margin': minimum_error_margin,
                              'minimum_frame_tail_margin': minimum_tail_margin,
                              'seconds': time.monotonic()-started}), flush=True)
    assert queries == counters['structured_queries']
    assert rolling.hexdigest() == record['query_and_profile_stream_sha256']
    for name, gauge in gauges.items():
        assert len(profile_cache[name]) == gauge['distinct_returned_profile_classes']
        assert domination[name]['dominated'] == gauge['dominated_queries']
        assert domination[name]['nondominated'] == gauge['nondominated_queries']

    # Retained witnesses and widely separated periodic samples additionally
    # use direct complete polynomial matrices and regular-frame certificates.
    literal_cases = {}
    for gauge in gauges.values():
        for witness in gauge['retained_dominated']+gauge['retained_nondominated']:
            key = json.dumps([witness['vectors'], witness['scalar_polynomial'], witness['coefficient']])
            literal_cases[key] = witness
    periodic = record['periodic_complete_matrix_controls']
    if periodic:
        for index in sorted({round(i*(len(periodic)-1)/11) for i in range(12)}):
            witness = periodic[index]
            key = json.dumps([witness['vectors'], witness['scalar_polynomial'], witness['coefficient']])
            literal_cases[key] = witness
    direct_profiles = endpoints = 0
    for witness in literal_cases.values():
        values = literal_values(base, witness['vectors'], witness['scalar_polynomial'], witness['coefficient'])
        assert digest(values) == witness['digest']
        if 'values' in witness:
            assert values == tuple(map(tuple, witness['values']))
        for name, frame in frames.items():
            profiles, _, _, _ = certified_profiles(frame, witness['vectors'], witness['scalar_polynomial'], witness['coefficient'])
            for side in (0, 1):
                matrix, scalar = literal_matrix(values, 5, 18, heights(gauges[name]['energies']), side)
                certify_profile(matrix, scalar, profiles[side])
                direct_profiles += 1
        endpoints += 625
    names = (Path(__file__).name, 'native_certified_rank_one_smith_f4.py',
             'native_certified_rank_one_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Complete historical pure-correction stream replay through explicit finite frames and the tracked-error profile engine; exact all-scale convolution comparison.',
              'record': str(args.record), 'record_sha256': sha256(args.record.read_bytes()).hexdigest(),
              'status': 'passed', 'historical_source_provenance': provenance,
              'checker_source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'completed_requested_schedule': record['completed_requested_schedule'],
              'certified_tensor_queries': queries, 'certified_complete_matrix_profiles': 4*queries,
              'exact_frame_identity_entries': frame_entries, 'regular_constant_frame_ranks': regular_ranks,
              'minimum_certified_pivot_error_margin': minimum_error_margin,
              'minimum_frame_tail_transfer_margin': minimum_tail_margin,
              'certified_engine_calls_by_work_precision': dict(attempts),
              'matching_parameter_profile_stream_sha256': rolling.hexdigest(),
              'gauges': {name: dict(counts) | {'distinct_profile_classes': len(profile_cache[name])}
                         for name, counts in domination.items()},
              'all_tested_profile_products_pointwise_no_better_than_reference':
                  all(not counts['nondominated'] for counts in domination.values()),
              'additional_literal_tensor_replays': len(literal_cases),
              'additional_complete_polynomial_frame_profiles': direct_profiles,
              'literal_endpoint_entries': endpoints, 'elapsed_seconds': time.monotonic()-started}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
