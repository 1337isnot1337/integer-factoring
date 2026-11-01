#!/usr/bin/env python3
"""Broad pure-F4 tensor corrections, pruned by exact convolution dominance.

Structured profiles are search outputs until independently replayed.
Every retained witness, every nondominated profile class, and periodic
controls are compared with the complete original polynomial matrices.
Numerical root optimization is used only after the exact prefix test
fails to exclude an improvement over the reference at every scale.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import random
import time

from convolution_prefix_dominance import convolution, dominated
from pure_five_symbol_frames import PureFiveSymbolFrames, literal_values
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles
from search_saturated_dual_directions import GAUGES, check_endpoint, complete_profiles


FROBENIUS = (0, 1, 3, 2)


def two_zero_vectors():
    return ((1, 0, 0, 0, 0), (0, 0, 0, 1, 0), (1, 0, 0, 1, 0),
            (1, 0, 0, 2, 0), (1, 0, 0, 3, 0))


def structured_schedule():
    vector_cases = []
    for vectors in product(two_zero_vectors(), repeat=4):
        for coefficient in (1, 2, 3):
            key = tuple(value for row in vectors for value in row)+(coefficient,)
            conjugate = tuple(FROBENIUS[value] for value in key)
            if key <= conjugate:
                vector_cases.append((vectors, coefficient))
    assert len(vector_cases) == 978
    pairs = list(combinations(range(19), 2))
    pairs.sort(key=lambda pair: (not (6 <= pair[0] < pair[1] <= 10), pair))
    for a, b in pairs:
        polynomial = (1 << a) ^ (1 << b)
        for vectors, coefficient in vector_cases:
            yield vectors, (polynomial, 0), coefficient, 'two_zero_subspace'


def random_case(rng):
    vectors = []
    for leg in range(4):
        size = rng.choice((1, 2, 2, 3, 4, 5))
        support = sorted(rng.sample(range(5), size))
        row = [0]*5
        row[support[0]] = 1
        for i in support[1:]:
            row[i] = rng.randrange(1, 4)
        vectors.append(tuple(row))
    if rng.randrange(4):
        a, b = rng.sample(range(19), 2)
        polynomial = ((1 << a) ^ (1 << b), 0)
    else:
        a, b, c = sorted(rng.sample(range(19), 3))
        # Nonbinary endpoint-zero scalar: coefficients1,omega,omega^2.
        polynomial = ((1 << a) ^ (1 << c), (1 << b) ^ (1 << c))
    return tuple(vectors), polynomial, rng.randrange(1, 4), 'random_projective'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--random-queries', type=int, default=25000)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=2026103217)
    parser.add_argument('--replay-every', type=int, default=1000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.random_queries < 0 or args.seconds <= 0 or args.replay_every < 1:
        parser.error('Use a nonnegative random count and positive time/replay budgets.')
    started = time.monotonic()
    rng = random.Random(args.seed)
    base = five_letter_source()[0]
    frames = {name: PureFiveSymbolFrames(base, energy) for name, energy in GAUGES.items()}
    original_reference = complete_profiles(base, GAUGES)
    gauges, cache, reference_convolutions = {}, {}, {}
    for name, energy in GAUGES.items():
        profile = frames[name].reference_profiles
        assert profile == original_reference[name]
        cost = root_for_profiles(profile, energy, energy)
        gauges[name] = {'energies': energy, 'reference_profiles': profile, 'reference_cost': cost,
                        'best_recorded_rho': cost['rho'], 'retained_nondominated': [],
                        'retained_dominated': [], 'dominated_queries': 0, 'nondominated_queries': 0}
        cache[name] = {}
        reference_convolutions[name] = convolution(profile)
    names = (Path(__file__).name, 'pure_five_symbol_frames.py', 'convolution_prefix_dominance.py',
             'native_rank_one_smith_f4.py', 'native_rank_one_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp', 'f4_recentered_lattices.py',
             'search_saturated_dual_directions.py', 'search_f4_recentered_crossings.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'scope': 'Finite structured-engine screen of pure endpoint-zero F4 corrections; all-scale convolution-prefix pruning on returned integer profiles; retained and periodic complete-matrix replay; no untested structured profile is certified.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'reference_values': base, 'reference_digest': digest(base), 'alphabet': 5, 'degree': 18,
              'frame_precision': 128, 'structured_work_precision': 512,
              'planned_two_zero_subspace_queries': 978*171, 'planned_random_queries': args.random_queries,
              'two_zero_projective_vectors': two_zero_vectors(),
              'binary_scalar_exponents': list(combinations(range(19), 2))}
    counters = Counter(structured_queries=0, two_zero_subspace_queries=0, random_projective_queries=0,
                       complete_matrix_profile_comparisons=4, endpoint_entries_replayed=625,
                       structured_fallback_queries=0)
    pivot_counts = Counter()
    periodic_controls = []
    rolling = sha256()

    def save(completed=False):
        for name in gauges:
            gauges[name]['distinct_returned_profile_classes'] = len(cache[name])
            gauges[name]['nondominated_profile_classes'] = sum(not row['dominated'] for row in cache[name].values())
        report.update(counters=dict(counters), gauges=gauges, structured_pivot_counts=dict(pivot_counts),
                      periodic_complete_matrix_controls=periodic_controls,
                      completed_requested_schedule=completed,
                      query_and_profile_stream_sha256=rolling.hexdigest(), elapsed_seconds=time.monotonic()-started,
                      strict_pressure_candidate_found=any(row['critical_pressure'] < -1e-9
                          for gauge in gauges.values() for row in gauge['retained_nondominated']))
        args.output.write_text(json.dumps(report, indent=2)+'\n')

    def evaluate(vectors, polynomial, coefficient, stage):
        profiles, replayed, values = {}, False, None
        try:
            for name, frame in frames.items():
                profiles[name], counts = frame.profiles(vectors, polynomial, coefficient)
                for side in counts:
                    pivot_counts.update(side)
        except ArithmeticError:
            values = literal_values(base, vectors, polynomial, coefficient)
            profiles = complete_profiles(values, GAUGES)
            counters['structured_fallback_queries'] += 1

        def ensure_replay():
            nonlocal replayed, values
            if not replayed:
                values = literal_values(base, vectors, polynomial, coefficient)
                check_endpoint(values, 18)
                actual = complete_profiles(values, GAUGES)
                assert actual == profiles, (vectors, polynomial, coefficient, actual, profiles)
                counters['complete_matrix_profile_comparisons'] += 4
                counters['endpoint_entries_replayed'] += 625
                replayed = True

        parameters = {'vectors': vectors, 'scalar_polynomial': polynomial,
                      'coefficient': coefficient, 'stage': stage}
        for name, energy in GAUGES.items():
            profile, gauge = profiles[name], gauges[name]
            if profile not in cache[name]:
                excluded = dominated(profile, reference_convolutions[name])
                item = {'dominated': excluded}
                if not excluded:
                    ensure_replay()
                    item.update(root_for_profiles(profile, energy, energy))
                cache[name][profile] = item
            item = cache[name][profile]
            if item['dominated']:
                gauge['dominated_queries'] += 1
                if len(gauge['retained_dominated']) < 12 and profile not in {
                        tuple(map(tuple, row['profiles'])) for row in gauge['retained_dominated']}:
                    ensure_replay()
                    gauge['retained_dominated'].append(parameters | {'profiles': profile, 'values': values,
                        'digest': digest(values), 'convolution_prefix_dominated': True})
            else:
                gauge['nondominated_queries'] += 1
                if item['rho'] < gauge['best_recorded_rho']-1e-11:
                    gauge['best_recorded_rho'] = item['rho']
                    print(json.dumps({'event': 'new_best', 'gauge': name, **parameters, **item}), flush=True)
                existing = {tuple(map(tuple, row['profiles'])) for row in gauge['retained_nondominated']}
                if profile not in existing and (len(existing) < 24
                        or item['rho'] < max(row['rho'] for row in gauge['retained_nondominated'])):
                    ensure_replay()
                    gauge['retained_nondominated'].append(parameters | {'profiles': profile, 'values': values,
                        'digest': digest(values), 'convolution_prefix_dominated': False, **item})
                    gauge['retained_nondominated'][:] = sorted(gauge['retained_nondominated'], key=lambda row: row['rho'])[:24]
        counters['structured_queries'] += 1
        counters[stage+'_queries'] += 1
        rolling.update(json.dumps([vectors, polynomial, coefficient, profiles], separators=(',', ':'), sort_keys=True).encode()+b'\n')
        if counters['structured_queries'] % args.replay_every == 0:
            ensure_replay()
            if len(periodic_controls) < 256:
                periodic_controls.append(parameters | {'query': counters['structured_queries'],
                    'profiles': profiles, 'digest': digest(values)})
        if counters['structured_queries'] % 2000 == 0:
            save()
            print(json.dumps({'event': 'checkpoint', **dict(counters),
                'classes': {name: len(rows) for name, rows in cache.items()},
                'nondominated': {name: gauge['nondominated_queries'] for name, gauge in gauges.items()},
                'best': {name: gauge['best_recorded_rho'] for name, gauge in gauges.items()},
                'seconds': time.monotonic()-started}), flush=True)

    for case in structured_schedule():
        if time.monotonic()-started >= args.seconds:
            break
        evaluate(*case)
    if counters['two_zero_subspace_queries'] == 978*171:
        for trial in range(args.random_queries):
            if time.monotonic()-started >= args.seconds:
                break
            evaluate(*random_case(rng))
    completed = (counters['two_zero_subspace_queries'] == 978*171
                 and counters['random_projective_queries'] == args.random_queries)
    save(completed)
    print(json.dumps({'event': 'finished', **dict(counters), 'completed': completed,
        'best': {name: gauge['best_recorded_rho'] for name, gauge in gauges.items()},
        'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
