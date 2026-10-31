#!/usr/bin/env python3
"""Search low-degree saturated dual-kernel perturbations of the q5 reference.

The one-coordinate scan and scheduled two-coordinate scan use one
representative of every coefficient Frobenius orbit over F4.  Later sparse
walks are stochastic.  Every objective uses complete native Smith profiles;
its scale optimization is numerical and needs a separate exact certificate.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import time

from native_matrix_smith_f4 import MatrixOracle
from saturated_dual_kernel_directions import literal_directions
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles


GAUGES = {'natural': (0, 10, 6, 2, 8), 'two_zero': (0, 12, 6, 0, 8)}
PAIR_ORBITS = ((1, 1), (1, 2), (2, 1), (2, 2), (2, 3))


def independent_directions(directions, degree):
    pivots = {}
    indices = []
    for index, direction in enumerate(directions):
        assert all(high == 0 for low, high in direction)
        value = sum(low << ((degree+1)*entry) for entry, (low, high) in enumerate(direction))
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                indices.append(index)
                break
            value ^= pivots[pivot]
    return tuple(indices)


def specialize(base, sparse_directions, assignment):
    values = [list(pair) for pair in base]
    for index, coefficient in assignment:
        assert coefficient in (1, 2, 3)
        for entry, polynomial in sparse_directions[index]:
            if coefficient & 1:
                values[entry][0] ^= polynomial
            if coefficient >> 1:
                values[entry][1] ^= polynomial
    return tuple(map(tuple, values))


def complete_profiles(values, energies):
    matrices = [[[(0, 0)]*25 for _ in range(25)] for _ in range(2)]
    for index, pair in enumerate(values):
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        matrices[0][5*x+u][5*y+v] = pair
        matrices[1][5*x+y][5*u+v] = pair
    oracles = [MatrixOracle(matrix, 25) for matrix in matrices]
    result = {}
    for name, energy in energies.items():
        half = tuple(value//2 for value in energy)
        infinity = oracles[0].profile(row_weights=[half[u] for x in range(5) for u in range(5)],
                                     column_weights=[half[v] for y in range(5) for v in range(5)],
                                     infinity=True)['orders']
        zero = oracles[1].profile(row_weights=[half[x]+half[y] for x in range(5) for y in range(5)])['orders']
        assert all(abs(value-round(value)) < 1e-8 for value in infinity+zero)
        result[name] = tuple(map(round, infinity)), tuple(map(round, zero))
    return result


def check_endpoint(values, degree):
    for index, (low, high) in enumerate(values):
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        assert max(low, high).bit_length() <= degree+1
        assert low.bit_count() % 2 == int(x == v and y == u)
        assert high.bit_count() % 2 == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directions', type=Path)
    parser.add_argument('--pair-queries', type=int, default=55875)
    parser.add_argument('--random-queries', type=int, default=100000)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=2026103199)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.pair_queries < 0 or args.random_queries < 0 or args.seconds <= 0:
        parser.error('Query counts must be nonnegative and the time budget positive.')
    started = time.monotonic()
    deadline = started+args.seconds
    rng = random.Random(args.seed)
    source = json.loads(args.directions.read_text())
    assert source['dimensions'] == [5]*4 and source['degree'] == 18
    directions = literal_directions(source)
    indices = independent_directions(directions, 18)
    sparse = [tuple((entry, low) for entry, (low, high) in enumerate(direction) if low)
              for direction in directions]
    base = five_letter_source()[0]
    check_endpoint(base, 18)
    initial_profiles = complete_profiles(base, GAUGES)
    cache = {name: {} for name in GAUGES}
    records = {}
    for name, energy in GAUGES.items():
        cost = root_for_profiles(initial_profiles[name], energy, energy)
        cache[name][initial_profiles[name]] = cost
        records[name] = {'energies': energy, 'reference_profiles': initial_profiles[name],
                         'reference_cost': cost, 'best_recorded_rho': cost['rho'], 'retained': []}
    population = [()]
    population_seen = {()}
    counters = {'single_queries': 0, 'pair_queries': 0, 'random_queries': 0, 'literal_endpoint_replays': 1}
    source_names = (Path(__file__).name, 'saturated_dual_kernel_directions.py', 'search_f4_recentered_crossings.py',
                    'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'scope': 'Exact endpoint-zero saturated directions; complete finite sparse scans and stochastic walks; exploratory numerical costs only.',
              'arguments': vars(args) | {'directions': str(args.directions), 'output': str(args.output)},
              'direction_record_sha256': sha256(args.directions.read_bytes()).hexdigest(),
              'reference_values': base, 'reference_digest': digest(base), 'alphabet': 5, 'degree': 18,
              'independent_direction_indices': indices, 'parameter_dimension_over_F4': len(indices),
              'single_coefficient_orbits': (1, 2), 'pair_coefficient_orbits': PAIR_ORBITS,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in source_names}}

    def save(completed=False):
        for name in records:
            records[name]['distinct_profile_classes'] = len(cache[name])
        report.update(counters=counters, gauges=records, continuation_population=population,
                      completed_requested_schedule=completed, elapsed_seconds=time.monotonic()-started,
                      strict_pressure_candidate_found=any(item['critical_pressure'] < -1e-9
                          for record in records.values() for item in record['retained']))
        args.output.write_text(json.dumps(report, indent=2)+'\n')

    def evaluate(assignment, stage):
        values = specialize(base, sparse, assignment)
        profiles = complete_profiles(values, GAUGES)
        keep_population = False
        for name, energy in GAUGES.items():
            profile = profiles[name]
            if profile not in cache[name]:
                cache[name][profile] = root_for_profiles(profile, energy, energy)
            cost = cache[name][profile]
            item = records[name]
            if cost['rho'] < item['best_recorded_rho']-1e-11:
                item['best_recorded_rho'] = cost['rho']
                print(json.dumps({'event': 'new_best', 'gauge': name, 'assignment': assignment, **cost}), flush=True)
            if cost['rho'] <= item['best_recorded_rho']+0.02:
                keep_population = True
            existing = {tuple(map(tuple, witness['profiles'])) for witness in item['retained']}
            if profile not in existing and (len(item['retained']) < 24 or cost['rho'] < max(w['rho'] for w in item['retained'])):
                check_endpoint(values, 18)
                counters['literal_endpoint_replays'] += 1
                item['retained'].append({'assignment': assignment, 'stage': stage, 'profiles': profile,
                                         'values': values, 'digest': digest(values), **cost})
                item['retained'][:] = sorted(item['retained'], key=lambda w: w['rho'])[:24]
        key = tuple(map(tuple, assignment))
        if keep_population and key not in population_seen:
            population_seen.add(key)
            if len(population) < 128:
                population.append(key)
            elif rng.randrange(4) == 0:
                population[rng.randrange(1, len(population))] = key
        counters[stage+'_queries'] += 1
        total = sum(counters[name+'_queries'] for name in ('single', 'pair', 'random'))
        if total % 100 == 0:
            save()
            print(json.dumps({'event': 'checkpoint', **counters,
                              'best': {name: item['best_recorded_rho'] for name, item in records.items()},
                              'profile_classes': {name: len(value) for name, value in cache.items()},
                              'seconds': time.monotonic()-started}), flush=True)

    for index in indices:
        for coefficient in (1, 2):
            if time.monotonic() >= deadline:
                break
            evaluate(((index, coefficient),), 'single')
        if time.monotonic() >= deadline:
            break
    pairs = list(combinations(indices, 2))
    rng.shuffle(pairs)
    planned_pairs = min(args.pair_queries, len(pairs)*len(PAIR_ORBITS))
    for first, second in pairs:
        for a, b in PAIR_ORBITS:
            if counters['pair_queries'] >= planned_pairs or time.monotonic() >= deadline:
                break
            evaluate(((first, a), (second, b)), 'pair')
        if counters['pair_queries'] >= planned_pairs or time.monotonic() >= deadline:
            break
    for trial in range(args.random_queries):
        if time.monotonic() >= deadline:
            break
        parameters = dict(rng.choice(population))
        for _ in range(rng.randrange(1, 4)):
            index = rng.choice(indices)
            parameters[index] = parameters.get(index, 0) ^ rng.randrange(1, 4)
            if not parameters[index]:
                del parameters[index]
        while len(parameters) > 6:
            del parameters[rng.choice(tuple(parameters))]
        evaluate(tuple(sorted(parameters.items())), 'random')
    completed = (counters['single_queries'] == 2*len(indices) and counters['pair_queries'] == planned_pairs
                 and counters['random_queries'] == args.random_queries)
    report['planned_pair_queries'] = planned_pairs
    save(completed)
    print(json.dumps({'event': 'finished', **counters, 'completed': completed,
                      'best': {name: item['best_recorded_rho'] for name, item in records.items()},
                      'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
