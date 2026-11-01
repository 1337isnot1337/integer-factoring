#!/usr/bin/env python3
"""Exact coefficient normal forms in a prescribed source-slab affine family.

Each trial chooses a priority order for complete polynomial coefficients,
then eliminates every available leading coefficient. This is finite binary
linear algebra, not an approximation to a Smith profile. Both complete
profiles are evaluated on every resulting literal tensor.
"""

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles
from search_saturated_dual_directions import GAUGES, check_endpoint, complete_profiles


def permuted_coefficients(row, positions, degree):
    result = 0
    for coordinate, value in enumerate(row):
        while value:
            bit = value & -value
            value ^= bit
            power = bit.bit_length() - 1
            result ^= 1 << positions[coordinate * (degree + 1) + power]
    return result


def echelon_with_witness(vectors):
    pivots = {}
    for index, value in enumerate(vectors):
        expression = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value, expression
                break
            other, coefficients = pivots[pivot]
            value ^= other
            expression ^= coefficients
    return pivots


def normal_form(value, pivots):
    remainder = coefficients = 0
    while value:
        pivot = value.bit_length() - 1
        if pivot in pivots:
            other, expression = pivots[pivot]
            value ^= other
            coefficients ^= expression
        else:
            bit = 1 << pivot
            remainder ^= bit
            value ^= bit
    return remainder, coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directions', type=Path)
    parser.add_argument('--trials', type=int, default=500)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--seed', type=int, default=90731)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    rng = random.Random(args.seed)
    source = json.loads(args.directions.read_text())
    degree = source['degree']
    assert source['dimensions'] == [5] * 4 and degree == 18
    descriptors = [(i, shift) for i, primitive in enumerate(source['primitives'])
                   for shift in range(degree - primitive['normalized_degree'] + 1)]
    directions = [[value << shift for value in source['primitives'][i]['normalized_slab']]
                  for i, shift in descriptors]
    base = five_letter_source()[0]
    baseline_profiles = complete_profiles(base, GAUGES)
    cache = {name: {baseline_profiles[name]: root_for_profiles(baseline_profiles[name], energy, energy)}
             for name, energy in GAUGES.items()}
    records = {name: {'energies': energy, 'baseline_profiles': baseline_profiles[name],
                      'baseline_cost': cache[name][baseline_profiles[name]],
                      'best_rho': cache[name][baseline_profiles[name]]['rho'], 'retained': []}
               for name, energy in GAUGES.items()}
    report = {'scope': 'Finite exact coefficient-priority normal forms in the supplied endpoint-zero module; numerical scale optimization only.',
              'arguments': vars(args) | {'directions': str(args.directions), 'output': str(args.output)},
              'source_record_sha256': sha256(args.directions.read_bytes()).hexdigest(),
              'reference_digest': digest(base), 'slab_direction_descriptors': descriptors,
              'independent_slab_directions': len(directions), 'completed_trials': 0,
              'profile_queries': 0, 'duplicate_tensors': 0, 'gauges': records}
    seen = {digest(base)}
    half = (0, 5, 3, 1, 4)
    terms = [(coordinate, power) for coordinate in range(125) for power in range(degree + 1)]

    def save():
        report['elapsed_seconds'] = time.monotonic() - started
        report['distinct_profile_classes'] = {name: len(values) for name, values in cache.items()}
        args.output.write_text(json.dumps(report, indent=2) + '\n')

    for trial in range(args.trials):
        if time.monotonic() - started >= args.seconds:
            break
        scale = math.exp(rng.uniform(-3.2, -0.5))
        offset = rng.uniform(-8, 8)
        jitter = [rng.random() for _ in range(125)]

        def priority(term):
            coordinate, power = term
            y, u, v = coordinate % 5, coordinate // 5 % 5, coordinate // 25
            infinity = power - half[u] - half[v]
            zero = -power - half[y]
            if trial == 0:
                score = power
            elif trial == 1:
                score = infinity
            elif trial == 2:
                score = zero
            else:
                a, b = scale * infinity, scale * zero + offset
                maximum = max(a, b)
                score = maximum + math.log(math.exp(a - maximum) + math.exp(b - maximum))
            return score, jitter[coordinate], coordinate, power

        ordering = sorted(terms, key=priority)
        positions = [0] * (125 * (degree + 1))
        for position, (coordinate, power) in enumerate(ordering):
            positions[coordinate * (degree + 1) + power] = position
        pivots = echelon_with_witness([permuted_coefficients(row, positions, degree) for row in directions])
        assert len(pivots) == len(directions)
        parameters, values = [], [list(pair) for pair in base]
        for x in range(5):
            row = [base[x + 5 * j][0] for j in range(125)]
            _, coefficients = normal_form(permuted_coefficients(row, positions, degree), pivots)
            parameters.append(coefficients)
            active = coefficients
            while active:
                bit = active & -active
                active ^= bit
                direction = directions[bit.bit_length() - 1]
                for j, polynomial in enumerate(direction):
                    values[x + 5 * j][0] ^= polynomial
        values = tuple(map(tuple, values))
        check_endpoint(values, degree)
        key = digest(values)
        report['completed_trials'] = trial + 1
        if key in seen:
            report['duplicate_tensors'] += 1
            if trial % 25 == 0:
                save()
            continue
        seen.add(key)
        profiles = complete_profiles(values, GAUGES)
        report['profile_queries'] += 1
        for name, energy in GAUGES.items():
            profile = profiles[name]
            if profile not in cache[name]:
                cache[name][profile] = root_for_profiles(profile, energy, energy)
            cost = cache[name][profile]
            record = records[name]
            if cost['rho'] < record['best_rho']:
                record['best_rho'] = cost['rho']
                print(json.dumps({'event': 'new_best', 'trial': trial, 'gauge': name, **cost}), flush=True)
            if len(record['retained']) < 12 or cost['rho'] < record['retained'][-1]['rho']:
                record['retained'].append({'trial': trial, 'priority_scale': scale, 'priority_offset': offset,
                    'slab_parameter_masks': parameters, 'values': values, 'digest': key,
                    'profiles': profile, **cost})
                record['retained'].sort(key=lambda item: item['rho'])
                del record['retained'][12:]
        if trial % 25 == 0:
            save()
            print(json.dumps({'event': 'progress', 'trial': trial, 'queries': report['profile_queries'],
                              'best': {name: item['best_rho'] for name, item in records.items()},
                              'seconds': time.monotonic() - started}), flush=True)
    report['completed_requested_trials'] = report['completed_trials'] == args.trials
    save()
    print(json.dumps({'event': 'finished', 'trials': report['completed_trials'],
                      'queries': report['profile_queries'], 'best': {name: item['best_rho'] for name, item in records.items()}}), flush=True)


if __name__ == '__main__':
    main()
