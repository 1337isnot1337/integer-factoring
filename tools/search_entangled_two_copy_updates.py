#!/usr/bin/env python3
"""Finite search of four physical vectors that correlate two crossing copies."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from entangled_two_copy_frames import EntangledTwoCopyFrames, FIELD
from f4_recentered_lattices import f4_rank
from native_matrix_smith_f4 import MatrixOracle
from search_f4_recentered_crossings import five_letter_source, root_for_profiles
from two_copy_crossing_couplings import CoupledCopies, endpoint_power


def compact(profiles):
    return [dict(sorted(Counter(side).items())) for side in profiles]


def convolution(profiles):
    first, second = map(Counter, profiles); result = Counter()
    for a, ca in first.items():
        for b, cb in second.items():
            result[a+b] += ca*cb
    return tuple(sorted(result.items()))


def prefix_dominated(candidate, reference):
    """Candidate ascending prefix sums <= reference gives no smaller Laplace sum."""
    assert sum(c for _, c in candidate) == sum(c for _, c in reference)
    i = j = 0; remaining_a = candidate[0][1]; remaining_b = reference[0][1]
    difference = 0
    while i < len(candidate) and j < len(reference):
        count = min(remaining_a, remaining_b)
        difference += count*(candidate[i][0]-reference[j][0])
        if difference > 0:
            return False
        remaining_a -= count; remaining_b -= count
        if not remaining_a:
            i += 1
            if i < len(candidate): remaining_a = candidate[i][1]
        if not remaining_b:
            j += 1
            if j < len(reference): remaining_b = reference[j][1]
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', type=int, default=12000)
    parser.add_argument('--seconds', type=float, default=5400.)
    parser.add_argument('--seed', type=int, default=311026216)
    parser.add_argument('--replay-every', type=int, default=100)
    parser.add_argument('--entry-floor-slack', type=int, default=0)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.cases < 1 or args.replay_every < 1 or args.entry_floor_slack < 0:
        parser.error('Use positive counts and nonnegative floor slack.')
    started = time.monotonic(); deadline = started+args.seconds
    rng = random.Random(args.seed)
    values, heights, ex, ey = five_letter_source()
    compiler = CoupledCopies(5, values, 18, heights, precision=128)
    frames = EntangledTwoCopyFrames(compiler)
    baseline_profiles = compiler.baseline_control()
    energies = ([a+b for a in ex for b in ex], [a+b for a in ey for b in ey])
    reference_spectrum = convolution(baseline_profiles)
    baseline = root_for_profiles(baseline_profiles, *energies)
    half_energy = [v//2 for v in energies[0]]
    minimum = [min(side) for side in baseline_profiles]
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'entangled_two_copy_frames.py', 'native_tensor_frame_f4.py',
             'native_tensor_frame_f4.cpp', 'native_rank_one_smith_f4.py', 'native_rank_one_smith_f4.cpp',
             'pure_two_copy_frames.py', 'two_copy_crossing_couplings.py', 'f4_recentered_lattices.py',
             'native_rectangular_smith_f4.cpp', 'search_f4_recentered_crossings.py')
    record = {'scope': 'Finite structured predictions on correlated physical vectors; only separately listed original-matrix replays are certified. Numerical costs are exploratory.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names},
              'physical_alphabet': 25, 'complete_degree_bound': 36,
              'baseline': baseline, 'baseline_profiles': compact(baseline_profiles),
              'paired_energies': energies, 'profile_classes': [], 'queries': [],
              'original_matrix_replays': [], 'unresolved': [], 'counters': {},
              'completed_requested_schedule': False}
    counters = Counter(); signatures = {}; best = {'kind': 'unchanged reference', **baseline}

    def save():
        record.update(elapsed_seconds=time.monotonic()-started, counters=dict(counters), best=best,
                      strict_pressure_candidate_found=best['critical_pressure'] < -1e-9)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(record, indent=2)+'\n'); temporary.replace(args.output)

    def vector():
        threshold = rng.choice((0, 0, 1, 1, 2, 2, 3, 4, 5, 6, 7, 8))
        ceiling = rng.choice((10, 10, min(10, threshold+rng.randrange(1, 5))))
        density = rng.choice((.3, .6, 1.))
        v = [rng.randrange(1, 4) if threshold <= h <= ceiling and rng.random() < density else 0
             for h in half_energy]
        if not any(v):
            v[rng.choice([i for i, h in enumerate(half_energy) if threshold <= h <= ceiling])] = 1
        inverse = (0, 1, 3, 2)[next(c for c in v if c)]
        return [FIELD[c][inverse] for c in v]

    save()
    for index in range(args.cases):
        if time.monotonic() >= deadline:
            record['time_budget_reached'] = True; break
        vectors = [vector() for _ in range(4)]
        symmetry = index % 8
        if symmetry == 0:
            vectors = [vectors[0]]*4
        elif symmetry == 1:
            vectors[3], vectors[2] = vectors[0], vectors[1]
        elif symmetry == 2:
            vectors[1], vectors[3] = vectors[0], vectors[2]
        ranks = [f4_rank([v[5*i:5*(i+1)] for i in range(5)]) for v in vectors]
        mins = [min(h for h, c in zip(half_energy, v) if c) for v in vectors]
        lower = max(0, minimum[1]-mins[0]-mins[1]-args.entry_floor_slack)
        top = min(36, mins[2]+mins[3]-minimum[0]+args.entry_floor_slack)
        span = top-lower
        if span <= 0:
            counters['empty_endpoint_polynomial_windows'] += 1; continue
        powers = sorted({span, *(p for p in (1, 2, 4, 8, 16, 32) if p <= span)})
        power = span if index % 3 else rng.choice(powers)
        shift = rng.choice((lower, top-power))
        polynomial = endpoint_power(power) << shift; coefficient = rng.choice((1, 1, 2, 3))
        recipe = {'index': index, 'monomial_shift': shift, 'endpoint_power': power,
                  'F4_coefficient': coefficient, 'physical_vector_matrix_ranks': ranks,
                  'vector_digest': sha256(bytes(c for v in vectors for c in v)).hexdigest()}
        counters['structured_profile_queries'] += 2
        try:
            profiles, _ = frames.profiles(vectors, polynomial, coefficient)
            signature = tuple(tuple(side) for side in profiles)
            if signature not in signatures:
                spectrum = convolution(profiles)
                dominated = prefix_dominated(spectrum, reference_spectrum)
                item = {'first_case': index, 'profiles': compact(profiles),
                        'predicted_prefix_dominated_by_reference': dominated}
                if not dominated:
                    item.update(root_for_profiles(profiles, *energies)); counters['numerical_optimizations'] += 1
                signatures[signature] = len(record['profile_classes']); record['profile_classes'].append(item)
            profile_class = signatures[signature]; item = record['profile_classes'][profile_class]
            dominated = item['predicted_prefix_dominated_by_reference']
            counters['predicted_dominated_queries' if dominated else 'predicted_nondominated_queries'] += 1
            improving = not dominated and item['rho'] < best['rho']-1e-11
            replayed = False
            if index % args.replay_every == 0 or improving:
                original = tuple(compiler.integer_orders(MatrixOracle(matrix, 625, precision=128).profile()['orders'],
                                                         2*compiler.scalars[side])
                                 for side, matrix in enumerate(frames.original_matrices(vectors, polynomial, coefficient)))
                counters['complete_original_matrix_profiles'] += 2
                if original != profiles:
                    record['unresolved'].append({**recipe, 'vectors': vectors,
                        'reason': 'Structured/original profile disagreement',
                        'predicted': compact(profiles), 'original': compact(original)})
                    save(); continue
                replayed = True
                record['original_matrix_replays'].append({**recipe, 'vectors': vectors,
                    'profile_class': profile_class, 'original_profiles': compact(original)})
            record['queries'].append({**recipe, 'profile_class': profile_class,
                                       'original_matrix_replayed': replayed})
            if improving:
                best = {'kind': 'replayed correlated pure correction', **recipe, 'vectors': vectors, **item}
                print(json.dumps({'improvement': best}), flush=True)
        except (ArithmeticError, ValueError, OverflowError) as error:
            record['unresolved'].append({**recipe, 'vectors': vectors, 'reason': str(error)})
        if index % 25 == 0:
            print(json.dumps({'case': index, 'classes': len(signatures), 'best': best['rho'],
                              'counters': dict(counters), 'seconds': time.monotonic()-started}), flush=True)
            save()
    else:
        record['completed_requested_schedule'] = True
    save()


if __name__ == '__main__':
    main()
