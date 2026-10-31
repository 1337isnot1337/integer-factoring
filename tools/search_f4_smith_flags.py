#!/usr/bin/env python3
"""Search exact five-symbol coefficient fibers with varying F4 Smith flags."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from f4_recentered_lattices import RecenteredCrossingLattice, RecenteredFiber
from f4_smith_flag_fibers import SmithFlagChannels
from search_f4_recentered_crossings import digest, five_letter_source, increment, root_for_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flags', type=int, default=300)
    parser.add_argument('--trials-per-flag', type=int, default=200)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--lower-window', type=int, default=4)
    parser.add_argument('--upper-window', type=int, default=5)
    parser.add_argument('--paired-energies', default='0,12,6,0,8')
    parser.add_argument('--seed', type=int, default=2026103196)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.flags, args.trials_per_flag, args.samples, args.lower_window, args.upper_window) < 1 or args.seconds <= 0:
        parser.error('Counts, windows, and the time budget must be positive.')
    energy = tuple(map(int, args.paired_energies.split(',')))
    if len(energy) != 5 or any(value % 2 for value in energy):
        parser.error('Use five even integer paired energies.')
    half, zero = tuple(value//2 for value in energy), (0,)*5
    heights = (zero, zero, half, half, half, half, zero, zero)
    values = five_letter_source()[0]
    started = time.monotonic()
    deadline = started+args.seconds
    rng = random.Random(args.seed)
    lattice = RecenteredCrossingLattice(values, 5, 18, heights, args.upper_window)
    floors = tuple(max(0, value-args.lower_window) for value in lattice.orders)
    fiber = RecenteredFiber(lattice, floors)
    channels = SmithFlagChannels(lattice, fiber)
    profiles = lattice.check_witness(values, lattice.orders)
    reference_cost = root_for_profiles(profiles, energy, energy)
    cost_cache = {profiles: reference_cost}
    classes = {}
    flags = []
    retained = []
    best = reference_cost['rho']
    queries = literal_replays = 0
    source_names = (Path(__file__).name, 'f4_smith_flag_fibers.py', 'f4_recentered_lattices.py',
                    'search_f4_recentered_crossings.py', 'native_matrix_smith_f4.py',
                    'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'kind': 'equal-order-F4-flag-search',
              'scope': 'Exact fixed-floor coefficient fibers with changing regular F4 row flags; finite searches and numerical roots only.',
              'arguments': vars(args) | {'output': str(args.output)},
              'alphabet': 5, 'degree': 18, 'reference_values': values, 'reference_digest': digest(values),
              'reference_profiles': profiles, 'reference_cost': reference_cost,
              'heights_infinity_then_zero': heights, 'paired_energies': (energy, energy),
              'floor_targets': floors, 'ceiling_targets': lattice.caps, 'reference_targets': lattice.orders,
              'fiber_variables': fiber.variables, 'fiber_compression': fiber.compression,
              'equal_order_blocks': channels.blocks, 'construction_seconds': time.monotonic()-started,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in source_names}}

    def save(completed=False):
        report.update(flags=flags, complete_profile_queries=queries, literal_witness_replays=literal_replays,
                      retained_witnesses=sorted(retained, key=lambda item: item['rho']),
                      profile_classes=[{'profiles': key, 'samples': value} for key, value in sorted(classes.items())],
                      best_recorded_rho=best, completed_requested_schedule=completed,
                      completed_target_walks=sum(flag['completed_trials'] for flag in flags),
                      strict_pressure_candidate_found=any(item['critical_pressure'] < -1e-9 for item in retained),
                      elapsed_seconds=time.monotonic()-started)
        args.output.write_text(json.dumps(report, indent=2)+'\n')

    print(json.dumps({'event': 'constructed', 'floor_variables': fiber.variables,
                      'seconds': time.monotonic()-started, 'rho': best}), flush=True)
    for flag_index in range(args.flags):
        if time.monotonic() >= deadline:
            break
        flag_started = time.monotonic()
        matrices = channels.identity_mixes() if flag_index == 0 else channels.random_mixes(rng)
        mixed = channels.mix(matrices)
        assert mixed.check_witness(values, lattice.orders) == profiles
        flag = {'flag': flag_index, 'matrices': matrices, 'completed_trials': 0,
                'nonbinary': any(c in (2, 3) for matrix in matrices for row in matrix for c in row),
                'construction_seconds': time.monotonic()-flag_started}
        flags.append(flag)
        parents = [lattice.orders]
        seen = set()
        flag_classes = set()
        for trial in range(args.trials_per_flag):
            if time.monotonic() >= deadline:
                break
            targets = list(lattice.orders if trial == 0 else rng.choice(parents))
            if trial:
                for row in rng.sample(range(50), rng.randrange(2, 26)):
                    targets[row] = max(floors[row], targets[row]-rng.randrange(1, args.lower_window+1))
            pivots = mixed.solve(targets)
            assert pivots is not None
            ordering = list(range(50))
            rng.shuffle(ordering)
            if trial % 4 == 1:
                ordering.sort(key=lambda row: targets[row]+rng.random()*args.lower_window)
            for row in ordering:
                while targets[row] < lattice.caps[row]:
                    proposal = increment(mixed, pivots, row, targets[row])
                    if proposal is None:
                        break
                    pivots = proposal
                    targets[row] += 1
            targets = tuple(targets)
            if targets not in seen:
                seen.add(targets)
                if len(parents) < 128:
                    parents.append(targets)
                for sample_index in range(args.samples):
                    sample = mixed.sample(pivots, rng)
                    realized = lattice.profiles(sample)
                    queries += 1
                    classes[realized] = classes.get(realized, 0)+1
                    flag_classes.add(realized)
                    if realized not in cost_cache:
                        cost_cache[realized] = root_for_profiles(realized, energy, energy)
                    cost = cost_cache[realized]
                    if realized not in {tuple(map(tuple, item['profiles'])) for item in retained}:
                        witness = {'flag': flag_index, 'trial': trial, 'values': sample, 'profiles': realized,
                                   'targets': targets, 'digest': digest(sample),
                                   'fiber_dimension': fiber.variables-len(pivots), **cost}
                        if len(retained) < 24 or cost['rho'] < max(item['rho'] for item in retained):
                            assert mixed.check_witness(sample, targets) == realized
                            literal_replays += 1
                            retained.append(witness)
                            retained[:] = sorted(retained, key=lambda item: item['rho'])[:24]
                    if cost['rho'] < best-1e-11:
                        best = cost['rho']
                        print(json.dumps({'event': 'new_best', 'flag': flag_index, 'trial': trial, **cost}), flush=True)
            flag['completed_trials'] = trial+1
            if trial % 50 == 0:
                flag.update(distinct_target_vectors=len(seen), distinct_profile_classes=len(flag_classes))
                save()
        flag.update(distinct_target_vectors=len(seen), distinct_profile_classes=len(flag_classes),
                    completed_requested_trials=flag['completed_trials'] == args.trials_per_flag,
                    elapsed_seconds=time.monotonic()-flag_started)
        save()
        print(json.dumps({'event': 'flag_complete', 'flag': flag_index, 'targets': len(seen),
                          'profile_classes': len(flag_classes), 'total_profile_classes': len(classes),
                          'best_rho': best, 'seconds': time.monotonic()-started}), flush=True)
    completed = len(flags) == args.flags and all(flag['completed_requested_trials'] for flag in flags)
    save(completed)
    print(json.dumps({'event': 'finished', 'flags': len(flags), 'queries': queries,
                      'best_rho': best, 'completed': completed, 'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
