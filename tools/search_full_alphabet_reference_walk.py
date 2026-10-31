#!/usr/bin/env python3
"""Exact coefficient target walks from an explicit binary reference and gauges."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice, CompressedFullAlphabetFiber
from finite_precision_smith_frames import build_lattice
from search_full_alphabet_two_cut_lattices import increment, root_for_profiles


def paired_energies(heights, q):
    return [
        [sum(heights[side+leg][i] for side in (0, 4) for leg in legs) for i in range(q)]
        for legs in ((0, 3), (1, 2))
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--lower-window', type=int, default=8)
    parser.add_argument('--upper-window', type=int, default=12)
    parser.add_argument('--trials', type=int, default=50000)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--seed', type=int, default=311843)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    reference = json.loads(args.reference.read_text())
    q, degree = reference['q'], reference['degree']
    heights = reference['heights_infinity_then_zero']
    values = tuple((value[0], 0) for value in reference['values'])
    energies = paired_energies(heights, q)
    lattice = build_lattice(values, q, degree, heights, raise_cap=args.upper_window)
    row_count = len(lattice.orders)
    floors = tuple(max(0, value-args.lower_window) for value in lattice.orders)
    fiber = CompressedFullAlphabetFiber(lattice, floors)
    profiles = lattice.profiles(values)
    reference_cost = root_for_profiles(profiles, *energies)
    print(json.dumps({'event': 'constructed', 'seconds': time.monotonic()-started,
                      'variables': lattice.variables, 'fiber_variables': fiber.variables,
                      'reference': reference_cost}), flush=True)
    rng = random.Random(args.seed)
    parents, seen = [lattice.orders], {lattice.orders}
    classes, retained = {}, []
    best = reference_cost['rho']
    queries = full_replays = trials = 0
    source_names = [Path(__file__).name, 'full_alphabet_two_cut_lattices.py', 'finite_precision_smith_frames.py',
                    'search_full_alphabet_two_cut_lattices.py', 'weyl_lattice_filter_synthesis.py',
                    'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp']
    record = {
        'status': 'Exact coefficient feasibility; numerical root optimization; finite target walk only.',
        'arguments': vars(args) | {'reference': str(args.reference), 'output': str(args.output)},
        'reference_sha256': sha256(args.reference.read_bytes()).hexdigest(),
        'reference': reference, 'binary_reference_values': values,
        'degree': degree, 'physical_alphabet': q, 'heights_infinity_then_zero': heights,
        'physical_paired_energies': energies, 'variables': lattice.variables,
        'fiber_variables': fiber.variables, 'reference_profiles': profiles,
        'reference_cost': reference_cost, 'row_target_floors': floors,
        'row_target_ceilings': lattice.caps, 'reference_row_targets': lattice.orders,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in source_names},
    }

    def save(completed=False):
        record.update(
            completed_requested_trials=completed, completed_trials=trials,
            distinct_target_vectors=len(seen), complete_profile_queries=queries,
            independent_full_tensor_replays=full_replays, distinct_profile_classes=len(classes),
            classes=[{'profiles': key, 'samples': value} for key, value in sorted(classes.items())],
            retained_witnesses=sorted(retained, key=lambda row: row['rho'])[:40],
            elapsed_seconds=time.monotonic()-started, best_recorded_rho=best,
            strict_crossing_candidate_found=any(row['critical_pressure'] < -1e-8 for row in retained),
        )
        args.output.write_text(json.dumps(record, indent=2)+'\n')

    for trial in range(args.trials):
        if time.monotonic()-started >= args.seconds:
            break
        trials = trial+1
        targets = list(lattice.orders if trial == 0 else rng.choice(parents))
        if trial:
            for index in rng.sample(range(row_count), rng.randrange(2, row_count//2+1)):
                targets[index] = max(floors[index], targets[index]-rng.randrange(1, args.lower_window+1))
        pivots = fiber.solve(targets)
        assert pivots is not None
        order = list(range(row_count))
        rng.shuffle(order)
        if trial % 4 == 1:
            order.sort(key=lambda i: targets[i]+rng.random()*args.lower_window)
        for row in order:
            while targets[row] < lattice.caps[row]:
                proposal = increment(fiber, pivots, row, targets[row])
                if proposal is None:
                    break
                pivots = proposal
                targets[row] += 1
        target = tuple(targets)
        if trial and target in seen:
            continue
        seen.add(target)
        if len(parents) < 256:
            parents.append(target)
        elif rng.randrange(4) == 0:
            parents[rng.randrange(1, len(parents))] = target
        for sample in range(args.samples):
            witness = fiber.sample(pivots, rng)
            profiles = lattice.profiles(witness)
            queries += 1
            if profiles in classes:
                classes[profiles] += 1
                continue
            classes[profiles] = 1
            cost = root_for_profiles(profiles, *energies)
            row = {'trial': trial, 'sample': sample, 'targets': target, 'profiles': profiles,
                   'values': witness, 'fiber_solution_dimension': fiber.variables-len(pivots), **cost}
            if len(retained) < 40 or cost['rho'] < max(item['rho'] for item in retained):
                retained.append(row)
                retained = sorted(retained, key=lambda item: item['rho'])[:40]
            if cost['rho'] < best-1e-10 or full_replays < 4:
                assert lattice.check_witness(witness, target) == profiles
                full_replays += 1
                if cost['rho'] < best-1e-10:
                    best = cost['rho']
                    print(json.dumps({'event': 'new_best', 'trial': trial, **cost}), flush=True)
        if trial % 100 == 0:
            save()
            print(json.dumps({'event': 'progress', 'trial': trial, 'targets': len(seen),
                              'classes': len(classes), 'rho': best,
                              'seconds': time.monotonic()-started}), flush=True)
    save(trials == args.trials)
    print(json.dumps({'event': 'finished', 'trials': trials, 'classes': len(classes),
                      'best_rho': best, 'completed_requested_trials': trials == args.trials}), flush=True)


if __name__ == '__main__':
    main()
