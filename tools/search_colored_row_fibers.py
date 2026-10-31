#!/usr/bin/env python3
"""Finite exact F4 row-frame walks with complete two-cut profile evaluation."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from colored_row_fibers import ColoredRowFiber
from typed_charge_lattice_synthesis import TypedChargeLattice, natural_heights
from search_full_alphabet_two_cut_lattices import root_for_profiles
from search_full_alphabet_reference_walk import paired_energies


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frames', type=int, default=40)
    parser.add_argument('--trials-per-frame', type=int, default=250)
    parser.add_argument('--operations-per-cut', type=int, default=2)
    parser.add_argument('--lower-window', type=int, default=8)
    parser.add_argument('--upper-window', type=int, default=12)
    parser.add_argument('--shift-window', type=int, default=3)
    parser.add_argument('--samples', type=int, default=2)
    parser.add_argument('--charge', choices=('zero', 'count', 'single-species'), default='zero')
    parser.add_argument('--negative-boundary', action='store_true')
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--seed', type=int, default=318070)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic(); rng = random.Random(args.seed)
    heights, _ = natural_heights(4)
    values = None
    reference = None
    if args.reference:
        reference = json.loads(args.reference.read_text())
        heights = reference['heights_infinity_then_zero']
        values = tuple((pair[0], 0) for pair in reference['values'])
    charges = {'zero': [0]*4, 'count': [1]*4, 'single-species': [2, 2, 0, 0]}[args.charge]
    lattice = TypedChargeLattice(4, heights, raise_cap=args.upper_window,
        charge_weights=charges, negative_boundary=args.negative_boundary, reference_values=values)
    reference_values = lattice.reference_values
    floors = tuple(max(0, order-args.lower_window) for order in lattice.orders)
    prepared = ColoredRowFiber(lattice, floors, ()).prepared
    energies = paired_energies(heights, 4)
    reference_profiles = lattice.profiles(reference_values)
    reference_cost = root_for_profiles(reference_profiles, *energies)
    source_names = [Path(__file__).name, 'colored_row_fibers.py', 'typed_charge_lattice_synthesis.py',
        'f4_lattice_linear.py', 'full_alphabet_two_cut_lattices.py', 'finite_precision_smith_frames.py',
        'search_full_alphabet_two_cut_lattices.py', 'search_full_alphabet_reference_walk.py',
        'weyl_lattice_filter_synthesis.py', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp']
    record = {'status': 'Exact finite colored-frame search; numerical full-cost roots; no all-family claim.',
        'arguments': vars(args) | {'reference': str(args.reference) if args.reference else None, 'output': str(args.output)},
        'reference': reference, 'reference_values': reference_values,
        'heights_infinity_then_zero': heights, 'charge_weights': charges,
        'degree': lattice.degree, 'endpoint_pole': lattice.pole, 'variables': lattice.variables,
        'fiber_variables': len(prepared[1]), 'reference_profiles': reference_profiles,
        'reference_cost': reference_cost, 'physical_paired_energies': energies,
        'row_target_floors': floors, 'row_target_ceilings': lattice.caps,
        'reference_row_targets': lattice.orders,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in source_names}}
    frames = []; classes = {}; retained = []; queries = trials = full_replays = 0
    best = reference_cost['rho']
    print(json.dumps({'event': 'constructed', 'variables': lattice.variables,
        'fiber_variables': len(prepared[1]), 'seconds': time.monotonic()-start}), flush=True)

    def save():
        record.update(frames=frames, completed_frames=sum(row['complete'] for row in frames),
            completed_requested_trials=(len(frames) == args.frames and all(row['complete'] for row in frames)),
            completed_trials=trials, complete_profile_queries=queries,
            independent_full_tensor_replays=full_replays, distinct_profile_classes=len(classes),
            classes=[{'profiles': key, 'samples': row['samples'], 'cost': row['cost']} for key, row in sorted(classes.items())],
            retained_witnesses=sorted(retained, key=lambda row: row['rho'])[:24],
            best_recorded_rho=best, elapsed_seconds=time.monotonic()-start,
            strict_crossing_candidate_found=any(row['critical_pressure'] < -1e-8 for row in retained))
        args.output.write_text(json.dumps(record, indent=2)+'\n')

    for frame_number in range(args.frames):
        if time.monotonic()-start >= args.seconds:
            break
        operations = []
        for side in range(2):
            for _ in range(args.operations_per_cut):
                i, j = rng.sample(range(lattice.size), 2)
                power = max(0, lattice.orders[side*lattice.size+i]-lattice.orders[side*lattice.size+j])
                power += rng.randrange(args.shift_window)
                operations.append((side, i, j, power, rng.choice((2, 3))))
        fiber = ColoredRowFiber(lattice, floors, operations, prepared=prepared)
        seen = set(); parents = [lattice.orders]
        frame = {'operations': operations, 'completed_trials': 0, 'complete': False}
        frames.append(frame)
        for trial in range(args.trials_per_frame):
            if time.monotonic()-start >= args.seconds:
                break
            frame['completed_trials'] = trial+1; trials += 1
            targets = list(lattice.orders if trial == 0 else rng.choice(parents))
            if trial:
                for index in rng.sample(range(len(targets)), rng.randrange(2, len(targets)//2+1)):
                    targets[index] = max(floors[index], targets[index]-rng.randrange(1, args.lower_window+1))
            pivots = fiber.solve(targets)
            assert pivots is not None
            order = list(range(len(targets))); rng.shuffle(order)
            if trial % 4 == 1:
                order.sort(key=lambda i: targets[i]+rng.random()*args.lower_window)
            for row in order:
                while targets[row] < lattice.caps[row]:
                    proposal = fiber.increment(pivots, row, targets[row])
                    if proposal is None:
                        break
                    pivots = proposal; targets[row] += 1
            target = tuple(targets)
            if target in seen:
                continue
            seen.add(target)
            if len(parents) < 128:
                parents.append(target)
            elif rng.randrange(4) == 0:
                parents[rng.randrange(1, len(parents))] = target
            for sample in range(args.samples):
                witness = fiber.sample(pivots, rng)
                profiles = lattice.profiles(witness); queries += 1
                if profiles in classes:
                    classes[profiles]['samples'] += 1
                    continue
                cost = root_for_profiles(profiles, *energies)
                classes[profiles] = {'samples': 1, 'cost': cost}
                row = {'frame': frame_number, 'trial': trial, 'sample': sample, 'targets': target,
                    'profiles': profiles, 'values': witness, 'fiber_solution_dimension': fiber.variables-len(pivots), **cost}
                if len(retained) < 24 or cost['rho'] < max(item['rho'] for item in retained):
                    retained.append(row); retained = sorted(retained, key=lambda item: item['rho'])[:24]
                if full_replays < 8 or cost['rho'] < best-1e-10:
                    assert fiber.check_witness(witness, targets) == profiles
                    full_replays += 1
                if cost['rho'] < best-1e-10:
                    best = cost['rho']
                    print(json.dumps({'event': 'new_best', 'frame': frame_number, 'trial': trial, **cost}), flush=True)
        frame['distinct_target_vectors'] = len(seen)
        frame['complete'] = frame['completed_trials'] == args.trials_per_frame
        save()
        print(json.dumps({'event': 'frame', 'frame': frame_number, 'trials': trials,
            'classes': len(classes), 'best_rho': best, 'seconds': time.monotonic()-start}), flush=True)
    save()
    print(json.dumps({'event': 'finished', 'frames': len(frames), 'trials': trials,
        'classes': len(classes), 'best_rho': best, 'complete': record['completed_requested_trials']}), flush=True)


if __name__ == '__main__':
    main()
