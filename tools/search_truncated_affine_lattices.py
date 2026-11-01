#!/usr/bin/env python3
"""Exact complete-entry affine walks starting from a degree-truncated crossing.

All coefficients in the chosen degree range are free over F4 subject to the
strict endpoint and two fixed regular Smith frames. This search does not
impose paired traces or derivative constraints. Its numerical pressure
values require a separate exact certificate before any score claim.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice, CompressedFullAlphabetFiber
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_full_alphabet_two_cut_lattices import increment, root_for_profiles
from search_fused_mixed_weighted import atomic_write as save
from truncate_reflection_crossings import truncate_even_crossing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, required=True)
    parser.add_argument('--energies', default='0,12,4,16')
    parser.add_argument('--lower-window', type=int, default=6)
    parser.add_argument('--upper-window', type=int, default=8)
    parser.add_argument('--trials', type=int, default=100000)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--seed', type=int, default=312414)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    began = time.monotonic()
    deadline = began + args.seconds
    energies = list(map(int, args.energies.split(',')))
    if len(energies) != 4 or any(value % 2 for value in energies):
        raise ValueError('Use four even integer paired energies.')
    source = truncate_even_crossing(tensor(4), args.degree)
    degree = source['degree']
    half, zero = [value // 2 for value in energies], [0] * 4
    heights = [zero, zero, half, half, half, half, zero, zero]
    lattice = FullAlphabetTwoCutLattice(source['values'], 4, degree, heights, raise_cap=args.upper_window)
    floors = tuple(max(0, value - args.lower_window) for value in lattice.orders)
    fiber = CompressedFullAlphabetFiber(lattice, floors)
    reference = lattice.profiles(lattice.seed(0))
    reference_cost = root_for_profiles(reference, energies, energies)
    # A literal original-tensor replay validates the complete chosen frames.
    assert lattice.check_witness(source['values'], lattice.orders) == reference
    names = (Path(__file__).name, 'full_alphabet_two_cut_lattices.py',
             'search_full_alphabet_two_cut_lattices.py', 'reflection_dual_rail_gaussian.py',
             'truncate_reflection_crossings.py', 'weyl_lattice_filter_synthesis.py',
             'search_fused_mixed_weighted.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'scope': __doc__, 'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'source_tensor_sha256': digest(source['values']), 'degree': degree, 'physical_alphabet': 4,
              'heights_infinity_then_zero': heights, 'variables': lattice.variables,
              'fiber_variables': fiber.variables, 'reference_profiles': reference,
              'reference_cost': reference_cost, 'row_target_floors': floors,
              'row_target_ceilings': lattice.caps, 'reference_row_targets': lattice.orders}
    print(json.dumps({'event': 'constructed', 'degree': degree, 'fiber_variables': fiber.variables,
                      'reference_rho': reference_cost['rho'], 'elapsed_seconds': time.monotonic() - began}), flush=True)
    rng = random.Random(args.seed)
    parents, seen = [lattice.orders], {lattice.orders}
    classes, retained = {}, []
    best = reference_cost['rho']
    queries = trials = 0
    full_replays = 1
    threshold = False

    def snapshot(completed=False):
        report.update(completed=completed, completed_requested_trials=trials == args.trials,
                      completed_trials=trials, distinct_target_vectors=len(seen),
                      complete_profile_queries=queries, independent_full_tensor_replays=full_replays,
                      distinct_profile_classes=len(classes),
                      classes=[{'profiles': key, 'samples': value} for key, value in sorted(classes.items())],
                      retained_witnesses=sorted(retained, key=lambda row: row['rho'])[:40],
                      elapsed_seconds=time.monotonic() - began, best_recorded_rho=best,
                      strict_numerical_candidate=threshold)
        save(args.output, report)

    snapshot()
    for trial in range(args.trials):
        if time.monotonic() >= deadline or threshold:
            break
        trials = trial + 1
        targets = list(lattice.orders if trial == 0 else rng.choice(parents))
        if trial:
            for index in rng.sample(range(32), rng.randrange(2, 17)):
                targets[index] = max(floors[index], targets[index] - rng.randrange(1, args.lower_window + 1))
        pivots = fiber.solve(targets)
        assert pivots is not None
        order = list(range(32))
        rng.shuffle(order)
        if trial % 4 == 1:
            order.sort(key=lambda i: targets[i] + rng.random() * args.lower_window)
        for row in order:
            while targets[row] < lattice.caps[row]:
                proposal = increment(fiber, pivots, row, targets[row])
                if proposal is None:
                    break
                pivots = proposal
                targets[row] += 1
        target = tuple(targets)
        if trial and target in seen:
            if trial % 100 == 0:
                snapshot()
            continue
        seen.add(target)
        if len(parents) < 256:
            parents.append(target)
        elif rng.randrange(4) == 0:
            parents[rng.randrange(1, len(parents))] = target
        for sample in range(args.samples):
            values = fiber.sample(pivots, rng)
            profiles = lattice.profiles(values)
            queries += 1
            if profiles in classes:
                classes[profiles] += 1
                continue
            classes[profiles] = 1
            cost = root_for_profiles(profiles, energies, energies)
            witness = {'trial': trial, 'sample': sample, 'targets': target, 'profiles': profiles,
                       'values': values, 'tensor_sha256': digest(values),
                       'fiber_solution_dimension': fiber.variables - len(pivots), **cost}
            if len(retained) < 40 or cost['rho'] < max(row['rho'] for row in retained):
                retained.append(witness)
                retained = sorted(retained, key=lambda row: row['rho'])[:40]
            if cost['rho'] < best - 1e-10 or full_replays < 5:
                assert lattice.check_witness(values, target) == profiles
                full_replays += 1
                if cost['rho'] < best - 1e-10:
                    best = cost['rho']
                    print(json.dumps({'event': 'new_best', 'trial': trial, **cost}), flush=True)
            threshold |= cost['critical_pressure'] < -1e-8
        if trial % 25 == 0 or threshold:
            snapshot()
            print(json.dumps({'event': 'progress', 'trial': trial, 'targets': len(seen),
                              'classes': len(classes), 'best_rho': best,
                              'elapsed_seconds': time.monotonic() - began}), flush=True)
    report['stopping_reason'] = ('numerical candidate needs exact independent replay' if threshold else
                                 ('requested finite walk completed' if trials == args.trials else
                                  'requested wall-clock budget completed'))
    snapshot(True)


if __name__ == '__main__':
    main()
