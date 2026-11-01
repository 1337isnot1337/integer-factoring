#!/usr/bin/env python3
"""Exact restricted row-lattice searches in saturated source directions.

Targets are sufficient fixed-frame conditions. Failure is scoped to the
displayed frames and finite schedule. Complete profile costs use numerical
physical-type optimization and are not certified minima.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from f4_recentered_lattices import RecenteredCrossingLattice
from restricted_f4_crossing_fiber import RestrictedCrossingFiber
from saturated_dual_kernel_directions import literal_directions
from search_f4_recentered_crossings import five_letter_source, root_for_profiles, increment, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directions', type=Path)
    parser.add_argument('--direction-degree', type=int, default=18)
    parser.add_argument('--gauge', choices=('natural', 'two-zero'), default='natural')
    parser.add_argument('--trials', type=int, default=1000)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--raise-cap', type=int, default=5)
    parser.add_argument('--lower-window', type=int, default=6)
    parser.add_argument('--seed', type=int, default=90631)
    parser.add_argument('--strategy', choices=('tradeoffs', 'forced-raises'), default='tradeoffs')
    parser.add_argument('--engine', choices=('projected', 'shifted'), default='projected')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    deadline = started + args.seconds
    rng = random.Random(args.seed)
    source = json.loads(args.directions.read_text())
    source['direction_descriptors'] = [descriptor for descriptor in source['direction_descriptors']
        if source['primitives'][descriptor[0]]['normalized_degree'] + descriptor[2] <= args.direction_degree]
    direction_count = len(source['direction_descriptors'])
    values, heights, ex, ey = five_letter_source()
    if args.gauge == 'two-zero':
        ex = ey = (0, 12, 6, 0, 8)
        half, zero = tuple(value // 2 for value in ex), (0,) * 5
        heights = (zero, zero, half, half, half, half, zero, zero)
    if args.engine == 'shifted':
        from source_module_crossing_fiber import SourceModuleCrossingFiber
        fiber = SourceModuleCrossingFiber(values, source, heights, args.raise_cap)
        lattice = fiber.full_lattice
    else:
        directions = literal_directions(source)
        lattice = RecenteredCrossingLattice(values, 5, 18, heights, args.raise_cap)
        fiber = RestrictedCrossingFiber(lattice, directions)
    baseline_profiles = lattice.check_witness(values, lattice.orders)
    baseline_cost = root_for_profiles(baseline_profiles, ex, ey)
    baseline_pivots = fiber.solve(fiber.orders)
    assert baseline_pivots is not None
    one_row_raises = []
    for row in range(50):
        solution = increment(fiber, baseline_pivots, row, fiber.orders[row])
        if solution is not None:
            one_row_raises.append({'row': row, 'fiber_dimension': fiber.variables - len(solution)})
    record = {'scope': 'Exact sufficient row-frame constraints in the complete displayed polynomial module; finite search with floating profile costs.',
              'arguments': vars(args) | {'directions': str(args.directions), 'output': str(args.output)},
              'directions_record_sha256': sha256(args.directions.read_bytes()).hexdigest(),
              'reference_digest': digest(values), 'reference_values': values,
              'heights_infinity_then_zero': heights, 'paired_energies': (ex, ey),
              'direction_count': direction_count, 'independent_direction_count': fiber.variables,
              'baseline_profiles': baseline_profiles, 'baseline_cost': baseline_cost,
              'reference_targets': fiber.orders, 'ceiling_targets': fiber.caps,
              'baseline_preserving_dimension': fiber.variables - len(baseline_pivots),
              'complete_single_row_raises': one_row_raises, 'completed_trials': 0,
              'relaxed_forced_raise_feasibility': [],
              'construction_seconds': time.monotonic() - started, 'retained_witnesses': []}
    parents = [fiber.orders]
    seen = set()
    classes = {}
    costs = {}
    best = baseline_cost['rho']
    queries = 0

    def save(complete=False):
        record.update(distinct_target_vectors=len(seen),
                      profile_classes=[{'profiles': profiles, 'samples': count} for profiles, count in sorted(classes.items())],
                      complete_profile_queries=queries, best_recorded_rho=best,
                      elapsed_seconds=time.monotonic() - started,
                      completed_requested_schedule=complete)
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(record, indent=2) + '\n')
        temporary.replace(args.output)

    save()
    print(json.dumps({'event': 'constructed', 'variables': fiber.variables,
                      'baseline_preserving_dimension': record['baseline_preserving_dimension'],
                      'single_row_raises': len(one_row_raises), 'rho': best,
                      'seconds': record['construction_seconds']}), flush=True)
    for trial in range(args.trials):
        if time.monotonic() >= deadline:
            break
        targets = list(fiber.orders if trial == 0 else rng.choice(parents))
        if trial and args.strategy == 'forced-raises':
            targets = [0] * 50
            forced = [(trial - 1) % 50]
            if trial > 50:
                forced = rng.sample(range(50), rng.randrange(1, 5))
            for row in forced:
                targets[row] = fiber.orders[row] + 1
        elif trial:
            for row in rng.sample(range(50), rng.randrange(2, 40)):
                targets[row] = max(0, targets[row] - rng.randrange(1, args.lower_window + 1))
        pivots = fiber.solve(targets)
        if trial and args.strategy == 'forced-raises':
            record['relaxed_forced_raise_feasibility'].append(
                {'trial': trial, 'rows': forced, 'feasible': pivots is not None})
        record['completed_trials'] = trial + 1
        if pivots is None:
            if trial % 10 == 0:
                save()
            continue
        ordering = list(range(50))
        rng.shuffle(ordering)
        if trial % 4 == 1:
            ordering.sort(key=lambda row: targets[row] + rng.random() * args.lower_window)
        for row in ordering:
            while targets[row] < fiber.caps[row]:
                proposed = increment(fiber, pivots, row, targets[row])
                if proposed is None:
                    break
                pivots = proposed
                targets[row] += 1
        targets = tuple(targets)
        if targets in seen:
            if trial % 10 == 0:
                save()
            continue
        seen.add(targets)
        if len(parents) < 128:
            parents.append(targets)
        elif rng.randrange(4) == 0:
            parents[rng.randrange(1, len(parents))] = targets
        sample = fiber.sample(pivots, rng)
        profiles = fiber.profiles(sample)
        queries += 1
        classes[profiles] = classes.get(profiles, 0) + 1
        if profiles not in costs:
            costs[profiles] = root_for_profiles(profiles, ex, ey)
        cost = costs[profiles]
        if len(record['retained_witnesses']) < 12 or cost['rho'] < max(w['rho'] for w in record['retained_witnesses']):
            assert fiber.check_witness(sample, targets) == profiles
            witness = {'values': sample, 'digest': digest(sample), 'profiles': profiles,
                       'targets': targets, 'trial': trial, 'fiber_dimension': fiber.variables - len(pivots), **cost}
            record['retained_witnesses'].append(witness)
            record['retained_witnesses'].sort(key=lambda value: value['rho'])
            del record['retained_witnesses'][12:]
        if cost['rho'] < best:
            best = cost['rho']
            print(json.dumps({'event': 'new_best', 'trial': trial, **cost}), flush=True)
        if trial % 10 == 0:
            save()
            print(json.dumps({'event': 'progress', 'trial': trial, 'targets': len(seen),
                              'classes': len(classes), 'best': best, 'seconds': time.monotonic() - started}), flush=True)
    save(record['completed_trials'] == args.trials)
    print(json.dumps({'event': 'finished', 'trials': record['completed_trials'],
                      'targets': len(seen), 'queries': queries, 'best': best}), flush=True)


if __name__ == '__main__':
    main()
