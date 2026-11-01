#!/usr/bin/env python3
"""Recenter complete-entry Smith-fiber synthesis on explicit F4 crossings.

The input tensors, source hashes, reference frames, every retained full
coefficient witness, and completed trial counts are recorded.  The
one-dimensional pressure minimization is numerical and is not a global
optimality certificate for this family or for physical gauges.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from affine_place_profile_oracle import primitive_seed
from f4_two_cut_smith_fibers import F4FullAlphabetLattice, CompressedF4Fiber
from reflection_laurent_physical_frames import digest
from search_full_alphabet_two_cut_lattices import root_for_profiles
from search_fused_mixed_weighted import atomic_write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--lower-window', type=int, default=4)
    parser.add_argument('--upper-window', type=int, default=6)
    parser.add_argument('--trials', type=int, default=20000)
    parser.add_argument('--samples', type=int, default=2)
    parser.add_argument('--heights', default='0,6,2,8')
    parser.add_argument('--seed', type=int, default=312691)
    args = parser.parse_args()
    began, rng = time.monotonic(), random.Random(args.seed)
    deadline = began + args.seconds
    original = args.input.read_bytes()
    data = json.loads(original)
    h, zero = list(map(int, args.heights.split(','))), [0]*4
    if len(h) != 4:
        raise ValueError('Supply four comma-separated integer physical heights.')
    heights = [zero, zero, h, h, h, h, zero, zero]
    names = (Path(__file__).name, 'f4_two_cut_smith_fibers.py', 'f4_regular_smith_frames.py',
             'finite_jet_smith_fibers.py', 'affine_place_profile_oracle.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'search_full_alphabet_two_cut_lattices.py')
    report = {'scope': __doc__, 'arguments': vars(args) | {'input': str(args.input), 'output': str(args.output)},
              'input_sha256': sha256(original).hexdigest(), 'heights_infinity_then_zero': heights,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'references': [], 'records': [], 'retained_witnesses': [], 'attempted_trials': 0,
              'completed_trials': 0, 'complete_profile_queries': 0, 'exact_full_replays': 0,
              'completed': False}
    workers = []
    for index, case in enumerate(data['cases']):
        seed = case['seed']
        assert digest(seed['values']) == case['tensor_sha256']
        assert seed['q_x'] == seed['q_u'] == 4
        lattice = F4FullAlphabetLattice(seed['values'], 4, seed['degree'], heights, args.upper_window)
        floors = tuple(max(0, value-args.lower_window) for value in lattice.orders)
        fiber = CompressedF4Fiber(lattice, floors)
        profiles = lattice.check_witness(seed['values'], lattice.orders)
        cost = root_for_profiles(profiles, [2*v for v in h], [2*v for v in h])
        reference = {'source_index': index, 'label': case['label'], 'degree': seed['degree'],
                     'tensor_sha256': case['tensor_sha256'], 'frame_orders': lattice.orders,
                     'floors': floors, 'caps': lattice.caps, 'fiber_variables': fiber.variables,
                     'profiles': profiles, 'cost': cost}
        report['references'].append(reference)
        report['exact_full_replays'] += 1
        workers.append({'case': case, 'lattice': lattice, 'fiber': fiber,
                        'parents': [(cost['rho'], lattice.orders)], 'seen': set(),
                        'profile_classes': set(), 'best': cost['rho']})
        report['elapsed_seconds'] = time.monotonic()-began
        atomic_write(args.output, report)
        print(json.dumps({'event': 'reference', **reference, 'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    stop = False
    for trial in range(args.trials):
        if time.monotonic() >= deadline or stop:
            break
        report['attempted_trials'] += 1
        index = trial % len(workers)
        worker = workers[index]
        lattice, fiber = worker['lattice'], worker['fiber']
        if trial < len(workers):
            targets = list(lattice.orders)
        elif trial % 7 == 0:
            targets = list(fiber.floors)
        else:
            if rng.randrange(4):
                parent = rng.choice(worker['parents'][:min(12, len(worker['parents']))])
            else:
                parent = rng.choice(worker['parents'])
            targets = list(parent[1])
            for row in rng.sample(range(32), rng.randrange(2, 21)):
                targets[row] = max(fiber.floors[row], targets[row]-rng.randrange(1, args.lower_window+1))
        pivots = fiber.solve(targets)
        assert pivots is not None
        order = list(range(32))
        rng.shuffle(order)
        if trial % 4 == 1:
            order.sort(key=lambda row: targets[row] + rng.random()*args.lower_window)
        for row in order:
            while targets[row] < lattice.caps[row]:
                proposal = fiber.increment(pivots, row, targets[row])
                if proposal is None:
                    break
                pivots = proposal
                targets[row] += 1
        targets = tuple(targets)
        if targets in worker['seen']:
            continue
        worker['seen'].add(targets)
        best_this_target = None
        for sample in range(args.samples):
            values = fiber.sample(pivots, None if sample == 0 else rng)
            seed = worker['case']['seed'] | {'values': values}
            primitive, common = primitive_seed(seed)
            profiles = lattice.profiles(primitive['values'])
            report['complete_profile_queries'] += 1
            if profiles in worker['profile_classes']:
                continue
            worker['profile_classes'].add(profiles)
            cost = root_for_profiles(profiles, [2*v for v in h], [2*v for v in h])
            best_this_target = min(cost['rho'], best_this_target if best_this_target is not None else cost['rho'])
            row = {'trial': trial, 'source_index': index, 'label': worker['case']['label'], 'sample': sample,
                   'targets': targets, 'fiber_solution_dimension': fiber.variables-len(pivots),
                   'raw_tensor_sha256': digest(values), 'tensor_sha256': digest(primitive['values']),
                   'actual_degree': primitive['degree'], 'common_factor': common, 'profiles': profiles,
                   'cost': cost}
            report['records'].append(row)
            retained = report['retained_witnesses']
            if len(retained) < 40 or cost['rho'] < max(item['cost']['rho'] for item in retained):
                lattice.check_witness(values, targets)
                report['exact_full_replays'] += 1
                retained.append(row | {'seed': primitive, 'raw_values': values})
                report['retained_witnesses'] = sorted(retained, key=lambda item: item['cost']['rho'])[:40]
            if cost['rho'] < worker['best']-1e-10:
                worker['best'] = cost['rho']
                print(json.dumps({'event': 'new_best', 'label': worker['case']['label'], **cost,
                                  'degree': primitive['degree'], 'trial': trial}), flush=True)
            if cost['rho'] < 2**-.5-1e-8:
                stop = True
        if best_this_target is not None:
            worker['parents'].append((best_this_target, targets))
            worker['parents'] = sorted(set(worker['parents']))[:128]
        report['completed_trials'] += 1
        report['elapsed_seconds'] = time.monotonic()-began
        report['best_by_source'] = [worker['best'] for worker in workers]
        report['distinct_target_vectors'] = [len(worker['seen']) for worker in workers]
        report['distinct_profile_classes'] = [len(worker['profile_classes']) for worker in workers]
        atomic_write(args.output, report)
        if trial % 20 == 0:
            print(json.dumps({'event': 'progress', 'trial': trial, 'best': report['best_by_source'],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic()-began,
                  stopping_reason='strict numerical threshold; independent certificate needed' if stop else
                  ('declared wall-clock cap' if time.monotonic() >= deadline else 'complete trial list'))
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
