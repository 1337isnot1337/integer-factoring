#!/usr/bin/env python3
"""Finite F4 physical-projector screens with complete endpoint controls.

The family enumeration and local polynomial arithmetic are exact. Optimizing
the Gibbs root and scale is numerical, so these outputs are search evidence.
"""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import time

from physical_split_crossings import (
    SplitProfile, field_dot, hyperplane_split, identity_split, pair_split, retract,
)
from reflection_dual_rail_gaussian import physical_bit_energies, tensor


def physical_energies(interval):
    weights = physical_bit_energies(interval)
    return [sum(weight for bit, weight in enumerate(weights) if state >> bit & 1)
            for state in range(1 << len(weights))]


def cases(suite):
    if suite == 'hyperplanes':
        base, energy = tensor(4), physical_energies(4)
        for functional in product(range(4), repeat=4):
            if not any(functional) or next(value for value in functional if value) != 1:
                continue
            for complement in product(range(4), repeat=4):
                if field_dot(functional, complement) != 1:
                    continue
                alphabet, left, right = hyperplane_split(functional, complement)
                yield {'functional': functional, 'complement': complement}, \
                    retract(base, left, right), [([energy[i] for i in alphabet], None)]
        return

    base, energy = tensor(8), physical_energies(8)
    if suite in ('energy', 'asymmetric-energy'):
        for right_pair in [(0, 1)] + [(1, c) for c in range(4)]:
            for left_pair in product(range(4), repeat=2):
                if field_dot(left_pair, right_pair) != 1:
                    continue
                alphabet, left, right = pair_split(16, (6, 9), left_pair, right_pair)
                ex = [energy[i] for i in alphabet]
                if suite == 'energy':
                    projected, ey = retract(base, left, right), None
                else:
                    identity = identity_split(16)
                    projected = retract(base, left, right, identity, identity)
                    ey = energy
                yield {'pair': (6, 9), 'left': left_pair, 'right': right_pair}, \
                    projected, [(ex, ey)]
        return

    pairs = sorted(combinations(range(16), 2), key=lambda pair: abs(energy[pair[0]] - energy[pair[1]]))
    for pair in pairs:
        for left_pair, right_pair in [((1, 0), (1, 1)), ((0, 1), (1, 1)), ((2, 3), (1, 1))]:
            alphabet, left, right = pair_split(16, pair, left_pair, right_pair)
            energies = []
            for replacement in sorted({energy[pair[0]], energy[pair[1]],
                                       (energy[pair[0]] + energy[pair[1]]) / 2}):
                energies.append(([replacement if i == pair[0] else energy[i]
                                  for i in alphabet], None))
            yield {'pair': pair, 'left': left_pair, 'right': right_pair}, \
                retract(base, left, right), energies


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', required=True,
                        choices=('hyperplanes', 'energy', 'asymmetric-energy', 'pairs'))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    start = time.time()
    records, catalogue, by_profile = [], [], {}
    endpoint_entries = seeds = 0
    source_paths = [Path(__file__), Path(__file__).with_name('physical_split_crossings.py'),
                    Path(__file__).with_name('reflection_dual_rail_gaussian.py'),
                    Path(__file__).with_name('native_matrix_smith_f4.py'),
                    Path(__file__).with_name('native_matrix_smith_f4.cpp'),
                    Path(__file__).with_name('native_rectangular_smith_f4.cpp')]
    hashes = {path.name: sha256(path.read_bytes()).hexdigest() for path in source_paths}

    def save():
        best = min(range(len(catalogue)), key=lambda i: catalogue[i]['rho'])
        report = {
            'suite': args.suite, 'seed_count': seeds, 'height_cases': len(records),
            'endpoint_entries': endpoint_entries, 'profile_classes': len(catalogue),
            'best_class': best, 'classes': catalogue, 'cases': records,
            'elapsed_seconds': time.time() - start, 'source_sha256': hashes,
            'scope': 'exact finite split maps and complete local profiles at the displayed heights; exploratory numerical pressure only',
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({'seeds': seeds, 'height_cases': len(records),
                          'profile_classes': len(catalogue), 'best_rho': catalogue[best]['rho'],
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)

    for description, projected, energy_cases in cases(args.suite):
        oracle = SplitProfile(projected)
        seeds += 1
        endpoint_entries += len(projected['values'])
        for energy_x, energy_y in energy_cases:
            if energy_y is None:
                energy_y = energy_x
            profiles = oracle.profiles(energy_x, energy_y)
            for profile in profiles:
                assert all(abs(value - round(value)) < 1e-9 for value in profile)
            exact = [[round(value) for value in profile] for profile in profiles]
            key = json.dumps([energy_x, energy_y, exact], separators=(',', ':'))
            if key not in by_profile:
                numeric = oracle.pressure(energy_x, energy_y, exact)
                by_profile[key] = len(catalogue)
                catalogue.append({'energy_x': energy_x, 'energy_y': energy_y,
                                  'profiles': [dict(Counter(profile)) for profile in exact],
                                  'rho': numeric['rho'], 'log_scale': numeric['log_scale']})
            records.append({**description, 'class': by_profile[key]})
        if seeds % (100 if args.suite == 'hyperplanes' else 10) == 0:
            save()
    expected = {'hyperplanes': 5440, 'energy': 20, 'asymmetric-energy': 20, 'pairs': 360}
    assert seeds == expected[args.suite]
    save()


if __name__ == '__main__':
    main()
