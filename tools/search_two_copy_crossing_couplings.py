#!/usr/bin/env python3
"""A bounded finite search of polynomial couplings on 25 physical letters."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import random
import time

from search_f4_recentered_crossings import five_letter_source, root_for_profiles
from two_copy_crossing_couplings import CoupledCopies, constant_slab_tensors, endpoint_power


def compact(profiles):
    return [dict(sorted(Counter(side).items())) for side in profiles]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directions', type=Path, default=Path('evidence/source-slab-module-directions-h16.json'))
    parser.add_argument('--seconds', type=float, default=3600.)
    parser.add_argument('--maximum-cases', type=int)
    parser.add_argument('--random-seed', type=int, default=31102691)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic(); deadline = start+args.seconds
    directions = json.loads(args.directions.read_text())
    generators = [row for row in directions['minimal_basis']['generators'] if max(row) <= 1]
    assert len(generators) == 3
    supports = constant_slab_tensors(generators, 5)
    values, heights, ex, ey = five_letter_source()
    compiler = CoupledCopies(5, values, 18, heights)
    energies = ([a+b for a in ex for b in ex], [a+b for a in ey for b in ey])
    baseline_profiles = compiler.baseline_control()
    baseline = root_for_profiles(baseline_profiles, *energies)
    minimum = [min(side) for side in baseline_profiles]
    bounds = [[min(item[2] for item in compiler.weighted_support(row['support'], side))
               for side in (0, 1)] for row in supports]
    main_cases, extra_cases, seen = [], [], set()
    for i, j in combinations_with_replacement(range(len(supports)), 2):
        a, b = (bounds[i][side]+bounds[j][side] for side in (0, 1))
        lower = max(0, minimum[1]-b)
        top = min(36, a-minimum[0])
        span = top-lower
        if span <= 0:
            continue
        for power in sorted({span, *(r for r in (1, 2, 4, 8, 16, 32) if r <= span)}, reverse=True):
            for shift in sorted({lower, top-power}):
                for coefficient in (1, 2):
                    case = (i, j, shift, power, coefficient)
                    if case in seen:
                        continue
                    seen.add(case)
                    (main_cases if power == span else extra_cases).append(case)
    rng = random.Random(args.random_seed)
    rng.shuffle(main_cases); rng.shuffle(extra_cases)
    cases = main_cases+extra_cases
    if args.maximum_cases is not None:
        cases = cases[:args.maximum_cases]
    source_root = Path(__file__).resolve().parent
    source_names = (Path(__file__).name, 'two_copy_crossing_couplings.py',
                    'search_f4_recentered_crossings.py', 'native_matrix_smith_f4.py',
                    'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Finite exact-coefficient integer-gauge two-copy profiles and exploratory floating pressures; no factoring improvement certified.',
              'arguments': vars(args) | {'directions': str(args.directions), 'output': str(args.output)},
              'source_sha256': {name: sha256((source_root/name).read_bytes()).hexdigest() for name in source_names},
              'direction_record_sha256': sha256(args.directions.read_bytes()).hexdigest(),
              'source_alphabet': 5, 'physical_alphabet': 25, 'complete_degree_bound': 36,
              'source_heights_infinity_then_zero': heights, 'paired_energies': energies,
              'constant_slabs': generators, 'supports': supports, 'baseline': baseline,
              'baseline_profiles': compact(baseline_profiles), 'catalogue_size': len(cases),
              'maximal_endpoint_power_cases': len(main_cases), 'completed_cases': [],
              'profile_classes': [], 'unresolved_cases': [], 'completed_requested_schedule': False}
    signatures = {}; best = {'kind': 'uncoupled reference', **baseline, 'profiles': compact(baseline_profiles)}

    def save():
        report.update(elapsed_seconds=time.monotonic()-start, best=best,
                      strict_pressure_candidate_found=best['critical_pressure'] < -1e-9)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(report, indent=2)+'\n')
        temporary.replace(args.output)

    save()
    for index, (i, j, shift, power, coefficient) in enumerate(cases):
        if time.monotonic() >= deadline:
            report['time_budget_reached'] = True
            break
        recipe = {'first_support': i, 'second_support': j, 'monomial_shift': shift,
                  'endpoint_power': power, 'F4_coefficient': coefficient}
        polynomial = endpoint_power(power) << shift
        try:
            profiles = compiler.profiles(supports[i]['support'], supports[j]['support'], polynomial, coefficient)
        except (ArithmeticError, ValueError, OverflowError) as error:
            report['unresolved_cases'].append({'index': index, **recipe, 'reason': str(error)})
            save()
            continue
        signature = tuple(tuple(side) for side in profiles)
        if signature not in signatures:
            cost = root_for_profiles(profiles, *energies)
            signatures[signature] = len(report['profile_classes'])
            report['profile_classes'].append({'first_case': index, **recipe, **cost, 'profiles': compact(profiles)})
        profile_class = signatures[signature]
        cost = report['profile_classes'][profile_class]
        report['completed_cases'].append({'index': index, **recipe, 'profile_class': profile_class,
                                          'rho': cost['rho'], 'critical_pressure': cost['critical_pressure']})
        if cost['rho'] < best['rho']-1e-11:
            best = {'kind': 'coupled candidate', 'case': index, **cost}
        if index % 5 == 0:
            print(json.dumps({'case': index, 'classes': len(signatures), 'best': best['rho'],
                              'candidate': cost['rho'], 'seconds': time.monotonic()-start}), flush=True)
        save()
    else:
        report['completed_requested_schedule'] = True
    save()


if __name__ == '__main__':
    main()
