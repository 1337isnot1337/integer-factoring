#!/usr/bin/env python3
"""Finite exact-tensor screening of unrestricted balanced q5 physical codes.

Profiles use complete finite-field polynomial matrices.  The scale/root
minimization remains exploratory floating point; no minimum is certified.
The complete deterministic and random code schedule is stored for replay.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import time

from balanced_four_particle_codes import (FOUR_PARTICLE_STATES, REFERENCE_CODE,
                                         canonical_code, construct, grade)
from reflection_dual_rail_gaussian import reflection_letters
from search_f4_recentered_crossings import digest, root_for_profiles
from search_saturated_dual_directions import complete_profiles


def schedule(low_states, random_codes, seed):
    rng = random.Random(seed)
    result, seen = [], set()

    def add(code, stage):
        code = canonical_code(code)
        if code not in seen:
            seen.add(code)
            result.append((code, stage))

    add(REFERENCE_CODE, 'reference')
    for removed in REFERENCE_CODE:
        kept = tuple(mask for mask in REFERENCE_CODE if mask != removed)
        for replacement in FOUR_PARTICLE_STATES:
            if replacement not in REFERENCE_CODE:
                add(kept+(replacement,), 'one_symbol_replacement')
    ordered = sorted(FOUR_PARTICLE_STATES, key=lambda mask: (grade(mask), mask))
    for code in combinations(ordered[:low_states], 5):
        add(code, 'complete_low_grade_shell')
    for trial in range(random_codes):
        # Half the choices concentrate on short low-grade windows, where
        # repeated minimum energies can occur without a large spread.
        if trial % 2:
            width = rng.randrange(7, 19)
            start = rng.randrange(len(ordered)-width+1)
            code = rng.sample(ordered[start:start+width], 5)
        else:
            code = rng.sample(ordered, 5)
        add(code, 'random_code')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--low-states', type=int, default=12)
    parser.add_argument('--random-codes', type=int, default=1000)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=2026103208)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 5 <= args.low_states <= 18 or args.random_codes < 0 or args.seconds <= 0:
        parser.error('Use 5..18 low states, a nonnegative random count, and a positive time budget.')
    started = time.monotonic()
    planned = schedule(args.low_states, args.random_codes, args.seed)
    cache, results, retained = {}, [], []
    reflection = set(reflection_letters(8))
    best = float('inf')
    names = (Path(__file__).name, 'balanced_four_particle_codes.py', 'selected_gaussian_coefficients.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'saturated_dual_kernel_directions.py', 'search_saturated_dual_directions.py',
             'search_f4_recentered_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'scope': 'Five-state strict Gaussian codes from all seventy four-particle masks; natural paired energies only; numerical cost screening.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'planned_schedule': planned, 'alphabet': 5, 'available_physical_states': FOUR_PARTICLE_STATES,
              'coefficient_field': 'F2', 'common_entry_gcd_removed_exactly': True}

    def save(completed=False):
        report.update(results=results, retained=retained, completed_requested_schedule=completed,
                      complete_tensor_queries=len(results), exact_endpoint_entries=625*len(results),
                      distinct_energy_profile_classes=len(cache), best_recorded_rho=best,
                      strict_pressure_candidate_found=any(row['critical_pressure'] < -1e-9 for row in results),
                      elapsed_seconds=time.monotonic()-started)
        args.output.write_text(json.dumps(report, indent=2)+'\n')

    for code, stage in planned:
        if time.monotonic()-started >= args.seconds:
            break
        tensor = construct(code)
        energy = tensor['paired_energies']
        profiles = complete_profiles(tensor['values'], {'natural': energy})['natural']
        key = energy, profiles
        if key not in cache:
            cache[key] = root_for_profiles(profiles, energy, energy)
        cost = cache[key]
        row = {key: value for key, value in tensor.items() if key != 'values'}
        row.update(stage=stage, profiles=profiles, digest=digest(tensor['values']),
                   nonreflection_states=sum(mask not in reflection for mask in code), **cost)
        results.append(row)
        if cost['rho'] < best-1e-10:
            best = cost['rho']
            print(json.dumps({'event': 'new_best', 'query': len(results), 'code': code,
                              'degree': tensor['degree'], 'energies': energy, **cost}), flush=True)
        if len(retained) < 24 or cost['rho'] < max(item['rho'] for item in retained):
            retained.append(row | {'values': tensor['values']})
            retained[:] = sorted(retained, key=lambda item: item['rho'])[:24]
        if len(results) % 100 == 0:
            save()
            print(json.dumps({'event': 'checkpoint', 'queries': len(results), 'planned': len(planned),
                              'profile_classes': len(cache), 'best_rho': best,
                              'seconds': time.monotonic()-started}), flush=True)
    save(len(results) == len(planned))
    print(json.dumps({'event': 'finished', 'queries': len(results), 'planned': len(planned),
                      'completed': len(results) == len(planned), 'best_rho': best,
                      'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
