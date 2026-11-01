#!/usr/bin/env python3
"""Exact finite screening of distinct paired four-particle alphabets."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from asymmetric_four_particle_codes import complete_profiles, construct
from balanced_four_particle_codes import FOUR_PARTICLE_STATES, REFERENCE_CODE, canonical_code
from reflection_dual_rail_gaussian import reflection_letters
from search_f4_recentered_crossings import digest, root_for_profiles


def schedule():
    reference = canonical_code(REFERENCE_CODE)
    result = [(reference, reference, 'reference')]
    for removed in reference:
        kept = tuple(mask for mask in reference if mask != removed)
        for replacement in FOUR_PARTICLE_STATES:
            if replacement not in reference:
                result.append((reference, canonical_code(kept + (replacement,)), 'one_paired_symbol_replacement'))
    assert len(result) == len(set((x, y) for x, y, _ in result)) == 326
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    planned = schedule()
    reflection = set(reflection_letters(8))
    cache, rows, retained = {}, [], []
    best = float('inf')
    names = (Path(__file__).name, 'asymmetric_four_particle_codes.py', 'balanced_four_particle_codes.py',
             'selected_gaussian_coefficients.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'saturated_dual_kernel_directions.py',
             'search_f4_recentered_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'scope': 'All325 one-symbol replacements on Y/U with the retained X/V alphabet fixed; exact tensors and complete profiles, exploratory numerical physical-type costs.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'planned_schedule': planned, 'alphabet': 5}

    def save():
        report.update(results=rows, retained=retained, complete_tensor_queries=len(rows),
                      completed_requested_schedule=len(rows) == len(planned),
                      distinct_energy_profile_classes=len(cache),
                      best_recorded_rho=best, elapsed_seconds=time.monotonic() - started,
                      strict_pressure_candidate_found=any(row['critical_pressure'] < -1e-9 for row in rows))
        args.output.write_text(json.dumps(report, indent=2) + '\n')

    for code_x, code_y, stage in planned:
        if time.monotonic() - started >= args.seconds:
            break
        tensor = construct(code_x, code_y)
        profiles = complete_profiles(tensor['values'], tensor['energies_x'], tensor['energies_y'])
        key = tensor['energies_x'], tensor['energies_y'], profiles
        if key not in cache:
            cache[key] = root_for_profiles(profiles, tensor['energies_x'], tensor['energies_y'])
        cost = cache[key]
        row = {key: value for key, value in tensor.items() if key != 'values'}
        row.update(stage=stage, profiles=profiles, digest=digest(tensor['values']),
                   nonreflection_y_states=sum(mask not in reflection for mask in code_y), **cost)
        rows.append(row)
        if cost['rho'] < best:
            best = cost['rho']
            print(json.dumps({'event': 'new_best', 'query': len(rows), 'code_y': code_y, **cost}), flush=True)
        if len(retained) < 24 or cost['rho'] < max(item['rho'] for item in retained):
            retained.append(row | {'values': tensor['values']})
            retained.sort(key=lambda item: item['rho'])
            del retained[24:]
        if len(rows) % 50 == 0:
            save()
            print(json.dumps({'event': 'progress', 'queries': len(rows), 'best': best,
                              'seconds': time.monotonic() - started}), flush=True)
    save()
    print(json.dumps({'event': 'finished', 'queries': len(rows), 'best': best,
                      'seconds': time.monotonic() - started}), flush=True)


if __name__ == '__main__':
    main()
