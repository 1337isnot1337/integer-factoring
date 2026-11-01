#!/usr/bin/env python3
"""Complete low-grade Y/U shell with the retained X/V alphabet fixed."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

from asymmetric_four_particle_codes import complete_profiles, construct
from balanced_four_particle_codes import FOUR_PARTICLE_STATES, REFERENCE_CODE, canonical_code, grade
from reflection_dual_rail_gaussian import reflection_letters
from search_f4_recentered_crossings import digest, root_for_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--states', type=int, default=16)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--reuse', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 5 <= args.states <= 20:
        parser.error('Use5..20 lowest-grade physical states.')
    started = time.monotonic()
    code_x = canonical_code(REFERENCE_CODE)
    ordered = sorted(FOUR_PARTICLE_STATES, key=lambda mask: (grade(mask), mask))[:args.states]
    planned = [canonical_code(code) for code in combinations(ordered, 5)]
    prior = json.loads(args.reuse.read_text())
    assert prior['completed_requested_schedule']
    prior_rows = {tuple(row['code_y']): row for row in prior['results'] if tuple(row['code_x']) == code_x}
    reused = [prior_rows[code] for code in planned if code in prior_rows]
    fresh = [code for code in planned if code not in prior_rows]
    cache, rows, retained = {}, [], []
    best = min(row['rho'] for row in reused)
    reflection = set(reflection_letters(8))
    names = (Path(__file__).name, 'asymmetric_four_particle_codes.py', 'balanced_four_particle_codes.py',
             'selected_gaussian_coefficients.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'saturated_dual_kernel_directions.py',
             'search_f4_recentered_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'scope': 'Complete low-grade asymmetric physical-code shell; exact tensors and full profiles, numerical cost screening.',
              'arguments': vars(args) | {'reuse': str(args.reuse), 'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'reused_record_sha256': sha256(args.reuse.read_bytes()).hexdigest(),
              'fixed_code_x': code_x, 'available_shell_states': ordered,
              'complete_shell_schedule': planned, 'reused_results': reused,
              'fresh_schedule': fresh}

    def save():
        report.update(results=rows, retained=retained, complete_new_tensor_queries=len(rows),
                      completed_requested_schedule=len(rows) == len(fresh),
                      covered_shell_codes=len(reused) + len(rows),
                      best_recorded_rho=best, elapsed_seconds=time.monotonic() - started,
                      strict_pressure_candidate_found=any(row['critical_pressure'] < -1e-9 for row in reused + rows))
        args.output.write_text(json.dumps(report, indent=2) + '\n')

    for code_y in fresh:
        if time.monotonic() - started >= args.seconds:
            break
        tensor = construct(code_x, code_y)
        profiles = complete_profiles(tensor['values'], tensor['energies_x'], tensor['energies_y'])
        key = tensor['energies_x'], tensor['energies_y'], profiles
        if key not in cache:
            cache[key] = root_for_profiles(profiles, tensor['energies_x'], tensor['energies_y'])
        cost = cache[key]
        row = {key: value for key, value in tensor.items() if key != 'values'}
        row.update(profiles=profiles, digest=digest(tensor['values']),
                   nonreflection_y_states=sum(mask not in reflection for mask in code_y), **cost)
        rows.append(row)
        if cost['rho'] < best:
            best = cost['rho']
            print(json.dumps({'event': 'new_best', 'query': len(rows), 'code_y': code_y, **cost}), flush=True)
        if len(retained) < 24 or cost['rho'] < max(item['rho'] for item in retained):
            retained.append(row | {'values': tensor['values']})
            retained.sort(key=lambda item: item['rho'])
            del retained[24:]
        if len(rows) % 200 == 0:
            save()
            print(json.dumps({'event': 'progress', 'new_queries': len(rows),
                              'covered': len(reused) + len(rows), 'total': len(planned),
                              'best': best, 'seconds': time.monotonic() - started}), flush=True)
    save()
    print(json.dumps({'event': 'finished', 'new_queries': len(rows),
                      'covered': len(reused) + len(rows), 'best': best,
                      'seconds': time.monotonic() - started}), flush=True)


if __name__ == '__main__':
    main()
