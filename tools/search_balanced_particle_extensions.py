#!/usr/bin/env python3
"""Screen every one-state extension of the retained balanced five-state code."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from balanced_particle_extensions import complete_profiles, construct, extension_cases
from search_f4_recentered_crossings import digest, root_for_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queries', type=int, default=196)
    parser.add_argument('--seconds', type=float, default=1800.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.queries < 1 or args.seconds <= 0:
        parser.error('Use positive query and time limits.')
    started = time.monotonic()
    names = (Path(__file__).name, 'balanced_particle_extensions.py', 'balanced_four_particle_codes.py',
             'selected_gaussian_coefficients.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'search_f4_recentered_crossings.py')
    report = {'scope': 'All65 one-state balanced extensions of the retained q5 code, with symmetric6x6 and both rectangular5x6 orientations. Complete natural integer-height profiles and numerical Gibbs/scale costs.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'planned_cases': 196, 'results': []}
    cache = {}
    for family, added, code_x, code_y in extension_cases():
        if len(report['results']) >= args.queries or time.monotonic()-started >= args.seconds:
            break
        seed = construct(code_x, code_y)
        profiles = complete_profiles(seed)
        key = profiles, seed['energies_x'], seed['energies_y']
        if key not in cache:
            cache[key] = root_for_profiles(profiles, seed['energies_x'], seed['energies_y'])
        cost = cache[key]
        row = {'family': family, 'added_mask': added, 'code_x': code_x, 'code_y': code_y,
               'dimensions': seed['dimensions'], 'degree': seed['degree'],
               'entry_gcd_in_z': seed['entry_gcd_in_z'], 'digest': digest(seed['values']),
               'energies_x': seed['energies_x'], 'energies_y': seed['energies_y'],
               'profiles': profiles, **cost}
        report['results'].append(row)
        report.update(completed_requested_schedule=len(report['results']) == min(args.queries,196),
                      completed_full_extension_schedule=len(report['results']) == 196,
                      complete_tensor_queries=len(report['results']),
                      endpoint_entries_checked=sum(item['dimensions'][0]**2*item['dimensions'][1]**2 for item in report['results']),
                      best_recorded_rho=min(item['rho'] for item in report['results']),
                      distinct_energy_profile_classes=len(cache), elapsed_seconds=time.monotonic()-started)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'event': 'extension', 'query': len(report['results']),
                          'family': family, 'added_mask': added, 'rho': cost['rho'],
                          'best': report['best_recorded_rho'], 'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
