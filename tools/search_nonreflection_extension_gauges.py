#!/usr/bin/env python3
"""Apply positive-mass full gauges to the best nonreflection q6 extensions.

This is a separate, explicitly selected schedule. The earlier top-six
natural selection consists entirely of reflection-compatible states.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from balanced_particle_extensions import construct
from search_balanced_extension_gauges import optimize
from search_f4_recentered_crossings import digest


def reflection_admissible(mask):
    return all(((mask >> i) & 1)+((mask >> (7-i)) & 1) == 1 for i in range(4))


def select(record, per_family):
    rows = []
    for family in ('symmetric_6x6', 'rectangular_5x6', 'rectangular_6x5'):
        candidates = sorted((row for row in record['results'] if row['family'] == family
                             and not reflection_admissible(row['added_mask'])), key=lambda row: row['rho'])
        seen = set()
        for candidate in candidates:
            if candidate['digest'] in seen:
                continue
            seen.add(candidate['digest'])
            rows.append(candidate)
            if len(seen) == per_family:
                break
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--per-family', type=int, default=3)
    parser.add_argument('--iterations', type=int, default=250)
    parser.add_argument('--floors', type=float, nargs='+', default=[0., .001, .01])
    parser.add_argument('--height-bound', type=float, default=30.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.per_family < 1 or args.iterations < 1 or any(not 0 <= value < .5 for value in args.floors):
        parser.error('Use positive schedule sizes and mass floors in[0,.5).')
    source = json.loads(args.source.read_text())
    assert source['completed_full_extension_schedule']
    selected = select(source, args.per_family)
    started = time.monotonic()
    names = (Path(__file__).name, 'search_balanced_extension_gauges.py',
             'balanced_particle_extensions.py', 'bounded_paired_gibbs.py',
             'balanced_four_particle_codes.py', 'selected_gaussian_coefficients.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'search_f4_recentered_crossings.py')
    report = {'scope': 'Finite complete physical-gauge optimization of selected nonreflection balanced extensions, including positive added-letter type lower bounds. No minimum or exact pressure-sign certificate.',
              'arguments': vars(args) | {'source': str(args.source), 'output': str(args.output)},
              'input_sha256': {str(args.source): sha256(args.source.read_bytes()).hexdigest()},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'planned_candidates': len(selected), 'candidates': []}
    for row in selected:
        seed = construct(row['code_x'], row['code_y'])
        assert digest(seed['values']) == row['digest']
        runs = optimize(seed, row, source['results'][0], args.floors, args.iterations, args.height_bound)
        valid = [run['best'] for run in runs if run['best'] is not None]
        assert valid
        best = min(valid, key=lambda item: item['rho'])
        report['candidates'].append({'source_row': row, 'literal_seed': seed, 'runs': runs, 'best': best})
        report.update(completed_requested_schedule=len(report['candidates']) == len(selected),
                      best_recorded_rho=min(item['best']['rho'] for item in report['candidates']),
                      elapsed_seconds=time.monotonic()-started)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'event': 'candidate', 'completed': len(report['candidates']),
                          'planned': len(selected), 'family': row['family'], 'added_mask': row['added_mask'],
                          'best': best['rho'], 'natural': row['rho'],
                          'queries': sum(run['profile_queries'] for run in runs),
                          'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
