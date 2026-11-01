#!/usr/bin/env python3
"""All physical gauges for selected q6 and rectangular balanced extensions.

Both unconstrained and positively bounded added-letter types are tested.
The exact tensors are retained literally; all gauge and type optimization
is finite floating-point search evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from balanced_particle_extensions import construct, start_params
from bounded_paired_gibbs import paired_value_gradient
from native_rectangular_smith_f4 import RectangularOracle
from search_f4_recentered_crossings import digest


def cases(record, per_family):
    result = [record['results'][0]]
    for family in ('symmetric_6x6', 'rectangular_5x6', 'rectangular_6x5'):
        rows = sorted((row for row in record['results'] if row['family'] == family), key=lambda row: row['rho'])
        seen = set()
        for row in rows:
            if row['digest'] in seen:
                continue
            seen.add(row['digest'])
            result.append(row)
            if len(seen) == per_family:
                break
    return result


def lower_bounds(seed, added_mask, mass):
    return tuple(tuple(mass if mask == added_mask else 0. for mask in seed[key])
                 for key in ('code_x', 'code_y'))


def interior_start(seed, reference, added, target_mass):
    """Embed the q5 natural heights and give each added symbol target Gibbs mass."""
    scale, rate = reference['scale'], 2*reference['rho']
    energies = []
    for key, energy_key in (('code_x', 'energies_x'), ('code_y', 'energies_y')):
        old = dict(zip(reference[key], reference[energy_key]))
        partition = sum(math.exp(-scale*value/rate) for value in old.values())
        extra = rate/scale*math.log((1-target_mass)/(target_mass*partition))
        energies.append([old.get(mask, extra) for mask in seed[key]])
    return start_params(seed, scale, *energies)


def optimize(seed, source, reference, floors, maximum_iterations, height_bound):
    import numpy as np
    from scipy.optimize import minimize
    a, b = seed['dimensions'][:2]
    oracle = RectangularOracle(seed['values'], a, b)
    starts = [('natural', start_params(seed, source['scale'])),
              ('interior-warm', interior_start(seed, reference, source['added_mask'], .02)),
              ('uniform', start_params(seed, .08, [0.]*a, [0.]*b))]
    rows = []
    for mass in floors:
        lower_x, lower_y = lower_bounds(seed, source['added_mask'], mass)
        for start_name, start in starts:
            best = None
            attempts = 0
            failures = Counter()

            def objective(parameters):
                nonlocal best, attempts
                attempts += 1
                try:
                    value, gradient, px, py, orders, gx, gy = paired_value_gradient(
                        oracle, parameters, lower_x, lower_y)
                    assert all(value+1e-12 >= bound for value, bound in zip(px, lower_x))
                    assert all(value+1e-12 >= bound for value, bound in zip(py, lower_y))
                    if best is None or value < best['rho']:
                        best = {'rho': value, 'parameters': list(map(float, parameters)),
                                'type_x': px, 'type_y': py, 'orders_times_log_s': orders,
                                'energy_x_times_log_s': gx, 'energy_y_times_log_s': gy,
                                'query': attempts}
                    return value, gradient
                except (ArithmeticError, ValueError, OverflowError) as error:
                    failures[str(error)] += 1
                    return 1e6, np.zeros(oracle.parameters)

            result = minimize(objective, start, method='SLSQP', jac=True,
                              bounds=[(-height_bound,height_bound)]*(oracle.parameters-1)+[(.0001,2.)],
                              options={'maxiter': maximum_iterations, 'ftol': 1e-11})
            objective(result.x)
            rows.append({'start': start_name, 'minimum_added_mass': mass,
                         'lower_x': lower_x, 'lower_y': lower_y, 'best': best,
                         'profile_queries': attempts, 'failed_profiles': dict(failures),
                         'optimizer_success': bool(result.success), 'optimizer_message': str(result.message),
                         'optimizer_evaluations': int(result.nfev)})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--per-family', type=int, default=6)
    parser.add_argument('--iterations', type=int, default=250)
    parser.add_argument('--floors', type=float, nargs='+', default=[0., .001, .01])
    parser.add_argument('--height-bound', type=float, default=30.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.per_family < 1 or args.iterations < 1 or any(not 0 <= value < .5 for value in args.floors):
        parser.error('Use positive schedule sizes and added-letter mass floors in[0,.5).')
    source = json.loads(args.source.read_text())
    assert source['completed_full_extension_schedule']
    selected = cases(source, args.per_family)
    started = time.monotonic()
    names = (Path(__file__).name, 'balanced_particle_extensions.py', 'bounded_paired_gibbs.py',
             'balanced_four_particle_codes.py', 'selected_gaussian_coefficients.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'search_f4_recentered_crossings.py')
    report = {'scope': 'Finite complete physical-gauge optimization of selected balanced extensions, including constrained positive added-letter type mass. No minimum or exact pressure-sign certificate.',
              'arguments': vars(args) | {'source': str(args.source), 'output': str(args.output)},
              'input_sha256': {str(args.source): sha256(args.source.read_bytes()).hexdigest()},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'planned_candidates': len(selected), 'candidates': []}
    for row in selected:
        seed = construct(row['code_x'], row['code_y'])
        assert digest(seed['values']) == row['digest']
        floors = (0.,) if row['added_mask'] is None else args.floors
        runs = optimize(seed, row, source['results'][0], floors, args.iterations, args.height_bound)
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
