#!/usr/bin/env python3
"""Optimize all physical heights for selected nonreflection Gaussian codes.

Complete polynomial tensors and field operations are exact. Valuation
weights, scale, Gibbs types, and SLSQP optimization are numerical search
evidence. A finite optimizer run is not a minimum certificate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from asymmetric_four_particle_codes import construct
from search_f4_recentered_crossings import digest


def optimize_observed(seed, starts, maximum_iterations):
    """Retain only candidates whose complete numerical oracle call succeeded.

    A failed trial can have a smaller arbitrary penalty than a valid large
    starting cost. Never turn that optimization penalty into a witness.
    """
    import numpy as np
    from scipy.optimize import minimize
    from native_rectangular_smith_f4 import RectangularOracle
    oracle = RectangularOracle(seed['values'], 5, 5)
    records = []
    for name, start in starts:
        observed = None
        initial_rho = None
        errors = Counter()

        def objective(params):
            nonlocal observed, initial_rho
            try:
                rho, px, py, orders, gx, gy = oracle.types(params)
                entropy = -sum(value*math.log(value) for value in px+py if value > 0)
                if entropy < 1e-9:
                    raise ArithmeticError('The paired type entropy is too small for a stable gradient.')
                if initial_rho is None:
                    initial_rho = rho
                if observed is None or rho < observed['rho']:
                    observed = {'rho': rho, 'params': list(map(float, params)), 'type_x': px, 'type_y': py,
                                'energy_x_times_log_s': gx, 'energy_y_times_log_s': gy,
                                'orders_times_log_s': orders}
                _, gradient = oracle.evaluate(params, px, py)
                return rho, np.asarray(gradient)/(2*entropy)
            except (ArithmeticError, ValueError) as error:
                errors[str(error)] += 1
                return 1e6, np.zeros(len(params))

        initial = objective(start)[0]
        result = minimize(objective, start, jac=True, method='SLSQP',
                          bounds=[(-20, 20)]*32+[(.0001, 2)],
                          options={'maxiter': maximum_iterations, 'ftol': 1e-11})
        # Inspect the final iterate, but choose the best successful full
        # profile evaluation seen anywhere in the run.
        objective(result.x)
        if observed is None:
            records.append({'start': name, 'status': 'no_successful_complete_profile',
                            'errors': dict(errors), 'evaluations': result.nfev,
                            'optimizer_success': bool(result.success), 'optimizer_message': str(result.message)})
            continue
        records.append(observed | {'start': name, 'initial_rho': initial_rho,
                       'initial_objective': initial, 'status': 'successful_profile_retained',
                       'evaluations': result.nfev, 'errors': dict(errors),
                       'optimizer_success': bool(result.success), 'optimizer_message': str(result.message),
                       'optimizer_reported_objective': float(result.fun),
                       'used_initial_candidate': observed['params'] == list(map(float, start))})
    valid = [row for row in records if row['status'] == 'successful_profile_retained']
    if not valid:
        raise ArithmeticError('No starting family yielded a complete numerical profile.')
    return min(valid, key=lambda row: row['rho']), records


def start_params(ex, ey, scale):
    import numpy as np
    params = np.zeros(33)
    vectors = {(0, 2): ey, (0, 3): ex, (1, 0): ex, (1, 1): ey}
    for (side, leg), energy in vectors.items():
        for symbol in range(1, 5):
            params[(4*side+leg)*4+symbol-1] = (energy[symbol]-energy[0])*scale/2
    params[-1] = scale
    return params


def selected_cases(symmetric, asymmetric, per_family):
    sources = [('symmetric', symmetric, 'nonreflection_states'),
               ('asymmetric', asymmetric, 'nonreflection_y_states')]
    cases = [('reference', asymmetric['results'][0])]
    seen = {cases[0][1]['digest']}
    for family, record, key in sources:
        count = 0
        for row in sorted(record['results'], key=lambda item: item['rho']):
            if not row[key] or row['digest'] in seen:
                continue
            seen.add(row['digest'])
            cases.append((family, row))
            count += 1
            if count == per_family:
                break
    return cases


def main():
    import numpy as np
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--symmetric', type=Path, default=Path('evidence/balanced-four-particle-code-search.json'))
    parser.add_argument('--asymmetric', type=Path, default=Path('evidence/asymmetric-four-particle-codes-q5.json'))
    parser.add_argument('--per-family', type=int, default=4)
    parser.add_argument('--iterations', type=int, default=250)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.per_family < 1 or args.iterations < 1:
        parser.error('Use positive candidate and iteration counts.')
    started = time.monotonic()
    sym, asym = json.loads(args.symmetric.read_text()), json.loads(args.asymmetric.read_text())
    assert sym['completed_requested_schedule'] and asym['completed_requested_schedule']
    cases = selected_cases(sym, asym, args.per_family)
    names = (Path(__file__).name, 'asymmetric_four_particle_codes.py', 'balanced_four_particle_codes.py',
             'selected_gaussian_coefficients.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'search_f4_recentered_crossings.py')
    report = {'scope': 'Finite all33-variable physical-height and scale optimization of selected exact physical-code seeds; no minimum or threshold certificate.',
              'arguments': vars(args) | {'symmetric': str(args.symmetric), 'asymmetric': str(args.asymmetric), 'output': str(args.output)},
              'input_sha256': {str(path): sha256(path.read_bytes()).hexdigest() for path in (args.symmetric, args.asymmetric)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'planned_candidates': len(cases), 'candidates': []}
    for family, initial in cases:
        x = initial.get('code_x', initial.get('code_masks'))
        y = initial.get('code_y', initial.get('code_masks'))
        seed = construct(x, y)
        assert digest(seed['values']) == initial['digest']
        seed['q_x'] = seed['q_u'] = 5
        ex, ey = seed['energies_x'], seed['energies_y']
        scale = initial['scale']
        target_x = np.zeros(5)
        target_y = np.zeros(5)
        target_x[np.argsort(ex)] = (0, 0, 2.1, 2.6, 4.1)
        target_y[np.argsort(ey)] = (0, 0, 2.1, 2.6, 4.1)
        target = start_params(target_x, target_y, 1.)
        target[-1] = .12
        starts = [('recorded-natural', start_params(ex, ey, scale)),
                  ('two-zero-target', target),
                  ('uniform', start_params((0,)*5, (0,)*5, .08))]
        best, runs = optimize_observed(seed, starts, args.iterations)
        report['candidates'].append({'family': family, 'source_row': initial,
                                     'literal_seed': seed, 'best': best, 'runs': runs})
        report.update(completed_requested_schedule=len(report['candidates']) == len(cases),
                      best_recorded_rho=min(row['best']['rho'] for row in report['candidates']),
                      elapsed_seconds=time.monotonic()-started)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'event': 'candidate', 'family': family, 'code_x': x, 'code_y': y,
                          'initial_rho': initial['rho'], 'best_rho': best['rho'],
                          'evaluations': sum(run['evaluations'] for run in runs),
                          'optimizer_successes': sum(run['optimizer_success'] for run in runs),
                          'completed': len(report['candidates']), 'planned': len(cases),
                          'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
