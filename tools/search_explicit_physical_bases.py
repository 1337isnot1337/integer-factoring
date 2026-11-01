#!/usr/bin/env python3
"""Bounded complete-basis search on explicit normalized crossing seeds.

Every constant physical basis and endpoint is exact over F4. Balanced
costs and positive supporting-line weights use complete native profiles
and a shared Gibbs-optimized pair of types. Floating results are search
records, not a factoring certificate.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_laurent_physical_frames import digest
from search_affine_coupled_reflection_frames import balanced
from search_fused_mixed_weighted import SearchDeadline, atomic_write
from search_reflection_physical_bases import (
    canonical_basis, change_basis, exact_controls, generate_basis_pair,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--trial-seconds', type=float, default=10)
    parser.add_argument('--iterations', type=int, default=180)
    parser.add_argument('--cases', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=312541)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    source = json.loads(raw)
    rng, noise = random.Random(args.seed), np.random.default_rng(args.seed)
    controls, seeds = [], source['cases']
    for case in seeds:
        seed = case['seed']
        assert seed['q_x'] == seed['q_u'] == 4
        assert digest(seed['values']) == case['tensor_sha256']
        params = np.asarray(case['initial_params'], float)
        value = balanced(AsymmetricPressure(seed, 1), params)[0]
        assert abs(value - case['initial_balanced_rho']) < 2e-9
        controls.append({'label': case['label'], 'initial_balanced_rho': value,
                         'exact_profile_controls': exact_controls(seed)})
    names = (Path(__file__).name, 'search_reflection_physical_bases.py',
             'search_affine_coupled_reflection_frames.py', 'asymmetric_crossing_pressure.py',
             'physical_split_crossings.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'scope': __doc__, 'arguments': vars(args) | {'input': str(args.input), 'output': str(args.output)},
              'input_sha256': sha256(raw).hexdigest(), 'initial_controls': controls,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'records': [], 'winners': {}, 'completed': False, 'completed_trials': 0,
              'counters': {'objective_calls': 0, 'unresolved_calls': 0, 'duplicate_bases': 0}}
    elites = {(index, mode): [] for index in range(len(seeds)) for mode in ('balanced', 'asymmetric')}
    seen = set()
    began = time.monotonic()
    deadline = began + args.seconds
    generation = [0] * len(seeds)
    attempts = 0
    while attempts < args.cases and time.monotonic() < deadline:
        index = attempts % len(seeds)
        case = seeds[index]
        number = generation[index]
        generation[index] += 1
        mode = 'asymmetric' if number % 3 == 2 else 'balanced'
        kind, left, right = generate_basis_pair(number, 4, rng, elites[index, mode])
        key = index, mode, canonical_basis(left), canonical_basis(right)
        attempts += 1
        if key in seen:
            report['counters']['duplicate_bases'] += 1
            continue
        seen.add(key)
        changed = change_basis(case['seed'], left, right)
        assert changed['degree'] == case['seed']['degree']
        oracle = AsymmetricPressure(changed, 1)
        trial_started = time.monotonic()
        trial_deadline = min(deadline, trial_started + args.trial_seconds)
        warm = np.asarray(case['initial_params'], float)
        balanced_elites = elites[index, 'balanced']
        if balanced_elites and number > 7 and number % 4:
            warm = np.asarray(balanced_elites[0]['best']['parameters'], float)
        starts = [warm]
        if number % 5 == 0:
            perturbed = warm.copy()
            perturbed[:-1] += noise.normal(0, .06, len(warm) - 1)
            perturbed[-1] *= .7 if number % 2 else 1.3
            starts.append(perturbed)
        best, outcomes, completed, error = None, [], True, None

        def objective(values):
            nonlocal best
            if time.monotonic() >= trial_deadline:
                raise SearchDeadline()
            report['counters']['objective_calls'] += 1
            try:
                if mode == 'balanced':
                    value, gradient, details = balanced(oracle, values)
                else:
                    value, gradient, details = oracle.evaluate(values)
                if best is None or value < best['value']:
                    best = {'value': float(value), 'parameters': list(map(float, values)),
                            'type_details': details}
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved_calls'] += 1
                return 1e6, np.zeros(len(values))

        try:
            for start in starts:
                angles = (math.pi / 4, .25, 1.32) if mode == 'asymmetric' else (None,)
                for angle in angles:
                    initial = start if angle is None else np.r_[start / math.sqrt(2), angle]
                    bounds = [(-32, 32)] * 24 + [(0.0005, 2)]
                    if angle is not None:
                        bounds.append((.03, math.pi / 2 - .03))
                    objective(initial)
                    fit = minimize(objective, initial, jac=True, method='SLSQP', bounds=bounds,
                                   options={'maxiter': args.iterations, 'ftol': 2e-11})
                    outcomes.append({'start_angle': angle, 'iterations': int(fit.nit),
                                     'evaluations': int(fit.nfev), 'success': bool(fit.success),
                                     'message': str(fit.message)})
        except SearchDeadline:
            completed, error = False, 'Declared trial or global time cap reached.'
        record = {'source_index': index, 'label': case['label'], 'generation': number,
                  'mode': mode, 'kind': kind, 'basis_x': left, 'basis_y': right,
                  'tensor_sha256': digest(changed['values']), 'degree': changed['degree'],
                  'endpoint_entries': 256, 'completed': completed, 'interruption': error,
                  'seconds': time.monotonic() - trial_started, 'outcomes': outcomes, 'best': best}
        if best is not None:
            winner_key = f'{index}:{mode}'
            if winner_key not in report['winners'] or best['value'] < report['winners'][winner_key]['best']['value']:
                report['winners'][winner_key] = record
            elites[index, mode].append(record)
            elites[index, mode].sort(key=lambda row: row['best']['value'])
            del elites[index, mode][12:]
        report['records'].append(record)
        report.update(elapsed_seconds=time.monotonic() - began, attempted_cases=attempts,
                      completed_trials=sum(row['completed'] for row in report['records']))
        atomic_write(args.output, report)
        print(json.dumps({'label': case['label'], 'mode': mode, 'generation': number,
                          'value': None if best is None else best['value'],
                          'best': report['winners'].get(f'{index}:{mode}', {}).get('best', {}).get('value'),
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        if best is not None and best['value'] < (1 / math.sqrt(2) if mode == 'balanced' else 1) - 1e-7:
            report['stopping_reason'] = 'Numerical threshold candidate requires exact independent replay.'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  stopping_reason=report.get('stopping_reason', 'Declared time or case budget completed.'))
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
