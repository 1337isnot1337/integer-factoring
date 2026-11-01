#!/usr/bin/env python3
"""Budgeted full-cost supporting-line search on explicit mixed-place seeds."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from affine_place_profile_oracle import AffinePlaceOracle, primitive_seed
from asymmetric_mixed_place_pressure import ScaledMixedSupportingPressure
from search_fused_mixed_weighted import SearchDeadline, atomic_write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--trial-seconds', type=float, default=15)
    parser.add_argument('--iterations', type=int, default=180)
    parser.add_argument('--starts', type=int, default=3)
    parser.add_argument('--seed', type=int, default=312443)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    source = json.loads(raw)
    rng = np.random.default_rng(args.seed)
    names = (Path(__file__).name, 'asymmetric_mixed_place_pressure.py',
             'affine_place_profile_oracle.py', 'fused_physical_basis_oracle.py',
             'search_fused_mixed_weighted.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Bounded floating complete-profile supporting-line search; no global optimum or factoring certificate.',
              'arguments': vars(args) | {'input': str(args.input), 'output': str(args.output)},
              'input_sha256': sha256(raw).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'records': [], 'winners': {}, 'completed_trials': 0, 'attempted_trials': 0,
              'budget_completed': False}
    constructed = []
    for case in source['cases']:
        seed, common = primitive_seed(case['seed'])
        constructed.append((AffinePlaceOracle(seed), common))
    tasks = [(i, j, start) for start in range(args.starts)
             for j in range(max(len(case['partitions']) for case in source['cases']))
             for i, case in enumerate(source['cases']) if j < len(case['partitions'])]
    began = time.monotonic()
    deadline = began + args.seconds
    for case_index, partition_index, restart in tasks:
        if time.monotonic() >= deadline:
            break
        case = source['cases'][case_index]
        oracle, content = constructed[case_index]
        # Inputs retain the balanced t convention. At theta=pi/4 the
        # supporting-line convention divides t and all actual heights by sqrt2.
        baseline = dict(case['partitions'][partition_index])
        baseline['outer_log_scale'] /= math.sqrt(2)
        trial_start = time.monotonic()
        pressure = ScaledMixedSupportingPressure(oracle, baseline, 'full', False, True,
                                                 min(deadline, trial_start + args.trial_seconds))
        values = pressure.initial(baseline, rng, restart)
        values[pressure.angle_index] = (math.pi / 4, .22, 1.35, .48, 1.09)[restart % 5]
        if restart == 0:
            for cut, groups in enumerate(pressure.groups):
                for index, group in enumerate(groups):
                    point = 'infinity' if group['point'] is None else str(group['point'][0])
                    if point in case.get('initial_heights', {}):
                        heights = np.asarray(case['initial_heights'][point], float)
                        values[pressure.height_slices[cut][index]] = heights[pressure.active_heights]
            if 'initial_types' in case:
                values[pressure.type_slice] = np.concatenate([np.log(np.asarray(p)[1:] / p[0])
                                                              for p in case['initial_types']])
        bounds = pressure.bounds()
        for slices in pressure.height_slices:
            for current in slices:
                bounds[current] = [(-48, 48)] * (current.stop - current.start)
        outcomes, completed, error = [], True, None
        try:
            pressure.evaluate(values)
            for smoothing in (.02, .002, .0002):
                optimized = minimize(lambda x: pressure.evaluate(x, smoothing), values, jac=True,
                                     method='L-BFGS-B', bounds=bounds,
                                     options={'maxiter': args.iterations, 'maxfun': 1600,
                                              'ftol': 2e-11, 'gtol': 3e-7, 'maxls': 25})
                outcomes.append({'smoothing': smoothing, 'success': bool(optimized.success),
                                 'iterations': int(optimized.nit), 'evaluations': int(optimized.nfev),
                                 'message': str(optimized.message)})
                values = np.asarray(pressure.best['optimizer_parameters'])
        except SearchDeadline:
            completed, error = False, 'Declared trial or global deadline reached.'
        except (ArithmeticError, ValueError) as exception:
            completed, error = False, f'{type(exception).__name__}: {exception}'
        record = {'case_index': case_index, 'label': case['label'], 'partition_index': partition_index,
                  'balanced_input_partition': case['partitions'][partition_index], 'restart': restart,
                  'degree': oracle.degree, 'common_factor': content, 'completed': completed,
                  'queries': pressure.queries, 'seconds': time.monotonic() - trial_start,
                  'outcomes': outcomes, 'interruption': error,
                  'best': {key: value for key, value in (pressure.best or {}).items() if key != 'place_records'}}
        report['records'].append(record)
        key = str(case_index)
        if pressure.best and (key not in report['winners'] or pressure.best['support_value'] < report['winners'][key]['best']['support_value']):
            report['winners'][key] = record | {'best': pressure.best}
        report.update(attempted_trials=len(report['records']),
                      completed_trials=sum(row['completed'] for row in report['records']),
                      elapsed_seconds=time.monotonic() - began, budget_completed=time.monotonic() >= deadline)
        atomic_write(args.output, report)
        print(json.dumps({'case': case['label'], 'partition': partition_index, 'restart': restart,
                          'support': pressure.best['support_value'] if pressure.best else None,
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(elapsed_seconds=time.monotonic() - began, budget_completed=time.monotonic() >= deadline,
                  trial_list_completed=len(report['records']) == len(tasks))
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
