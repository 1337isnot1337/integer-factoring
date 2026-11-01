#!/usr/bin/env python3
"""Exact paired subalphabet restrictions of the two primitive fiber cores.

All two/three-letter subsets are retained, including unequal paired
alphabet sizes. Endpoint/content checks are exact. Complete native
weighted pressures are exploratory and are never used as certificates.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from affine_place_profile_oracle import primitive_seed
from asymmetric_crossing_pressure import AsymmetricPressure
from physical_split_crossings import retract
from reflection_laurent_physical_frames import digest
from search_affine_coupled_reflection_frames import balanced
from search_fused_mixed_weighted import SearchDeadline, atomic_write


def restriction(seed, first, second):
    def inclusion(subset):
        return [[int(index == symbol) for index in range(4)] for symbol in subset]
    left, right = inclusion(first), inclusion(second)
    restricted = retract(seed, left, left, right, right)
    result, content = primitive_seed(restricted)
    return result, content


def restricted_parameters(params, first, second):
    params = list(params)
    result = []
    for side in range(2):
        for leg, subset in enumerate((first, second, second, first)):
            begin = side * 12 + leg * 3
            heights = [0] + params[begin:begin + 3]
            result.extend(heights[symbol] - heights[subset[0]] for symbol in subset[1:])
    result.append(params[-1])
    return np.asarray(result, dtype=float)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('evidence/deflated-finite-fiber-input.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--trial-seconds', type=float, default=2.5)
    parser.add_argument('--iterations', type=int, default=180)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    inputs = json.loads(raw)
    subsets = [list(value) for size in (2, 3) for value in combinations(range(4), size)]
    # Begin with the hot support and its nearest three-letter extensions.
    subsets.sort(key=lambda s: (not {0, 2}.issubset(s), len(s), s))
    names = (Path(__file__).name, 'affine_place_profile_oracle.py',
             'asymmetric_crossing_pressure.py', 'physical_split_crossings.py',
             'reflection_laurent_physical_frames.py',
             'search_affine_coupled_reflection_frames.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py')
    began = time.monotonic()
    deadline = began + args.seconds
    report = {
        'scope': __doc__, 'input_path': str(args.input),
        'input_sha256': sha256(raw).hexdigest(),
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'arguments': {'seconds': args.seconds, 'trial_seconds': args.trial_seconds,
                      'iterations': args.iterations},
        'records': [], 'winners': {}, 'completed': False,
        'counters': {'objective_calls': 0, 'unresolved_calls': 0},
    }
    for first in subsets:
        for second in subsets:
            for source_index, case in enumerate(inputs['cases']):
                if time.monotonic() >= deadline:
                    break
                assert digest(case['seed']['values']) == case['tensor_sha256']
                seed, content = restriction(case['seed'], first, second)
                oracle = AsymmetricPressure(seed, 1)
                warm = restricted_parameters(case['initial_params'], first, second)
                best = None
                trial_started = time.monotonic()
                trial_deadline = min(deadline, trial_started + args.trial_seconds)
                outcomes = []
                completed = True

                def objective(values):
                    nonlocal best
                    if time.monotonic() >= trial_deadline:
                        raise SearchDeadline()
                    report['counters']['objective_calls'] += 1
                    try:
                        value, gradient, details = balanced(oracle, values)
                        if best is None or value < best['value']:
                            best = {'value': float(value),
                                    'parameters': list(map(float, values)),
                                    'details': details}
                        return value, gradient
                    except (ArithmeticError, ValueError):
                        report['counters']['unresolved_calls'] += 1
                        return 1e6, np.zeros(len(values))

                try:
                    for scale in (1.0, .6):
                        start = warm * scale
                        objective(start)
                        fit = minimize(objective, start, jac=True, method='SLSQP',
                                       bounds=[(-32, 32)] * (len(warm) - 1) + [(.0005, 2)],
                                       options={'maxiter': args.iterations, 'ftol': 2e-11})
                        outcomes.append({'success': bool(fit.success),
                                         'iterations': int(fit.nit),
                                         'message': str(fit.message)})
                except SearchDeadline:
                    completed = False
                record = {
                    'source_index': source_index, 'source_label': case['label'],
                    'subset_x': first, 'subset_y': second,
                    'seed': seed, 'tensor_sha256': digest(seed['values']),
                    'scalar_content': content, 'degree': seed['degree'],
                    'strict_endpoint_entries': len(seed['values']),
                    'completed': completed, 'outcomes': outcomes,
                    'seconds': time.monotonic() - trial_started, 'best': best,
                }
                report['records'].append(record)
                if best is not None:
                    key = f'{source_index}:{len(first)}x{len(second)}'
                    previous = report['winners'].get(key)
                    if previous is None or best['value'] < previous['best']['value']:
                        report['winners'][key] = record
                report['elapsed_seconds'] = time.monotonic() - began
                report['completed_trials'] = sum(row['completed'] for row in report['records'])
                atomic_write(args.output, report)
                print(json.dumps({'source': source_index, 'x': first, 'y': second,
                                  'degree': seed['degree'],
                                  'rho': None if best is None else best['value'],
                                  'elapsed': report['elapsed_seconds']}), flush=True)
                if best is not None and best['value'] < 1 / math.sqrt(2) - 1e-7:
                    report['stopping_reason'] = 'Numerical threshold candidate needs exact replay.'
                    deadline = 0
            if time.monotonic() >= deadline:
                break
        if time.monotonic() >= deadline:
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  stopping_reason=report.get('stopping_reason', 'Declared case or time budget completed.'),
                  expected_restrictions=2 * len(subsets) ** 2)
    atomic_write(args.output, report)
    print(json.dumps({'completed': True, 'restrictions': len(report['records']),
                      'best': min((r['best']['value'] for r in report['records']
                                   if r['best'] is not None), default=None),
                      'elapsed': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
