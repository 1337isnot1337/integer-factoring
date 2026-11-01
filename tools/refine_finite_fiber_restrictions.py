#!/usr/bin/env python3
"""Bounded supporting-line refinement of explicit restricted crossings."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_laurent_physical_frames import digest
from search_affine_coupled_reflection_frames import balanced
from search_fused_mixed_weighted import SearchDeadline, atomic_write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trial-seconds', type=float, default=6)
    parser.add_argument('--iterations', type=int, default=300)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    cases = json.loads(raw)['cases']
    names = (Path(__file__).name, 'asymmetric_crossing_pressure.py',
             'search_affine_coupled_reflection_frames.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Complete native weighted supporting-line search on explicit selected restrictions; no exact threshold certificate.',
              'input_sha256': sha256(raw).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'arguments': {'trial_seconds': args.trial_seconds, 'iterations': args.iterations},
              'records': [], 'winners': {}, 'completed': False}
    began = time.monotonic()
    for case in cases:
        seed = case['seed']
        assert digest(seed['values']) == case['tensor_sha256']
        oracle = AsymmetricPressure(seed, 1)
        warm = np.asarray(case['initial_params'])
        assert abs(balanced(oracle, warm)[0] - case['initial_balanced_rho']) < 2e-9
        for angle in (.25, math.pi / 4, 1.32):
            started = time.monotonic()
            deadline = started + args.trial_seconds
            best = None
            queries = unresolved = 0
            completed = True
            outcome = None

            def objective(values):
                nonlocal best, queries, unresolved
                if time.monotonic() >= deadline:
                    raise SearchDeadline()
                queries += 1
                try:
                    value, gradient, details = oracle.evaluate(values)
                    if best is None or value < best['support_value']:
                        best = {'support_value': float(value),
                                'parameters': list(map(float, values)), 'details': details}
                    return value, gradient
                except (ArithmeticError, ValueError):
                    unresolved += 1
                    return 1e6, np.zeros(len(values))

            try:
                start = np.r_[warm / math.sqrt(2), angle]
                objective(start)
                fit = minimize(objective, start, jac=True, method='SLSQP',
                               bounds=[(-32, 32)] * (len(warm) - 1)
                                      + [(.0005, 2), (.015, math.pi / 2 - .015)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                outcome = {'success': bool(fit.success), 'iterations': int(fit.nit),
                           'message': str(fit.message)}
            except SearchDeadline:
                completed = False
            record = {'label': case['label'], 'tensor_sha256': case['tensor_sha256'],
                      'initial_angle': angle, 'best': best, 'completed': completed,
                      'outcome': outcome, 'queries': queries, 'unresolved': unresolved,
                      'seconds': time.monotonic() - started}
            report['records'].append(record)
            if best and (case['label'] not in report['winners']
                         or best['support_value'] < report['winners'][case['label']]['best']['support_value']):
                report['winners'][case['label']] = record
            report['elapsed_seconds'] = time.monotonic() - began
            atomic_write(args.output, report)
            print(json.dumps({'label': case['label'], 'initial_angle': angle,
                              'value': None if best is None else best['support_value'],
                              'seconds': report['elapsed_seconds']}), flush=True)
            if best and best['support_value'] < 1 - 1e-7:
                report['stopping_reason'] = 'Numerical threshold candidate requires exact replay.'
                report['completed'] = True
                atomic_write(args.output, report)
                return
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  stopping_reason='Declared finite trial list completed.')
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
