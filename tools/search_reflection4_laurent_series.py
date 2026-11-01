#!/usr/bin/env python3
"""Exhaust the q4 two-coefficient upper-block Laurent-series family.

The 65,536 words encode eight F4 coefficients, including the zero frame.
Every reconstructed tensor is checked before complete natural-gauge costs
are evaluated.  Retained words receive bounded full gauge/type refinement.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from laurent_series_coupled_frames import construct
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, digest
from search_affine_coupled_reflection_frames import balanced
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import Deadline, alarm
from search_fused_mixed_weighted import atomic_write as save


def decode(index):
    if not 0 <= index < 4 ** 8:
        raise ValueError('The q4 series word is an integer from0 through65,535.')
    result = {}
    for pole in (1, 2):
        matrix = [[0] * 4 for _ in range(4)]
        for i in range(2):
            for j in range(2, 4):
                matrix[i][j], index = index % 4, index // 4
        if any(value for row in matrix for value in row):
            result[pole] = matrix
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--screens', type=int, default=65536)
    parser.add_argument('--refine', type=int, default=32)
    parser.add_argument('--trial-seconds', type=float, default=20)
    parser.add_argument('--iterations', type=int, default=250)
    args = parser.parse_args()
    if not 1 <= args.screens <= 65536 or min(args.seconds, args.refine, args.trial_seconds, args.iterations) <= 0:
        parser.error('Use positive budgets and at most65,536 screens.')
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(4)
    bell_data = bell_normalize(base)
    original = AsymmetricPressure(base, 1)

    def natural(scale):
        return start_params(np.asarray([0, 12, 4, 16]) * scale, scale)

    fit = minimize_scalar(lambda s: balanced(original, natural(s))[0], bounds=(.005, 1),
                          method='bounded', options={'xatol': 1e-10})
    initial = natural(float(fit.x))
    names = (Path(__file__).name, 'laurent_series_coupled_frames.py', 'square_zero_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py',
             'asymmetric_crossing_pressure.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'search_affine_coupled_reflection_frames.py')
    report = {'scope': 'Declared finite exact family with floating complete-cost screening and bounded refinement; no optimum or factoring certificate.',
              'arguments': vars(args) | {'output': str(args.output)}, 'alphabet': 4, 'degree': 24,
              'coefficient_slots': [[pole, i, j] for pole in (1, 2) for i in range(2) for j in range(2, 4)],
              'encoding': 'Least significant base-four digit corresponds to the first coefficient slot.',
              'baseline_rho': float(fit.fun), 'baseline_scale': float(fit.x),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'screen_records': [], 'retained': [], 'refinements': [], 'unresolved_calls': 0,
              'completed_budget': False, 'completed_full_family': False}
    started = time.monotonic()
    for index in range(args.screens):
        if time.monotonic() - started >= .8 * args.seconds:
            break
        seed = base if index == 0 else construct(base, decode(index), bell_data=bell_data)
        oracle = AsymmetricPressure(seed, 1)
        try:
            rho = min(balanced(oracle, initial * factor)[0] for factor in (1, .5))
        except (ArithmeticError, ValueError):
            rho = None
            report['unresolved_calls'] += 1
        report['screen_records'].append([index, rho])
        if rho is not None:
            retained = report['retained']
            if len(retained) < args.refine or rho < retained[-1]['rho']:
                retained.append({'index': index, 'rho': rho, 'tensor_sha256': digest(seed['values'])})
                retained.sort(key=lambda row: row['rho'])
                del retained[args.refine:]
        if len(report['screen_records']) % 256 == 0:
            report.update(elapsed_seconds=time.monotonic() - started,
                          completed_screens=len(report['screen_records']))
            save(args.output, report)
            print(json.dumps({'screens': len(report['screen_records']), 'seconds': report['elapsed_seconds'],
                              'best_natural_rho': report['retained'][0]['rho']}), flush=True)
    report['completed_screens'] = len(report['screen_records'])
    report['completed_full_family'] = len(report['screen_records']) == 65536
    warm = initial.copy()
    for selected in report['retained']:
        if time.monotonic() - started >= args.seconds:
            break
        index = selected['index']
        seed = base if index == 0 else construct(base, decode(index), bell_data=bell_data)
        assert digest(seed['values']) == selected['tensor_sha256']
        oracle = AsymmetricPressure(seed, 1)
        best = [math.inf, None, None]

        def objective(parameters):
            try:
                value, gradient, details = balanced(oracle, parameters)
            except (ArithmeticError, ValueError):
                report['unresolved_calls'] += 1
                return 1e6, np.zeros(len(parameters))
            if value < best[0]:
                best[:] = [value, np.asarray(parameters).copy(), details]
            return value, gradient

        begin = time.monotonic()
        outcomes, complete = [], True
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, args.seconds - (begin - started))))
        try:
            for point in (warm.copy(), initial.copy()):
                fitted = minimize(objective, point, jac=True, method='SLSQP',
                                  bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                                  options={'maxiter': args.iterations, 'ftol': 2e-11})
                outcomes.append({'success': bool(fitted.success), 'iterations': int(fitted.nit),
                                 'evaluations': int(fitted.nfev), 'message': str(fitted.message)})
        except Deadline:
            complete = False
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        row = dict(selected, optimizer_completed=complete, outcomes=outcomes,
                   seconds=time.monotonic() - begin)
        if best[1] is not None:
            warm = best[1].copy()
            row.update(best_rho=best[0], best_parameters=best[1].tolist(), details=best[2])
        report['refinements'].append(row)
        report['elapsed_seconds'] = time.monotonic() - started
        save(args.output, report)
        print(json.dumps({'refined_index': index, 'rho': best[0], 'seconds': report['elapsed_seconds']}), flush=True)
    report.update(elapsed_seconds=time.monotonic() - started,
                  completed_budget=time.monotonic() - started >= args.seconds,
                  completed_refinements=len(report['refinements']))
    save(args.output, report)


if __name__ == '__main__':
    main()
