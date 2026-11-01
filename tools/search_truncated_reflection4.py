#!/usr/bin/env python3
"""Full-gauge, common-type, and asymmetric search of all even q4 truncations.

Every candidate is a complete strict crossing constructed by the exact
truncation module. Numerical minima are exploratory, never certificates.
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
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import Deadline, alarm, save
from truncate_reflection_crossings import truncate_even_crossing


def balanced(oracle, params):
    value, gradient, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
    return value / math.sqrt(2), gradient[:-1] / 2, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--trial-seconds', type=float, default=55)
    parser.add_argument('--iterations', type=int, default=350)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(4)
    source_names = (Path(__file__).name, 'truncate_reflection_crossings.py',
                    'asymmetric_crossing_pressure.py', 'reflection_dual_rail_gaussian.py',
                    'search_physical_code_gauges.py', 'native_rectangular_smith_f4.py',
                    'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'source_sha256': {
        name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in source_names},
        'requested_seconds': args.seconds, 'trial_seconds': args.trial_seconds,
        'iteration_limit': args.iterations, 'planned_cases': 24, 'cases': [],
        'completed': False, 'best_rho': None, 'best_support_value': None}
    began = time.monotonic()
    best_rho = best_support = math.inf
    threshold = False
    for bound in (24, 16, 12, 14, 18, 20, 22, 10, 8, 6, 4, 2):
        for cut in (1, 0):
            if time.monotonic() - began >= args.seconds or threshold:
                break
            seed = truncate_even_crossing(base, bound)
            oracle = AsymmetricPressure(seed, cut)
            record = {'degree_bound': bound, 'actual_degree': seed['degree'], 'cut_at_infinity': cut,
                      'tensor_sha256': digest(seed['values']), 'screens': [], 'optimizations': [],
                      'calls': 0, 'unresolved': 0}
            valid = {'balanced': [math.inf, None, None], 'asymmetric': [math.inf, None, None]}
            phase = 'balanced'
            case_began = time.monotonic()

            def objective(params):
                record['calls'] += 1
                try:
                    value, gradient, details = balanced(oracle, params) if phase == 'balanced' else oracle.evaluate(params)
                    if value < valid[phase][0]:
                        valid[phase][:] = [value, np.asarray(params).copy(), details]
                    return value, gradient
                except (ArithmeticError, ValueError):
                    record['unresolved'] += 1
                    return 1e6, np.zeros(len(params))

            signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds,
                                max(.01, args.seconds - (case_began - began))))
            try:
                starts = []
                for multiplier in (1, .75, .5, .25, 0):
                    natural = lambda scale: start_params(np.asarray([0, 12, 4, 16]) * scale * multiplier, scale)
                    fit = minimize_scalar(lambda scale: objective(natural(scale))[0], bounds=(.002, 1.5),
                                          method='bounded', options={'xatol': 1e-9})
                    starts.append((float(fit.fun), natural(fit.x)))
                    record['screens'].append({'height_multiplier': multiplier,
                                              'scale': float(fit.x), 'rho': float(fit.fun)})
                for _, initial in sorted(starts, key=lambda row: row[0])[:2]:
                    fit = minimize(objective, initial, jac=True, method='SLSQP',
                                   bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                                   options={'maxiter': args.iterations, 'ftol': 2e-11})
                    record['optimizations'].append({'phase': phase, 'success': bool(fit.success),
                                                    'iterations': int(fit.nit), 'message': str(fit.message)})
                initial = valid['balanced'][1]
                if initial is None:
                    initial = start_params(np.zeros(4), .05)
                phase = 'asymmetric'
                for angle in (math.pi / 4, .25, math.pi / 2 - .25):
                    fit = minimize(objective, np.r_[initial / math.sqrt(2), angle], jac=True, method='SLSQP',
                                   bounds=[(-24, 24)] * 24 + [(.0005, 2), (.015, math.pi / 2 - .015)],
                                   options={'maxiter': args.iterations, 'ftol': 2e-11})
                    record['optimizations'].append({'phase': phase, 'success': bool(fit.success),
                                                    'iterations': int(fit.nit), 'message': str(fit.message)})
                record['status'] = 'all scalar, full-gauge, and asymmetric loops completed'
            except Deadline:
                record['status'] = 'trial deadline during ' + phase + '; resolved calls retained'
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            for mode, (value, params, details) in valid.items():
                if params is not None:
                    record[mode] = {'value': value, 'params': params.tolist(), 'details': details}
            best_rho = min(best_rho, valid['balanced'][0])
            best_support = min(best_support, valid['asymmetric'][0])
            record['elapsed_seconds'] = time.monotonic() - case_began
            report['cases'].append(record)
            report.update(completed_cases=len(report['cases']), elapsed_seconds=time.monotonic() - began,
                          best_rho=None if math.isinf(best_rho) else best_rho,
                          best_support_value=None if math.isinf(best_support) else best_support)
            save(args.output, report)
            print(json.dumps({k: report[k] for k in ('completed_cases', 'elapsed_seconds',
                                                    'best_rho', 'best_support_value')}), flush=True)
            threshold = best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  queue_exhausted=len(report['cases']) == report['planned_cases'],
                  stopping_reason='numerical candidate needs independent exact replay' if threshold else
                  ('all declared candidates completed' if len(report['cases']) == report['planned_cases'] else
                   'requested wall-clock budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
