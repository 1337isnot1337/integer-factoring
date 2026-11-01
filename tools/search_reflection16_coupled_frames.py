#!/usr/bin/env python3
"""Bounded complete-profile search on coupled degree-preserving Laurent frames.

The input catalogue is exact. Every selected four-leg frame is expanded and
its final polynomiality and endpoint are checked again. Numerical pressure
records are candidates only, never exact exponent or optimality claims.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from reflection_laurent_physical_frames import (
    bell_normalize, bell_restore, coupled_transformed, digest,
)
from search_physical_code_gauges import start_params


class Deadline(Exception):
    pass


def alarm_handler(signum, frame):
    raise Deadline('per-trial wall-time budget reached')


def save(path, data):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, separators=(',', ':')) + '\n')
    temporary.replace(path)


def scaled_to_actual(variables):
    result = np.asarray(variables, float).copy()
    result[:-2] *= result[-2]
    return result


def actual_to_scaled(variables):
    result = np.asarray(variables, float).copy()
    result[:-2] /= result[-2]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalogue', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--trial-seconds', type=float, default=70)
    parser.add_argument('--iterations', type=int, default=80)
    parser.add_argument('--seed', type=int, default=312119)
    parser.add_argument('--trials', type=int, default=100000)
    args = parser.parse_args()
    if min(args.seconds, args.trial_seconds, args.iterations, args.trials) <= 0:
        parser.error('All budgets must be positive.')
    signal.signal(signal.SIGALRM, alarm_handler)
    raw = args.catalogue.read_bytes()
    catalogue = json.loads(raw)
    if catalogue['alphabet'] != 16 or catalogue['interval'] != 8:
        parser.error('Use the exact sixteen-letter interval-eight catalogue.')
    base = tensor(8)
    bell, factors = bell_normalize(base)
    assert bell_restore(bell, factors)['values'] == base['values']
    bit_energies = physical_bit_energies(8)
    energy = np.asarray([sum(value for i, value in enumerate(bit_energies) if state >> i & 1)
                         for state in range(16)], float)
    baseline = AsymmetricPressure(base, 1)

    def natural(scale):
        return np.r_[start_params(energy * scale, scale) / math.sqrt(2), math.pi / 4]

    fitted = minimize_scalar(lambda scale: baseline.evaluate(natural(scale), gradient=False)[0],
                             bounds=(.005, .3), method='bounded', options={'xatol': 1e-10})
    initial_actual = natural(float(fitted.x))
    baseline_value, baseline_details = baseline.evaluate(initial_actual, gradient=False)
    initial_scaled = actual_to_scaled(initial_actual)
    ordered = sorted(range(16), key=lambda x: (energy[x], x))
    priority = {symbol: i for i, symbol in enumerate(ordered)}
    choices = []
    for row in catalogue['records']:
        for permitted in row['allowed']:
            choices.append((row['target'], row['source'], permitted['k']))
    choices.sort(key=lambda c: (priority[c[0]] + priority[c[1]] + .5 * c[2],
                                c[2], priority[c[0]], priority[c[1]]))
    rng = random.Random(args.seed)
    shuffled = choices.copy()
    rng.shuffle(shuffled)
    report = {'scope': 'Bounded floating complete-cost search on exactly polynomial coupled frames; no certified support inequality, global optimum, or factoring improvement.',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'arguments': vars(args) | {'catalogue': str(args.catalogue), 'output': str(args.output)},
              'catalogue_sha256': sha256(raw).hexdigest(),
              'alphabet': 16, 'degree': base['degree'], 'baseline_support': baseline_value,
              'baseline_balanced_rho': baseline_value / math.sqrt(2),
              'baseline_params': initial_actual.tolist(), 'baseline_types':
              {'x': baseline_details['type_x'], 'y': baseline_details['type_y']},
              'catalogue_candidates': len(choices), 'records': [], 'best_candidate': None,
              'budget_completed': False, 'completed_trials': 0, 'attempted_trials': 0,
              'counters': {'objective_calls': 0, 'unresolved_calls': 0, 'rejected_poles': 0},
              'height_parameterization': 'Physical heights h=tH; optimized H is bounded, while one positive t scales every local order.'}
    names = (Path(__file__).name, 'reflection_laurent_physical_frames.py',
             'reflection_dual_rail_gaussian.py', 'physical_split_crossings.py',
             'asymmetric_crossing_pressure.py', 'search_physical_code_gauges.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    directory = Path(__file__).parent
    report['source_sha256'] = {name: sha256((directory / name).read_bytes()).hexdigest() for name in names}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    trial = 0
    while time.monotonic() - started < args.seconds and trial < args.trials:
        if trial < 2 * len(choices):
            target, source, pole = choices[trial // 2]
            scalar = 1
        else:
            target, source, pole = shuffled[(trial // 2) % len(shuffled)]
            scalar = 2 if trial % 4 < 2 else 3
        mode = 'balanced' if trial % 2 == 0 else 'asymmetric'
        begin = time.monotonic()
        try:
            changed = coupled_transformed(bell, target, source, -pole, scalar)
            seed = bell_restore(changed, factors)
        except ArithmeticError:
            report['counters']['rejected_poles'] += 1
            report['records'].append({'trial': trial, 'frame': [target, source, -pole, scalar],
                                      'mode': mode, 'status': 'exact final pole rejected'})
            trial += 1
            continue
        assert seed['degree'] <= base['degree']
        oracle = AsymmetricPressure(seed, 1)
        valid = [math.inf, None, None]
        before = dict(report['counters'])

        def objective(variables):
            report['counters']['objective_calls'] += 1
            actual = scaled_to_actual(variables)
            try:
                value, derivative, details = oracle.evaluate(actual)
            except (ArithmeticError, ValueError):
                report['counters']['unresolved_calls'] += 1
                return 1e6, np.zeros(len(variables))
            if value < valid[0]:
                valid[:] = [value, np.asarray(variables, float).copy(), details]
            gradient = derivative.copy()
            gradient[:-2] *= variables[-2]
            gradient[-2] += float(np.asarray(variables[:-2]) @ derivative[:-2])
            return value, gradient

        initial = initial_scaled.copy()
        previous = report['best_candidate']
        if previous is not None and trial % 4 >= 2:
            initial = np.asarray(previous['best_scaled_parameters'], float)
        if mode == 'balanced':
            initial[-1] = math.pi / 4
        bounds = [(-48, 48)] * (len(initial) - 2) + [(.001, .8)]
        bounds.append((math.pi / 4, math.pi / 4) if mode == 'balanced' else (.025, math.pi / 2 - .025))
        completed, optimizer = True, {}
        remaining = args.seconds - (time.monotonic() - started)
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, remaining)))
        try:
            initial_value = objective(initial)[0]
            result = minimize(objective, initial, jac=True, method='SLSQP', bounds=bounds,
                              options={'maxiter': args.iterations, 'ftol': 2e-10})
            optimizer = {'initial_support': initial_value, 'success': bool(result.success),
                         'message': str(result.message), 'iterations': int(result.nit),
                         'function_evaluations': int(result.nfev)}
        except Deadline:
            completed = False
            optimizer = {'message': 'per-trial or remaining global budget interrupted optimization'}
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        record = {'trial': trial, 'frame': [target, source, -pole, scalar], 'mode': mode,
                  'matrix_sha256': digest(seed['values']), 'degree': seed['degree'],
                  'optimizer_completed': completed, 'optimizer': optimizer,
                  'elapsed_seconds': time.monotonic() - begin,
                  'counters': {k: report['counters'][k] - before[k] for k in before}}
        if valid[1] is not None:
            record.update(best_support=valid[0], best_scaled_parameters=valid[1].tolist(),
                          best_parameters=scaled_to_actual(valid[1]).tolist(),
                          type_entropy=valid[2]['entropy_sum'], weights=valid[2]['weights'],
                          type_x=valid[2]['type_x'], type_y=valid[2]['type_y'])
            if report['best_candidate'] is None or valid[0] < report['best_candidate']['best_support']:
                report['best_candidate'] = record
        report['records'].append(record)
        report['attempted_trials'] = trial + 1
        report['completed_trials'] += int(completed)
        report['elapsed_seconds'] = time.monotonic() - started
        save(args.output, report)
        print(json.dumps({'attempted': report['attempted_trials'], 'completed': report['completed_trials'],
                          'best_support': report['best_candidate']['best_support'] if report['best_candidate'] else None,
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        trial += 1
        if valid[0] < 1 - 1e-6:
            report['stopping_reason'] = 'strict numerical supporting-line candidate requires independent exact replay'
            break
    report['elapsed_seconds'] = time.monotonic() - started
    report['budget_completed'] = report['elapsed_seconds'] >= args.seconds
    report['trial_limit_completed'] = trial >= args.trials
    report['stopping_reason'] = report.get('stopping_reason', 'declared time or trial budget completed')
    save(args.output, report)


if __name__ == '__main__':
    main()
