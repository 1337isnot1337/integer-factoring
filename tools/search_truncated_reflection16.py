#!/usr/bin/env python3
"""Complete-profile screen and bounded refinement of q16 degree truncations."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import logsumexp, softmax

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from reflection_laurent_physical_frames import digest
from search_fused_mixed_weighted import atomic_write
from search_physical_code_gauges import start_params
from search_reflection16_coupled_frames import Deadline, alarm_handler, scaled_to_actual
from truncate_reflection_crossings import truncate_even_crossing


def support(orders, energy, scale, angle, gradient=False):
    weights = (math.cos(angle), math.sin(angle))
    scaled_orders = scale * orders
    logs = [float(logsumexp(-scaled_orders[side] / weights[side])) for side in range(2)]
    probabilities = [softmax(-scaled_orders[side] / weights[side]) for side in range(2)]
    raw = sum(weight * value for weight, value in zip(weights, logs))
    physical = scale * energy

    def equation(rate):
        return raw - 2 * rate * logsumexp(-physical / rate)

    if equation(1e-10) <= 0:
        raise ArithmeticError('No resolved positive type root.')
    upper = 4.0
    while equation(upper) > 0:
        upper *= 2
        if upper > 1e6:
            raise ArithmeticError('Type root not bracketed.')
    rate = brentq(equation, 1e-10, upper, xtol=2e-12)
    types = softmax(-physical / rate)
    entropy = -2 * float(types @ np.log(types))
    if not gradient:
        return float(rate), types.tolist()
    ds = (-sum(float(p @ row) for p, row in zip(probabilities, orders))
          + 2 * float(types @ energy)) / entropy
    channel_entropies = [value + float(p @ row) / weight
                         for value, p, row, weight in zip(logs, probabilities, scaled_orders, weights)]
    da = (-weights[1] * channel_entropies[0] + weights[0] * channel_entropies[1]) / entropy
    return float(rate), np.asarray([ds, da])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degrees', type=int, nargs='+', default=[96, 64, 80, 48, 88, 72, 32, 68, 92, 16, 8, 2])
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--trial-seconds', type=float, default=60)
    parser.add_argument('--refine', type=int, default=6)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm_handler)
    base = tensor(8)
    energies = physical_bit_energies(8)
    energy = np.asarray([sum(value for bit, value in enumerate(energies) if letter >> bit & 1)
                         for letter in range(16)], float)
    source_names = (Path(__file__).name, 'truncate_reflection_crossings.py',
                    'asymmetric_crossing_pressure.py', 'reflection_dual_rail_gaussian.py',
                    'search_reflection16_coupled_frames.py', 'search_physical_code_gauges.py',
                    'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Exact complete source tensors and fixed rational gauges; numerical type/scale/support refinement only, no factoring certificate.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in source_names},
              'screens': [], 'refinements': [], 'completed_budget': False}
    began, best_cases = time.monotonic(), []
    deadline = began + args.seconds
    for bound in args.degrees:
        if time.monotonic() >= deadline:
            break
        seed = truncate_even_crossing(base, bound)
        oracle = AsymmetricPressure(seed, cut=1)
        for multiplier in (1.0, .75, .5, .25, 0.0):
            if time.monotonic() >= deadline:
                break
            unit = start_params(energy * multiplier, 1.0)
            native = oracle.oracle.prefixes(unit, 1)
            orders = np.asarray(native['orders']).reshape(2, 256)
            assert np.max(np.abs(orders * 4 - np.rint(orders * 4))) < 1e-10
            candidates = []
            for angle in (math.pi / 4, .22, 1.35):
                balanced = minimize_scalar(lambda value: support(orders, energy * multiplier, value, angle)[0],
                                           bounds=(.001, .6), method='bounded', options={'xatol': 2e-10})
                initial = np.asarray([balanced.x, angle])
                optimized = minimize(lambda values: support(orders, energy * multiplier, *values, gradient=True),
                                     initial, jac=True, method='SLSQP', bounds=[(.001, .8), (.015, math.pi / 2 - .015)],
                                     options={'maxiter': 180, 'ftol': 2e-11})
                values = optimized.x if optimized.fun < balanced.fun else initial
                value, types = support(orders, energy * multiplier, *values)
                candidates.append({'support_value': value, 'scale': float(values[0]), 'angle': float(values[1]),
                                   'common_type': types, 'optimizer_success': bool(optimized.success),
                                   'fixed_angle': angle, 'fixed_angle_support': float(balanced.fun),
                                   'fixed_angle_scale': float(balanced.x)})
            selected = min(candidates, key=lambda row: row['support_value'])
            actual = np.r_[unit * selected['scale'], selected['angle']]
            direct, _ = oracle.evaluate(actual, gradient=False)
            assert abs(direct - selected['support_value']) < 2e-9
            record = {'requested_degree': bound, 'actual_degree': seed['degree'],
                      'tensor_sha256': digest(seed['values']), 'height_multiplier': multiplier,
                      'complete_unit_profiles': orders.tolist(), 'natural_scalar_searches': candidates,
                      'best': selected, 'native_objective_control_error': abs(direct - selected['support_value']),
                      'scaled_parameters': np.r_[unit[:-1], selected['scale'], selected['angle']].tolist()}
            report['screens'].append(record)
            best_cases.append(record)
            report.update(elapsed_seconds=time.monotonic() - began, completed_screens=len(report['screens']))
            atomic_write(args.output, report)
            print(json.dumps({'degree': bound, 'multiplier': multiplier, 'support': selected['support_value'],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    selected_cases = []
    for row in sorted(best_cases, key=lambda current: current['best']['support_value']):
        if row['requested_degree'] not in [case['requested_degree'] for case in selected_cases]:
            selected_cases.append(row)
        if len(selected_cases) >= args.refine:
            break
    for start in selected_cases:
        if time.monotonic() >= deadline:
            break
        seed = truncate_even_crossing(base, start['requested_degree'])
        oracle = AsymmetricPressure(seed, 1)
        initial = np.asarray(start['scaled_parameters'])
        best = [math.inf, None, None]
        calls, unresolved = 0, 0

        def objective(values):
            nonlocal calls, unresolved
            calls += 1
            try:
                value, derivative, details = oracle.evaluate(scaled_to_actual(values))
            except (ArithmeticError, ValueError):
                unresolved += 1
                return 1e6, np.zeros(len(values))
            if value < best[0]:
                best[:] = [value, np.asarray(values).copy(), details]
            result = derivative.copy()
            result[:-2] *= values[-2]
            result[-2] += float(np.asarray(values[:-2]) @ derivative[:-2])
            return value, result

        trial_start = time.monotonic()
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, deadline - trial_start)))
        completed, outcome = True, {}
        try:
            initial_value = objective(initial)[0]
            fit = minimize(objective, initial, jac=True, method='SLSQP',
                           bounds=[(-64, 64)] * (len(initial) - 2) + [(.001, .8), (.015, math.pi / 2 - .015)],
                           options={'maxiter': 200, 'ftol': 2e-10})
            outcome = {'initial_support': initial_value, 'success': bool(fit.success),
                       'message': str(fit.message), 'iterations': int(fit.nit), 'evaluations': int(fit.nfev)}
        except Deadline:
            completed = False
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        report['refinements'].append({'requested_degree': start['requested_degree'],
                                      'completed': completed, 'calls': calls, 'unresolved_calls': unresolved,
                                      'seconds': time.monotonic() - trial_start, 'outcome': outcome,
                                      'best_support': best[0] if best[1] is not None else None,
                                      'best_scaled_parameters': best[1].tolist() if best[1] is not None else None,
                                      'best_details': best[2]})
        report.update(elapsed_seconds=time.monotonic() - began)
        atomic_write(args.output, report)
    report.update(elapsed_seconds=time.monotonic() - began, completed_budget=time.monotonic() >= deadline,
                  completed_screen_list=len(report['screens']) == 5 * len(args.degrees),
                  completed_refinement_list=len(report['refinements']) == len(selected_cases))
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
