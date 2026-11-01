#!/usr/bin/env python3
"""Bounded complete-cost search on arbitrary-rank q16 coupled frames."""
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
from f16_identity_operator_profiles import inverse_matrix
from physical_split_crossings import F4_PRODUCT, identity_split
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from reflection_laurent_physical_frames import bell_normalize, digest
from search_physical_code_gauges import start_params
from search_reflection16_coupled_frames import (
    Deadline, actual_to_scaled, alarm_handler, save, scaled_to_actual,
)
from square_zero_coupled_reflection_frames import construct, matching_generator, product


def random_generator(index, rng, energy):
    classes = [[i for i in range(16) if i.bit_count() % 2 == parity] for parity in (0, 1)]
    permutation = list(range(16))
    if index < 3:
        mask = (3, 12, 15)[index]
        permutation = [i ^ mask for i in range(16)]
        name = 'parity-preserving-xor-' + str(mask)
    else:
        for group in classes:
            if index == 3:
                group.sort(key=lambda i: (energy[i], i))
            else:
                rng.shuffle(group)
            stop = 4 if index % 5 == 4 else 8
            for i in range(0, stop, 2):
                a, b = group[i:i + 2]
                permutation[a], permutation[b] = b, a
        name = 'energy-adjacent-matching' if index == 3 else 'random-parity-matching'
    generator = matching_generator(permutation)
    if index >= 6 and index % 3:
        change = identity_split(16)
        for _ in range(rng.randrange(4, 17)):
            group = rng.choice(classes)
            target, source = rng.sample(group, 2)
            scalar = rng.choice((1, 2, 3))
            change[target] = [a ^ F4_PRODUCT[scalar][b] for a, b in zip(change[target], change[source])]
        generator = product(product(change, generator), [list(row) for row in inverse_matrix(change)])
        name += '-F4-conjugated'
    return {'name': name, 'generator': generator, 'coordinates': 'bell', 'pole': 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generators', type=Path, help='optional exact queue with cases containing generator, coordinates, pole')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--trial-seconds', type=float, default=90)
    parser.add_argument('--iterations', type=int, default=120)
    parser.add_argument('--trials', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=312136)
    args = parser.parse_args()
    if min(args.seconds, args.trial_seconds, args.iterations, args.trials) <= 0:
        parser.error('All budgets must be positive.')
    signal.signal(signal.SIGALRM, alarm_handler)
    rng = random.Random(args.seed)
    base = tensor(8)
    bell_data = bell_normalize(base)
    energy = np.asarray([sum(value for i, value in enumerate(physical_bit_energies(8)) if state >> i & 1)
                         for state in range(16)], float)
    baseline = AsymmetricPressure(base, 1)

    def natural(scale):
        return np.r_[start_params(energy * scale, scale) / math.sqrt(2), math.pi / 4]

    fit = minimize_scalar(lambda scale: baseline.evaluate(natural(scale), gradient=False)[0],
                          bounds=(.005, .3), method='bounded', options={'xatol': 1e-10})
    natural_actual = natural(float(fit.x))
    natural_scaled = actual_to_scaled(natural_actual)
    queue = None
    if args.generators:
        raw = args.generators.read_bytes()
        queue = json.loads(raw)['cases']
        if not queue:
            parser.error('The supplied generator queue is empty.')
    report = {'scope': 'Bounded floating complete-cost search; every full coefficient tensor is checked exactly. No support, optimum, or factoring improvement is certified.',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'arguments': vars(args) | {'output': str(args.output), 'generators': str(args.generators) if args.generators else None},
              'baseline_balanced_support': float(fit.fun), 'baseline_balanced_rho': float(fit.fun) / math.sqrt(2),
              'baseline_parameters': natural_actual.tolist(), 'alphabet': 16, 'source_degree': base['degree'],
              'queue_sha256': sha256(raw).hexdigest() if args.generators else None,
              'records': [], 'best_by_mode': {}, 'attempted_trials': 0, 'completed_optimizations': 0,
              'completed_budget': False,
              'counters': {'objective_calls': 0, 'unresolved_calls': 0, 'rejected_poles': 0}}
    names = (Path(__file__).name, 'square_zero_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py',
             'physical_split_crossings.py', 'asymmetric_crossing_pressure.py',
             'search_reflection16_coupled_frames.py', 'search_physical_code_gauges.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    begun = time.monotonic()
    trial, candidate = 0, None
    while trial < args.trials and time.monotonic() - begun < args.seconds:
        if trial % 2 == 0:
            index = trial // 2
            candidate = dict(queue[index % len(queue)]) if queue else random_generator(index, rng, energy)
        mode = 'balanced' if trial % 2 == 0 else 'asymmetric'
        started = time.monotonic()
        try:
            seed = construct(base, candidate['generator'], candidate.get('coordinates', 'bell'),
                             pole=candidate.get('pole', 1), bell_data=bell_data)
        except ArithmeticError:
            report['counters']['rejected_poles'] += 1
            report['records'].append({'trial': trial, 'candidate': candidate, 'mode': mode,
                                      'status': 'exact complete pole test rejected'})
            trial += 1
            continue
        oracle = AsymmetricPressure(seed, 1)
        best = [math.inf, None, None]
        counts = dict(report['counters'])

        def objective(variables):
            report['counters']['objective_calls'] += 1
            try:
                value, derivative, details = oracle.evaluate(scaled_to_actual(variables))
            except (ArithmeticError, ValueError):
                report['counters']['unresolved_calls'] += 1
                return 1e6, np.zeros(len(variables))
            if value < best[0]:
                best[:] = [value, np.asarray(variables, float).copy(), details]
            gradient = derivative.copy()
            gradient[:-2] *= variables[-2]
            gradient[-2] += float(np.asarray(variables[:-2]) @ derivative[:-2])
            return value, gradient

        initial = natural_scaled.copy()
        if trial >= 4 and trial % 4 >= 2 and mode in report['best_by_mode']:
            initial = np.asarray(report['best_by_mode'][mode]['best_scaled_parameters'], float)
        if mode == 'balanced':
            initial[-1] = math.pi / 4
        bounds = [(-48, 48)] * (len(initial) - 2) + [(.001, .8)]
        bounds.append((math.pi / 4, math.pi / 4) if mode == 'balanced' else (.025, math.pi / 2 - .025))
        completed, optimizer = True, {}
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, args.seconds - (time.monotonic() - begun))))
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
        record = {'trial': trial, 'candidate': candidate, 'mode': mode,
                  'rank': seed['square_zero_frame']['rank'], 'degree': seed['degree'],
                  'matrix_sha256': digest(seed['values']), 'optimizer_completed': completed,
                  'optimizer': optimizer, 'seconds': time.monotonic() - started,
                  'counters': {key: report['counters'][key] - counts[key] for key in counts}}
        if best[1] is not None:
            record.update(best_support=best[0], best_parameters=scaled_to_actual(best[1]).tolist(),
                          best_scaled_parameters=best[1].tolist(), weights=best[2]['weights'],
                          type_entropy=best[2]['entropy_sum'], type_x=best[2]['type_x'], type_y=best[2]['type_y'])
            if mode not in report['best_by_mode'] or best[0] < report['best_by_mode'][mode]['best_support']:
                report['best_by_mode'][mode] = record
        report['records'].append(record)
        report['attempted_trials'] = trial + 1
        report['completed_optimizations'] += int(completed)
        report['elapsed_seconds'] = time.monotonic() - begun
        save(args.output, report)
        print(json.dumps({'attempted': report['attempted_trials'], 'completed': report['completed_optimizations'],
                          'best': {key: row['best_support'] for key, row in report['best_by_mode'].items()},
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        trial += 1
        if best[0] < 1 - 1e-6:
            report['stopping_reason'] = 'strict numerical supporting-line candidate requires independent exact replay'
            break
    report['elapsed_seconds'] = time.monotonic() - begun
    report['completed_budget'] = report['elapsed_seconds'] >= args.seconds
    report['completed_trial_limit'] = trial >= args.trials
    report['stopping_reason'] = report.get('stopping_reason', 'declared time or trial budget completed')
    save(args.output, report)


if __name__ == '__main__':
    main()
