#!/usr/bin/env python3
"""Complete independent F3 charge masks, then bounded full-cost refinements.

The mask f(x)-f(v)+g(y)-g(u)=0 is the exact sum of three paired diagonal
physical frames over F4. Both label functions vary independently. All
729 normalized pairs are screened; identical tensors are deduplicated.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, degree, scale, valuation
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import Deadline, alarm, save


def masked(base, first, second):
    values = []
    for index, value in enumerate(base['values']):
        x, y, u, v = ((index // 4 ** j) % 4 for j in range(4))
        charge = (first[x] - first[v] + second[y] - second[u]) % 3
        values.append(value if charge == 0 else (0, 0))
        # An independent literal three-character sum includes all coefficients.
        total = (0, 0)
        for j in range(3):
            total = add(total, scale(value, (1, 2, 3)[j * charge % 3]))
        assert total == values[-1]
        if x == v and y == u:
            assert charge == 0
        assert ((values[-1][0].bit_count() & 1), (values[-1][1].bit_count() & 1)) == (int(x == v and y == u), 0)
    actual = max(map(degree, values))
    # The invariant paired trace already forces this degree lower bound.
    assert actual >= 16
    return dict(base, values=values, degree=actual,
                monomial_content=min(valuation(v) for v in values if v != (0, 0)))


def balanced(oracle, params):
    value, gradient, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
    return value / math.sqrt(2), gradient[:-1] / 2, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--refine', type=int, default=80)
    parser.add_argument('--trial-seconds', type=float, default=18)
    parser.add_argument('--iterations', type=int, default=260)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(4)
    initial_oracle = AsymmetricPressure(base, 1)
    natural = lambda s: start_params(np.asarray([0, 12, 4, 16]) * s, s)
    optimal_scale = minimize_scalar(lambda s: balanced(initial_oracle, natural(s))[0], bounds=(.005, 1),
                                    method='bounded', options={'xatol': 1e-10}).x
    initial = natural(optimal_scale)
    baseline = balanced(initial_oracle, initial)[0]
    names = (Path(__file__).name, 'asymmetric_crossing_pressure.py', 'reflection_dual_rail_gaussian.py',
             'reflection_laurent_physical_frames.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'search_physical_code_gauges.py')
    report = {'schema_version': 1, 'scope': __doc__, 'requested_seconds': args.seconds,
              'requested_refinements': args.refine, 'trial_seconds': args.trial_seconds,
              'iteration_limit': args.iterations, 'baseline_rho': baseline,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'screen_records': [], 'cases': [], 'counters': {'calls': 0, 'unresolved': 0}}
    seen, degrees = set(), Counter()
    began = time.monotonic()
    for index, tail in enumerate(product(range(3), repeat=6)):
        first, second = (0,) + tail[:3], (0,) + tail[3:]
        changed = masked(base, first, second)
        digest = sha256(json.dumps(changed['values'], separators=(',', ':')).encode()).hexdigest()
        degrees[changed['degree']] += 1
        if digest in seen:
            continue
        seen.add(digest)
        oracle = AsymmetricPressure(changed, 1)
        values = []
        for multiplier in (1, .5, .2):
            report['counters']['calls'] += 1
            try:
                values.append(balanced(oracle, initial * multiplier)[0])
            except (ArithmeticError, ValueError):
                report['counters']['unresolved'] += 1
                values.append(math.inf)
        report['screen_records'].append({'proposal': index, 'first_charge': first, 'second_charge': second,
                                         'degree': changed['degree'], 'seed_sha256': digest,
                                         'natural_rho': min(values) if min(values) < math.inf else None})
    report.update(proposals_screened=729, unique_tensors=len(seen), screen_complete=True,
                  degree_histogram=dict(degrees), independent_character_entries=729 * 256)
    selected = sorted((r for r in report['screen_records'] if r['natural_rho'] is not None),
                      key=lambda r: r['natural_rho'])[:args.refine]
    best_rho, best_support = baseline, math.inf
    warm = initial.copy()
    threshold = False
    for candidate in selected:
        if time.monotonic() - began >= args.seconds or threshold:
            break
        changed = masked(base, candidate['first_charge'], candidate['second_charge'])
        oracle = AsymmetricPressure(changed, 1)
        valid = {'balanced': [math.inf, None, None], 'asymmetric': [math.inf, None, None]}
        phase = 'balanced'
        record = dict(candidate, optimizations=[])
        start_time = time.monotonic()

        def objective(params):
            report['counters']['calls'] += 1
            try:
                value, gradient, details = balanced(oracle, params) if phase == 'balanced' else oracle.evaluate(params)
                if value < valid[phase][0]:
                    valid[phase][:] = [value, np.asarray(params).copy(), details]
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved'] += 1
                return 1e6, np.zeros(len(params))

        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, args.seconds - (start_time - began))))
        try:
            for start in (warm.copy(), initial * .25):
                result = minimize(objective, start, jac=True, method='SLSQP',
                                  bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                                  options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'success': bool(result.success),
                                                'iterations': int(result.nit), 'message': str(result.message)})
            phase = 'asymmetric'
            start = valid['balanced'][1] if valid['balanced'][1] is not None else initial
            for angle in (.35, 1.22):
                result = minimize(objective, np.r_[start / math.sqrt(2), angle], jac=True, method='SLSQP',
                                  bounds=[(-24, 24)] * 24 + [(.0005, 2), (.025, math.pi / 2 - .025)],
                                  options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'success': bool(result.success),
                                                'iterations': int(result.nit), 'message': str(result.message)})
            record['status'] = 'all four optimization loops completed'
        except Deadline:
            record['status'] = 'trial deadline in ' + phase + '; best resolved calls retained'
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        for mode, (value, params, details) in valid.items():
            if params is not None:
                record[mode] = {'value': value, 'params': params.tolist(), 'details': details}
        if valid['balanced'][0] < best_rho:
            best_rho, warm = valid['balanced'][0], valid['balanced'][1].copy()
        best_support = min(best_support, valid['asymmetric'][0])
        record['elapsed_seconds'] = time.monotonic() - start_time
        report['cases'].append(record)
        report.update(completed=False, completed_cases=len(report['cases']), elapsed_seconds=time.monotonic() - began,
                      best_rho=best_rho, best_support_value=None if math.isinf(best_support) else best_support)
        save(args.output, report)
        print(json.dumps({k: report[k] for k in ('completed_cases', 'elapsed_seconds', 'best_rho', 'best_support_value')}), flush=True)
        threshold = best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  refinement_queue_exhausted=len(report['cases']) == len(selected),
                  stopping_reason='numerical threshold candidate needs exact replay' if threshold else
                  ('selected refinement queue completed' if len(report['cases']) == len(selected) else 'requested time budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
