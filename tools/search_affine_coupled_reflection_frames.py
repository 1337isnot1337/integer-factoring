#!/usr/bin/env python3
"""Bounded affine sums of pole-cancelling four-letter reflection frames.

T + sum c_i (T_i - T) keeps the strict endpoint and leading coefficient.
The resulting path need not stay in the physical change-of-basis orbit.
Every pressure call uses the complete resulting polynomial coefficients.
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
from reflection_laurent_physical_frames import (
    add, bell_normalize, bell_restore, coupled_admissible_poles, coupled_terms,
    coupled_transformed, degree, scale, valuation,
)
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import Deadline, alarm, save


def catalogue(base):
    bell, factors = bell_normalize(base)
    records = []
    for target in range(4):
        for source in range(4):
            if target == source:
                continue
            for scalar in (1, 2, 3):
                terms = coupled_terms(bell, target, source, scalar)
                for pole in coupled_admissible_poles(terms, base['degree']):
                    transformed = coupled_transformed(bell, target, source, -pole, scalar, terms)
                    changed = bell_restore(transformed, factors)
                    delta = [add(a, b) for a, b in zip(base['values'], changed['values'])]
                    assert max(degree(v) for v in delta) < base['degree']
                    records.append({'target': target, 'source': source, 'scalar': scalar,
                                    'pole': pole, 'delta': delta})
    return records


def affine_seed(base, frames, word):
    values = list(base['values'])
    for index, coefficient in word:
        values = [add(a, scale(b, coefficient)) for a, b in zip(values, frames[index]['delta'])]
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('Affine sum changed a crossing endpoint entry.')
    actual = max(degree(value) for value in values)
    if actual != base['degree']:
        raise ArithmeticError('Affine sum changed the common leading coefficient.')
    return dict(base, values=values, degree=actual,
                monomial_content=min(valuation(v) for v in values if v != (0, 0)))


def balanced(oracle, params):
    rate, gradient, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
    return rate / math.sqrt(2), gradient[:-1] / 2, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--screen', type=int, default=12000)
    parser.add_argument('--refine', type=int, default=120)
    parser.add_argument('--trial-seconds', type=float, default=20)
    parser.add_argument('--iterations', type=int, default=260)
    parser.add_argument('--seed', type=int, default=312331)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.seconds, args.screen, args.refine, args.trial_seconds, args.iterations) <= 0:
        parser.error('positive budgets are required')
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(4)
    frames = catalogue(base)
    baseline = AsymmetricPressure(base, 1)
    natural = lambda s: start_params(np.asarray([0, 12, 4, 16]) * s, s)
    optimal_scale = minimize_scalar(lambda s: balanced(baseline, natural(s))[0], bounds=(.005, 1),
                                    method='bounded', options={'xatol': 1e-10}).x
    initial = natural(optimal_scale)
    baseline_rho = balanced(baseline, initial)[0]
    rng = np.random.default_rng(args.seed)
    names = (Path(__file__).name, 'reflection_laurent_physical_frames.py',
             'reflection_dual_rail_gaussian.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_reflection4_coupled_frames.py', 'search_physical_code_gauges.py')
    report = {'schema_version': 1, 'scope': 'exact affine strict paths; bounded numerical pressure search',
              'seed': args.seed, 'requested_seconds': args.seconds, 'requested_screens': args.screen,
              'requested_refinements': args.refine, 'trial_seconds': args.trial_seconds,
              'iteration_limit': args.iterations, 'baseline_rho': baseline_rho,
              'frame_catalogue': [{k: v for k, v in entry.items() if k != 'delta'} for entry in frames],
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'screen_records': [], 'cases': [], 'proposals_attempted': 0,
              'counters': {'calls': 0, 'unresolved': 0, 'duplicate_proposals': 0}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    started = time.monotonic()
    for index in range(args.screen):
        if time.monotonic() - started >= args.seconds * .4:
            break
        report['proposals_attempted'] = index + 1
        length = 2 if index % 3 else int(rng.integers(3, 6))
        chosen = rng.choice(len(frames), length, replace=False).tolist()
        word = sorted((int(i), int(rng.integers(1, 4))) for i in chosen)
        changed = affine_seed(base, frames, word)
        digest = sha256(json.dumps(changed['values'], separators=(',', ':')).encode()).hexdigest()
        if digest in seen:
            report['counters']['duplicate_proposals'] += 1
            continue
        seen.add(digest)
        oracle = AsymmetricPressure(changed, 1)
        best = math.inf
        for factor in (1, .5):
            report['counters']['calls'] += 1
            try:
                best = min(best, balanced(oracle, initial * factor)[0])
            except (ArithmeticError, ValueError):
                report['counters']['unresolved'] += 1
        report['screen_records'].append({'word': word, 'seed_sha256': digest,
                                         'degree': changed['degree'], 'monomial_content': changed['monomial_content'],
                                         'natural_rho': None if math.isinf(best) else best})
        if len(report['screen_records']) % 200 == 0:
            report.update(completed=False, elapsed_seconds=time.monotonic() - started,
                          completed_screens=len(report['screen_records']))
            save(args.output, report)
            print(json.dumps({'screened': len(report['screen_records']), 'elapsed_seconds': report['elapsed_seconds'],
                              'best_natural_rho': min(r['natural_rho'] for r in report['screen_records'] if r['natural_rho'] is not None)}), flush=True)
    selected = sorted((r for r in report['screen_records'] if r['natural_rho'] is not None),
                      key=lambda r: r['natural_rho'])[:args.refine]
    report['screen_phase_complete'] = report['proposals_attempted'] == args.screen
    report['completed_screens'] = len(report['screen_records'])
    best_rho, best_support = baseline_rho, math.inf
    warm = initial.copy()
    threshold = False
    for candidate in selected:
        if time.monotonic() - started >= args.seconds or threshold:
            break
        changed = affine_seed(base, frames, candidate['word'])
        oracle = AsymmetricPressure(changed, 1)
        phase = 'balanced'
        valid = {'balanced': [math.inf, None, None], 'asymmetric': [math.inf, None, None]}
        record = dict(candidate, optimizations=[])
        begin = time.monotonic()

        def objective(params):
            report['counters']['calls'] += 1
            try:
                value, gradient, detail = balanced(oracle, params) if phase == 'balanced' else oracle.evaluate(params)
                if value < valid[phase][0]:
                    valid[phase][:] = [value, np.asarray(params).copy(), detail]
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved'] += 1
                return 1e6, np.zeros(len(params))

        budget = min(args.trial_seconds, max(.01, args.seconds - (begin - started)))
        signal.setitimer(signal.ITIMER_REAL, budget)
        try:
            for start in (warm.copy(), initial * .25):
                fit = minimize(objective, start, jac=True, method='SLSQP',
                               bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'success': bool(fit.success),
                                                'iterations': int(fit.nit), 'evaluations': int(fit.nfev),
                                                'message': str(fit.message)})
            phase = 'asymmetric'
            start = valid['balanced'][1] if valid['balanced'][1] is not None else initial
            for angle in (.35, 1.22):
                fit = minimize(objective, np.r_[start / math.sqrt(2), angle], jac=True, method='SLSQP',
                               bounds=[(-24, 24)] * 24 + [(.0005, 2), (.025, math.pi / 2 - .025)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'success': bool(fit.success),
                                                'iterations': int(fit.nit), 'evaluations': int(fit.nfev),
                                                'message': str(fit.message)})
            record['status'] = 'all four optimization loops completed'
        except Deadline:
            record['status'] = 'trial deadline in ' + phase + '; best resolved calls retained'
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        for mode, (value, params, detail) in valid.items():
            if params is not None:
                record[mode] = {'value': value, 'params': params.tolist(), 'details': detail}
        if valid['balanced'][0] < best_rho:
            best_rho, warm = valid['balanced'][0], valid['balanced'][1].copy()
        best_support = min(best_support, valid['asymmetric'][0])
        record['elapsed_seconds'] = time.monotonic() - begin
        report['cases'].append(record)
        threshold = best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6
        report.update(completed=False, completed_cases=len(report['cases']),
                      elapsed_seconds=time.monotonic() - started, best_rho=best_rho,
                      best_support_value=None if math.isinf(best_support) else best_support)
        save(args.output, report)
        print(json.dumps({k: report[k] for k in ('completed_cases', 'elapsed_seconds', 'best_rho', 'best_support_value')}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  refinement_queue_exhausted=len(report['cases']) == len(selected),
                  stopping_reason='threshold candidate requires independent exact replay' if threshold else
                  ('selected refinement queue completed' if len(report['cases']) == len(selected) else 'requested time budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
