#!/usr/bin/env python3
"""Bounded rank-two square-zero deformation search with complete profile costs.

The square-zero matrices, polynomiality, degree, and endpoints are exact.
The pressure optimization is exploratory. Every resolved best call is kept;
timeouts and missing native pivots are explicit and never certify a rate.
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
from physical_split_crossings import identity_split
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from search_physical_code_gauges import start_params
from search_reflection_physical_bases import inverse, row_shear
from search_reflection4_coupled_frames import Deadline, alarm, save
from square_zero_coupled_reflection_frames import (
    construct, first_jet_identity, matching_generator, product, rank,
    square_zero_matrix,
)


def generators(rng):
    """Structured rank-two cases followed by exact random conjugates."""
    for permutation in ([1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]):
        yield 'perfect-matching', matching_generator(permutation)
    jordan = [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]]
    yield 'two-jordan-blocks', jordan
    # Each summand stays inside one parity eigenspace of the first jet.
    for even_direction in ((0, 3), (3, 0)):
        for odd_direction in ((1, 2), (2, 1)):
            for first_scalar in (1, 2, 3):
                for second_scalar in (1, 2, 3):
                    g = [[0] * 4 for _ in range(4)]
                    g[even_direction[0]][even_direction[1]] = first_scalar
                    g[odd_direction[0]][odd_direction[1]] = second_scalar
                    yield 'parity-blocks', g
    index = 0
    while True:
        if index % 2:
            basis = identity_split(4)
            for _ in range(1 + index % 9):
                i, j = rng.choice(4, 2, replace=False)
                basis = row_shear(basis, int(i), int(j), int(rng.integers(1, 4)))
        else:
            while True:
                basis = rng.integers(0, 4, (4, 4)).tolist()
                try:
                    inverse(basis)
                    break
                except ValueError:
                    pass
        yield 'conjugated-rank-two', product(product(basis, jordan), inverse(basis))
        index += 1


def balanced(oracle, params):
    value, derivative, detail = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
    return value / math.sqrt(2), derivative[:-1] / 2, detail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--trial-seconds', type=float, default=10)
    parser.add_argument('--iterations', type=int, default=240)
    parser.add_argument('--seed', type=int, default=312244)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.seconds, args.trial_seconds, args.iterations) <= 0:
        parser.error('positive budgets are required')
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(4)
    bell_data = bell_normalize(base)
    controls = first_jet_identity(base, bell_data)
    native = AsymmetricPressure(base, 1)
    natural = lambda s: start_params(np.asarray([0, 12, 4, 16]) * s, s)
    scale = minimize_scalar(lambda s: balanced(native, natural(s))[0], bounds=(.005, 1),
                            method='bounded', options={'xatol': 1e-10}).x
    initial = natural(scale)
    baseline, _, baseline_detail = balanced(native, initial)
    rng = np.random.default_rng(args.seed)
    names = (Path(__file__).name, 'square_zero_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_reflection4_coupled_frames.py', 'search_reflection_physical_bases.py',
             'search_physical_code_gauges.py', 'physical_split_crossings.py',
             'reflection_dual_rail_gaussian.py')
    report = {'schema_version': 1, 'scope': 'bounded rank-two square-zero frames; numerical pressure only',
              'seed': args.seed, 'requested_seconds': args.seconds,
              'trial_seconds': args.trial_seconds, 'iteration_limit': args.iterations,
              'first_jet_controls': controls, 'baseline_rho': baseline,
              'baseline_details': baseline_detail, 'cases': [], 'admissibility': [],
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'counters': {'calls': 0, 'unresolved': 0, 'rejected_poles': 0}}
    best_rho, best_support = baseline, math.inf
    seen = set()
    warm = initial.copy()
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    threshold = False
    for label, generator in generators(rng):
        if time.monotonic() - started >= args.seconds or threshold:
            break
        key = tuple(x for row in generator for x in row)
        if key in seen:
            continue
        seen.add(key)
        assert square_zero_matrix(generator) and rank(generator) == 2
        candidates = []
        entry = {'generator': generator, 'family': label, 'tested_poles': [1, 2, 3, 4], 'allowed_poles': []}
        for pole in (1, 2, 3, 4):
            try:
                changed = construct(base, generator, pole=pole, bell_data=bell_data)
                candidates.append((pole, changed))
                entry['allowed_poles'].append(pole)
            except ArithmeticError:
                report['counters']['rejected_poles'] += 1
                if pole == 1:
                    raise
        report['admissibility'].append(entry)
        for pole, seed in candidates:
            if time.monotonic() - started >= args.seconds or threshold:
                break
            oracle = AsymmetricPressure(seed, 1)
            phase = 'balanced' if len(report['cases']) % 3 != 2 else 'asymmetric'
            valid = [math.inf, None, None]
            record = dict(entry, pole=pole, phase=phase, degree=seed['degree'],
                          monomial_content=seed['monomial_content'], endpoint_entries=256,
                          seed_sha256=sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest(),
                          optimizations=[])
            begin = time.monotonic()
            before = report['counters'].copy()

            def objective(params):
                report['counters']['calls'] += 1
                try:
                    value, gradient, detail = balanced(oracle, params) if phase == 'balanced' else oracle.evaluate(params)
                    if value < valid[0]:
                        valid[:] = [value, np.asarray(params).copy(), detail]
                    return value, gradient
                except (ArithmeticError, ValueError):
                    report['counters']['unresolved'] += 1
                    return 1e6, np.zeros(len(params))

            budget = min(args.trial_seconds, max(.01, args.seconds - (begin - started)))
            signal.setitimer(signal.ITIMER_REAL, budget)
            try:
                starts = [warm.copy(), initial * .25,
                          np.r_[rng.normal(0, .08, 24), scale * .25]]
                if phase == 'asymmetric':
                    starts = [np.r_[value / math.sqrt(2), angle]
                              for value, angle in zip(starts, (.4, 1.17, math.pi / 4))]
                bounds = [(-24, 24)] * 24 + [(.0005, 2)]
                if phase == 'asymmetric':
                    bounds += [(.025, math.pi / 2 - .025)]
                for start in starts:
                    fit = minimize(objective, start, jac=True, method='SLSQP', bounds=bounds,
                                   options={'maxiter': args.iterations, 'ftol': 2e-11})
                    record['optimizations'].append({'success': bool(fit.success), 'iterations': int(fit.nit),
                                                    'evaluations': int(fit.nfev), 'message': str(fit.message)})
                record['status'] = 'all three optimization loops completed'
            except Deadline:
                record['status'] = 'trial deadline; best resolved calls retained'
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            if valid[1] is not None:
                record['best'] = {'value': valid[0], 'params': valid[1].tolist(), 'details': valid[2]}
                if phase == 'balanced' and valid[0] < best_rho:
                    best_rho, warm = valid[0], valid[1].copy()
                if phase == 'asymmetric':
                    best_support = min(best_support, valid[0])
            record['elapsed_seconds'] = time.monotonic() - begin
            record['counters'] = {k: report['counters'][k] - before[k] for k in before}
            report['cases'].append(record)
            threshold = best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6
            report.update(completed=False, completed_cases=len(report['cases']),
                          elapsed_seconds=time.monotonic() - started, best_rho=best_rho,
                          best_support_value=None if math.isinf(best_support) else best_support)
            save(args.output, report)
            print(json.dumps({k: report[k] for k in ('completed_cases', 'elapsed_seconds', 'best_rho', 'best_support_value')}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  stopping_reason='threshold candidate requires independent exact replay' if threshold else 'requested time budget completed')
    save(args.output, report)


if __name__ == '__main__':
    main()
