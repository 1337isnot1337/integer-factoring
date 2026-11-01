#!/usr/bin/env python3
"""Bounded search for a two-cut profile below a quarter-circle support line.

The basis maps, polynomial seed, and endpoints are exact over F4. Pressure,
types, heights, and supporting-line weights are numerical. A value below one
requires exact replay and a complete rectangular-recursion transfer.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import tensor
from search_physical_code_gauges import start_params
from search_reflection_physical_bases import canonical_basis, change_basis, generate_basis_pair


def save(path, report):
    """Keep each completed case on one line; source records remain inspectable."""
    metadata = {key: value for key, value in report.items() if key != 'cases'}
    head = json.dumps(metadata, indent=2)
    body = ',\n'.join(json.dumps(case, separators=(',', ':')) for case in report['cases'])
    value = head[:-2] + ',\n  "cases": [\n' + body + '\n  ]\n}\n'
    temporary = path.with_suffix(path.suffix + '.tmp')
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(value)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=7200)
    parser.add_argument('--iterations', type=int, default=200)
    parser.add_argument('--seed', type=int, default=311801)
    parser.add_argument('--case-limit', type=int, default=1000000)
    args = parser.parse_args()
    raw = args.source.read_bytes()
    source = json.loads(raw)
    if source['alphabet'] != 4 or source['interval'] != 4:
        raise ValueError('this search requires the four-letter reflection source')
    base = tensor(4)
    natural = start_params(np.asarray([0, 12, 4, 16]) * 0.2, 0.2)
    rng = random.Random(args.seed)
    noise = np.random.default_rng(args.seed)
    source_cases = sorted(source['cases'], key=lambda record: record['best']['rho'])
    source_cursor = generation = 0
    seen, elites = set(), []
    directory = Path(__file__).parent
    names = (Path(__file__).name, 'asymmetric_crossing_pressure.py',
             'search_reflection_physical_bases.py', 'physical_split_crossings.py',
             'reflection_dual_rail_gaussian.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    report = {'schema_version': 1, 'scope': 'numerical support-line search; no exact pressure or factoring claim',
              'seed': args.seed, 'requested_seconds': args.seconds,
              'source_snapshot_sha256': sha256(raw).hexdigest(), 'source_cases_available': len(source_cases),
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'angle_range': [0.03, math.pi / 2 - 0.03], 'iteration_limit': args.iterations,
              'counters': {'objective_calls': 0, 'unresolved_calls': 0}, 'cases': []}
    started, best_case = time.monotonic(), None
    while time.monotonic() - started < args.seconds and len(report['cases']) < args.case_limit:
        if source_cursor < len(source_cases):
            reference = source_cases[source_cursor]
            source_cursor += 1
            label, left, right = 'resolved-isotropic-source', reference['basis_x'], reference['basis_y']
            start = np.asarray(reference['best']['params'], float) / math.sqrt(2)
            cut = reference['best']['cut']
            source_hash = sha256(json.dumps(reference, sort_keys=True).encode()).hexdigest()
        else:
            label, left, right = generate_basis_pair(generation, 4, rng, elites)
            generation += 1
            start, cut, source_hash = natural / math.sqrt(2), 1, None
        key = (canonical_basis(left), canonical_basis(right), cut)
        if key in seen:
            continue
        seen.add(key)
        seed = change_basis(base, left, right)
        oracle = AsymmetricPressure(seed, cut)
        angles = [math.pi / 4, 0.3, math.pi / 2 - 0.3]
        if len(report['cases']) % 4:
            angles = [0.12 + noise.random() * 1.3, math.pi / 4, 1.3]
        runs, local_best = [], None
        for angle in angles:
            initial = np.r_[start, angle]
            valid = [math.inf, None, None]

            def objective(variables):
                report['counters']['objective_calls'] += 1
                try:
                    value, gradient, details = oracle.evaluate(variables)
                    if value < valid[0]:
                        valid[:] = [value, np.asarray(variables, float).copy(), details]
                    return value, gradient
                except (ArithmeticError, ValueError):
                    report['counters']['unresolved_calls'] += 1
                    return 1000000.0, np.zeros(len(variables))

            initial_value = objective(initial)[0]
            outcome = minimize(objective, initial, jac=True, method='SLSQP',
                               bounds=[(-24, 24)] * 24 + [(0.0005, 2), (0.03, math.pi / 2 - 0.03)],
                               options={'maxiter': args.iterations, 'ftol': 1e-11})
            runs.append({'initial_angle': float(angle), 'initial_value': float(initial_value),
                         'best_value': None if valid[1] is None else float(valid[0]),
                         'iterations': int(outcome.nit), 'evaluations': int(outcome.nfev),
                         'success': bool(outcome.success), 'message': str(outcome.message)})
            if valid[1] is not None and (local_best is None or valid[0] < local_best['support_value']):
                local_best = dict(valid[2], variables=valid[1].tolist())
        if local_best is None:
            record = {'kind': label, 'basis_x': left, 'basis_y': right, 'cut': cut,
                      'status': 'no resolved complete profile', 'runs': runs}
        else:
            record = {'kind': label, 'basis_x': left, 'basis_y': right, 'cut': cut,
                      'status': 'resolved', 'best': local_best, 'runs': runs,
                      'source_case_sha256': source_hash, 'endpoint_entries': 256}
            if best_case is None or local_best['support_value'] < best_case['best']['support_value']:
                best_case = record
            elites.append(record)
            elites.sort(key=lambda r: r['best']['support_value'])
            del elites[12:]
        report['cases'].append(record)
        report.update(completed=False, completed_cases=len(report['cases']),
                      elapsed_seconds=time.monotonic() - started,
                      best_support_value=None if best_case is None else best_case['best']['support_value'],
                      best_candidate=best_case)
        save(args.output, report)
        print(json.dumps({'cases': len(report['cases']), 'best_support_value': report['best_support_value'],
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        if best_case is not None and best_case['best']['support_value'] < 1 - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate needs independent exact replay'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'requested time or case budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
