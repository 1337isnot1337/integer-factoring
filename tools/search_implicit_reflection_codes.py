#!/usr/bin/env python3
"""Bounded exact-code construction and complete natural-gauge pressure search.

All choices at the final natural-energy level are kept. The tensor and
scalar-content removal are exact; the optimized pressure is numerical.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar

from affine_place_profile_oracle import primitive_seed
from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_affine_coupled_reflection_frames import balanced
from search_fused_mixed_weighted import SearchDeadline, atomic_write
from search_physical_code_gauges import start_params


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=16)
    parser.add_argument('--sizes', type=int, nargs='+', default=[8, 12, 16, 20, 24])
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--trial-seconds', type=float, default=120)
    parser.add_argument('--iterations', type=int, default=40)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    weights = physical_bit_energies(args.interval)
    energy = lambda symbol: sum(weight for bit, weight in enumerate(weights) if symbol >> bit & 1)
    labels = sorted(range(1 << (args.interval // 2)), key=lambda symbol: (energy(symbol), symbol))
    cases = []
    for size in args.sizes:
        if not 2 <= size <= len(labels):
            raise ValueError('Each code size must be between two and the ambient alphabet.')
        level = energy(labels[size - 1])
        lower = [symbol for symbol in labels if energy(symbol) < level]
        boundary = [symbol for symbol in labels if energy(symbol) == level]
        for chosen in combinations(boundary, size - len(lower)):
            cases.append((size, level, lower + list(chosen)))
    names = (Path(__file__).name, 'restrict_implicit_reflection.py',
             'reflection_dual_rail_gaussian.py', 'check_schur_gaussian_reduction.py',
             'affine_place_profile_oracle.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_affine_coupled_reflection_frames.py', 'search_physical_code_gauges.py')
    report = {'scope': __doc__, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'arguments': {key: value for key, value in vars(args).items() if key != 'output'},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'expected_cases': len(cases), 'records': [], 'completed': False}
    began, deadline = time.monotonic(), time.monotonic() + args.seconds
    best_value = float('inf')
    for size, level, letters in cases:
        if time.monotonic() >= deadline:
            break
        started = time.monotonic()
        raw = restricted_tensor(args.interval, letters)
        seed, common = primitive_seed(raw)
        record = {'interval': args.interval, 'size': size, 'level': level, 'letters': letters,
                  'energies': [energy(v) for v in letters], 'raw_degree': raw['degree'],
                  'primitive_degree': seed['degree'], 'scalar_content': common,
                  'tensor_sha256': digest(seed['values']), 'construction_counts': raw['construction_counts'],
                  'construction_seconds': time.monotonic() - started, 'pressure_queries': [],
                  'completed': False, 'unresolved': 0}
        report['records'].append(record)
        oracle = AsymmetricPressure(seed, 1)
        trial_deadline = min(deadline, time.monotonic() + args.trial_seconds)

        def objective(scale):
            nonlocal best_value
            if time.monotonic() >= trial_deadline:
                raise SearchDeadline()
            params = start_params(np.asarray(record['energies']) * scale, scale)
            try:
                value, _, details = balanced(oracle, params)
            except (ArithmeticError, ValueError):
                record['unresolved'] += 1
                return 1e6
            record['pressure_queries'].append({'scale': float(scale), 'rho': float(value)})
            if 'best' not in record or value < record['best']['rho']:
                record['best'] = {'rho': float(value), 'scale': float(scale),
                                  'params': params.tolist(), 'details': details}
            if value < best_value:
                best_value = value
                report['best_record_index'] = len(report['records']) - 1
            report['elapsed_seconds'] = time.monotonic() - began
            atomic_write(args.output, report)
            return value

        try:
            fit = minimize_scalar(objective, method='bounded', bounds=(.0001, 1),
                                  options={'xatol': 1e-8, 'maxiter': args.iterations})
            record['optimization'] = {'success': bool(fit.success), 'calls': int(fit.nfev),
                                      'message': str(fit.message)}
            record['completed'] = True
        except SearchDeadline:
            record['stopping_reason'] = 'Declared local or total deadline reached; all resolved calls retained.'
        record['elapsed_seconds'] = time.monotonic() - started
        report['elapsed_seconds'] = time.monotonic() - began
        atomic_write(args.output, report)
        print(json.dumps({'size': size, 'letters': letters, 'degree': seed['degree'],
                          'best_rho': record.get('best', {}).get('rho'),
                          'complete': record['completed'], 'elapsed': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - began,
                  case_list_exhausted=len(report['records']) == len(cases),
                  completed_cases=sum(record['completed'] for record in report['records']))
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
