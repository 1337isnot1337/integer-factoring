#!/usr/bin/env python3
"""Enumerate all infinity-regular upper-block Laurent series at q16.

For the fixed natural infinity heights, solve the two physical filtration
conditions coefficientwise, together with uniform polynomiality. The
result has two F4 parameters at pole2 and one at pole4. Every complete
natural infinity profile is independently compared with the source seed.
"""
from __future__ import annotations

import argparse
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
from dense_coupled_reflection_frames import digest, save
from laurent_series_coupled_frames import construct
from native_rectangular_smith_f4 import RectangularOracle
from reflection_dual_rail_gaussian import tensor, physical_bit_energies
from reflection_laurent_physical_frames import bell_normalize
from search_dense_coupled_reflection_frames import profile_rho
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import Deadline, alarm
from square_zero_coupled_reflection_frames import product as matmul


def filtered_basis(first, inverse, heights, pole):
    q = len(first)
    positions = [(i, j) for i in range(q // 2) for j in range(q // 2, q)
                 if heights[i] - heights[j] + pole >= 0 and
                 (pole < 3 or (i.bit_count() + j.bit_count()) % 2 == 0)]
    pivots = {}
    equations = []
    for a in range(q):
        for b in range(q):
            if heights[a] - heights[b] <= pole:
                continue
            row = sum((first[a][i] & inverse[j][b]) << k for k, (i, j) in enumerate(positions))
            equations.append(row)
            while row:
                pivot = (row & -row).bit_length() - 1
                if pivot in pivots:
                    row ^= pivots[pivot]
                else:
                    pivots[pivot] = row
                    break
    free = [i for i in range(len(positions)) if i not in pivots]
    basis = []
    for coordinate in free:
        vector = 1 << coordinate
        for pivot in sorted(pivots, reverse=True):
            if (pivots[pivot] & vector).bit_count() % 2:
                vector ^= 1 << pivot
        assert all((row & vector).bit_count() % 2 == 0 for row in equations)
        matrix = [[0] * q for _ in range(q)]
        for k, (i, j) in enumerate(positions):
            matrix[i][j] = (vector >> k) & 1
        conjugate = matmul(matmul(first, matrix), inverse)
        assert all(not matrix[i][j] or heights[i] - heights[j] + pole >= 0 for i in range(q) for j in range(q))
        assert all(not conjugate[i][j] or heights[i] - heights[j] <= pole for i in range(q) for j in range(q))
        basis.append(matrix)
    return {'pole': pole, 'allowed_entries': len(positions), 'constraint_rank': len(pivots),
            'dimension': len(basis), 'basis': basis}


def seed_from_parameters(base, bell_data, basis, coefficients):
    matrices = {}
    cursor = 0
    for item in basis:
        matrix = [[0] * 16 for _ in range(16)]
        for direction in item['basis']:
            c = coefficients[cursor]
            cursor += 1
            for i in range(16):
                for j in range(16):
                    if direction[i][j]:
                        matrix[i][j] ^= c
        if any(c for row in matrix for c in row):
            matrices[item['pole']] = matrix
    return construct(base, matrices, bell_data=bell_data) if matrices else base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--refine', type=int, default=12)
    parser.add_argument('--trial-seconds', type=float, default=65)
    parser.add_argument('--iterations', type=int, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(8)
    bell_data = bell_normalize(base)
    factors = bell_data[1]
    bit_weights = physical_bit_energies(8)
    energy = np.asarray([sum(v for j, v in enumerate(bit_weights) if i >> j & 1) for i in range(16)], float)
    heights = list(map(int, energy / 2))
    basis = [filtered_basis(factors['first_factor'], factors['inverse_factor'], heights, pole) for pole in range(1, 5)]
    assert [r['dimension'] for r in basis] == [0, 2, 0, 1]
    natural = lambda s: start_params(energy * s, s)
    original = RectangularOracle(base['values'], 16, 16)
    baseline_orders = original.prefixes(natural(1.0), 1)['orders']
    names = (Path(__file__).name, 'laurent_series_coupled_frames.py', 'reflection_dual_rail_gaussian.py',
             'reflection_laurent_physical_frames.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp', 'search_dense_coupled_reflection_frames.py')
    report = {'schema_version': 1, 'scope': __doc__, 'requested_seconds': args.seconds,
              'trial_seconds': args.trial_seconds, 'requested_refinements': args.refine, 'iteration_limit': args.iterations,
              'physical_infinity_heights': heights, 'linear_spaces': basis,
              'original_natural_profiles': baseline_orders, 'screen': [], 'refinements': [],
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'counters': {'calls': 0, 'unresolved': 0}}
    began, seen = time.monotonic(), set()
    for coefficients in product(range(4), repeat=3):
        changed = seed_from_parameters(base, bell_data, basis, coefficients)
        oracle = RectangularOracle(changed['values'], 16, 16)
        profiles = oracle.prefixes(natural(1.0), 1)['orders']
        assert profiles[:256] == baseline_orders[:256]
        fit = minimize_scalar(lambda s: profile_rho(profiles, energy, s), bounds=(.005, .4),
                              method='bounded', options={'xatol': 1e-10})
        hashed = digest(changed['values'])
        report['screen'].append({'coefficients': coefficients, 'seed_sha256': hashed,
                                'duplicate_tensor': hashed in seen, 'natural_profiles': profiles,
                                'natural_rho': float(fit.fun), 'natural_scale': float(fit.x)})
        seen.add(hashed)
        report.update(completed=False, screened_frames=len(report['screen']),
                      distinct_tensors=len(seen), elapsed_seconds=time.monotonic() - began)
        save(args.output, report)
        print(json.dumps({'screened_frames': len(report['screen']), 'best_natural_rho': min(r['natural_rho'] for r in report['screen']),
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(screen_complete=True, unchanged_infinity_profile_orders=64 * 256)
    selected = sorted((r for r in report['screen'] if not r['duplicate_tensor']), key=lambda r: r['natural_rho'])[:args.refine]
    best = min(r['natural_rho'] for r in report['screen'])
    for candidate in selected:
        if time.monotonic() - began >= args.seconds:
            break
        changed = seed_from_parameters(base, bell_data, basis, candidate['coefficients'])
        oracle = AsymmetricPressure(changed, 1)
        valid = [math.inf, None, None]
        record = dict(candidate, optimizations=[])
        start_time = time.monotonic()

        def objective(params):
            report['counters']['calls'] += 1
            try:
                value, gradient, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
                value, gradient = value / math.sqrt(2), gradient[:-1] / 2
                if value < valid[0]:
                    valid[:] = [value, np.asarray(params).copy(), details]
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved'] += 1
                return 1e6, np.zeros(len(params))

        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, args.seconds - (start_time - began))))
        try:
            for multiplier in (1, .5):
                result = minimize(objective, natural(candidate['natural_scale'] * multiplier), jac=True,
                                  method='SLSQP', bounds=[(-24, 24)] * 120 + [(.0005, 1)],
                                  options={'maxiter': args.iterations, 'ftol': 1e-10})
                record['optimizations'].append({'success': bool(result.success), 'iterations': int(result.nit),
                                                'evaluations': int(result.nfev), 'message': str(result.message)})
            record['status'] = 'both optimization loops completed'
        except Deadline:
            record['status'] = 'trial deadline; best resolved call retained'
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if valid[1] is not None:
            record['best'] = {'rho': valid[0], 'params': valid[1].tolist(), 'details': valid[2]}
            best = min(best, valid[0])
        record['elapsed_seconds'] = time.monotonic() - start_time
        report['refinements'].append(record)
        report.update(completed=False, elapsed_seconds=time.monotonic() - began, best_rho=best)
        save(args.output, report)
        if best < 1 / math.sqrt(2) - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate needs independent exact replay'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  refinement_queue_exhausted=len(report['refinements']) == len(selected), best_rho=best,
                  stopping_reason=report.get('stopping_reason', 'declared refinement or time budget reached'))
    save(args.output, report)


if __name__ == '__main__':
    main()
