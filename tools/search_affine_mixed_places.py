#!/usr/bin/env python3
"""Bounded full-cost mixed-place searches of frozen affine reflection words.

All cuts come from reconstructed complete polynomial tensors.  A single
pair of physical types is shared by every place.  The degree charge,
physical gauge means, and all local Smith channels are included.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from affine_place_profile_oracle import AffinePlaceOracle, primitive_seed
from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import tensor
from search_affine_coupled_reflection_frames import affine_seed, catalogue
from search_fused_mixed_weighted import SearchDeadline, atomic_write
from search_fused_physical_bases import ScaledHeightPressure
from search_physical_code_gauges import start_params


PARTITIONS = (([0], [None]), ([0], [None, 2, 3]), ([0, 2, 3], [None]),
              ([0, 2], [None, 3]), ([0, 3], [None, 2]), ([0], [2, 3]),
              ([2, 3], [None]), ([0], [None, 2]), ([0], [None, 3]),
              ([0, 2], [None]), ([0, 3], [None]))


def baseline(partition, scale):
    row, column = partition

    def allocations(places):
        if len(places) == 1:
            return [1.0]
        if 0 in places or None in places:
            return [0.7] + [0.3 / (len(places) - 1)] * (len(places) - 1)
        return [1 / len(places)] * len(places)

    return {'row_places': row, 'column_places': column,
            'row_allocations': allocations(row), 'column_allocations': allocations(column),
            'outer_log_scale': scale}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rank_screen', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--trial-seconds', type=float, default=30)
    parser.add_argument('--maxiter', type=int, default=160)
    parser.add_argument('--candidates', type=int, default=24)
    parser.add_argument('--trials', type=int, default=1000000)
    parser.add_argument('--seed', type=int, default=312351)
    args = parser.parse_args()
    if min(args.seconds, args.trial_seconds, args.maxiter, args.candidates, args.trials) <= 0:
        parser.error('All budgets must be positive.')
    raw = args.rank_screen.read_bytes()
    screen = json.loads(raw)
    directory = Path(__file__).parent
    for name, expected in screen['source_sha256'].items():
        if sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ArithmeticError('A rank-screen source changed: ' + name)
    base = tensor(4)
    frames = catalogue(base)
    ranked = sorted(screen['retained'], key=lambda row: (min(v for pair in row['ranks'] for v in pair),
                                                        sum(v for pair in row['ranks'] for v in pair),
                                                        row['natural_rho'] or math.inf))
    natural = sorted(screen['retained'], key=lambda row: row['natural_rho'] or math.inf)
    selected, seen = [], set()
    for rows in (ranked[:args.candidates], natural[:max(4, args.candidates // 2)]):
        for row in rows:
            if row['seed_sha256'] not in seen:
                selected.append(row)
                seen.add(row['seed_sha256'])
    # Keep an undeformed two-place baseline in the same adapter.
    selected.insert(0, {'word': [], 'seed_sha256': sha256(json.dumps(base['values'], separators=(',', ':')).encode()).hexdigest(),
                        'ranks': [[16, 16], [16, 16]], 'natural_rho': None, 'snapshot_index': None})
    original = AsymmetricPressure(base, 1)
    energy = np.asarray([0, 12, 4, 16], float)

    def natural_parameters(scale):
        return np.r_[start_params(energy * scale, scale) / math.sqrt(2), math.pi / 4]

    fit = minimize_scalar(lambda scale: original.evaluate(natural_parameters(scale), gradient=False)[0],
                          bounds=(.005, 1), method='bounded', options={'xatol': 1e-10})
    _, type_data = original.evaluate(natural_parameters(float(fit.x)), gradient=False)
    base_type = [type_data['type_x'], type_data['type_y']]
    rng = np.random.default_rng(args.seed)
    names = (Path(__file__).name, 'affine_place_profile_oracle.py', 'search_fused_mixed_weighted.py',
             'search_fused_physical_bases.py', 'fused_physical_basis_oracle.py',
             'search_affine_coupled_reflection_frames.py', 'reflection_laurent_physical_frames.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Bounded floating full-cost multipoint search; no global optimality or factoring certificate.',
              'arguments': vars(args) | {'rank_screen': str(args.rank_screen), 'output': str(args.output)},
              'rank_screen_sha256': sha256(raw).hexdigest(), 'selected_candidates': selected,
              'baseline_rho': float(fit.fun) / math.sqrt(2), 'baseline_log_scale': float(fit.x),
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'records': [], 'winners': {}, 'budget_completed': False, 'attempted_trials': 0, 'completed_trials': 0}
    oracles = {}
    started = time.monotonic()
    deadline = started + args.seconds
    for trial in range(args.trials):
        if time.monotonic() >= deadline:
            break
        candidate_index = trial % len(selected)
        round_index = trial // len(selected)
        partition_index = round_index % len(PARTITIONS)
        row = selected[candidate_index]
        if candidate_index not in oracles:
            source = affine_seed(base, frames, row['word'])
            digest = sha256(json.dumps(source['values'], separators=(',', ':')).encode()).hexdigest()
            if digest != row['seed_sha256']:
                raise ArithmeticError('Rebuilt full tensor disagrees with its frozen hash.')
            primitive, content = primitive_seed(source)
            oracles[candidate_index] = AffinePlaceOracle(primitive), content
        oracle, content = oracles[candidate_index]
        current = baseline(PARTITIONS[partition_index], float(fit.x))
        begin = time.monotonic()
        pressure = ScaledHeightPressure(oracle, current, 'full', False, True,
                                        min(deadline, begin + args.trial_seconds))
        values = pressure.initial(current, rng, 0)
        for cut, groups in enumerate(pressure.groups):
            for index, group in enumerate(groups):
                heights = np.zeros(16)
                if group['point'] == (0, 1):
                    heights[0:4] = heights[4:8] = energy / 2
                elif group['point'] is None:
                    heights[8:12] = heights[12:16] = energy / 2
                elif round_index >= len(PARTITIONS):
                    heights = rng.normal(0, .3, 16)
                values[pressure.height_slices[cut][index]] = heights[pressure.active_heights]
        values[pressure.type_slice] = np.concatenate([np.log(np.asarray(p)[1:] / p[0]) for p in base_type])
        bounds = pressure.bounds()
        for slices in pressure.height_slices:
            for current_slice in slices:
                bounds[current_slice] = [(-48, 48)] * (current_slice.stop - current_slice.start)
        outcomes, completed, interruption = [], True, None
        try:
            pressure.evaluate(values)
            for smoothing in (.025, .003, .0003):
                fit_local = minimize(lambda x: pressure.evaluate(x, smoothing), values, jac=True,
                                     method='L-BFGS-B', bounds=bounds,
                                     options={'maxiter': args.maxiter, 'maxfun': 1600,
                                              'ftol': 2e-11, 'gtol': 3e-7, 'maxls': 25})
                outcomes.append({'smoothing': smoothing, 'success': bool(fit_local.success),
                                 'iterations': int(fit_local.nit), 'evaluations': int(fit_local.nfev),
                                 'message': str(fit_local.message)})
                values = np.asarray(pressure.best['optimizer_parameters'])
        except SearchDeadline:
            completed, interruption = False, 'Declared per-trial or global deadline reached.'
        except (ArithmeticError, ValueError) as error:
            completed, interruption = False, f'{type(error).__name__}: {error}'
        record = {'trial': trial, 'candidate_index': candidate_index, 'word': row['word'],
                  'source_sha256': row['seed_sha256'], 'partition_index': partition_index,
                  'row_places': current['row_places'], 'column_places': current['column_places'],
                  'degree': oracle.degree, 'common_factor': content, 'queries': pressure.queries,
                  'elapsed_seconds': time.monotonic() - begin, 'completed': completed,
                  'outcomes': outcomes, 'interruption': interruption,
                  'best': {key: value for key, value in (pressure.best or {}).items() if key != 'place_records'}}
        report['records'].append(record)
        key = str(candidate_index)
        if pressure.best and (key not in report['winners'] or pressure.best['rho'] < report['winners'][key]['best']['rho']):
            report['winners'][key] = record | {'best': pressure.best}
        report.update(elapsed_seconds=time.monotonic() - started, attempted_trials=len(report['records']),
                      completed_trials=sum(row['completed'] for row in report['records']),
                      budget_completed=time.monotonic() >= deadline)
        atomic_write(args.output, report)
        print(json.dumps({'trial': trial, 'candidate': candidate_index, 'partition': partition_index,
                          'rho': pressure.best['rho'] if pressure.best else None,
                          'queries': pressure.queries, 'completed': completed,
                          'elapsed': report['elapsed_seconds']}), flush=True)
    report.update(elapsed_seconds=time.monotonic() - started,
                  budget_completed=time.monotonic() >= deadline,
                  trial_limit_completed=len(report['records']) == args.trials)
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
