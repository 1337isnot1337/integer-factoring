#!/usr/bin/env python3
"""Full-gauge closed-place refinement of exact paired restrictions.

Candidates are selected by complete adaptive unweighted profiles and by
the independently recorded two-place weighted pressure. Every selected
seed is explicit; all conjugate capacities, physical types, and degrees
are charged. Results are numerical search records, not certificates.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import log
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from f4_closed_place_oracle import F4ClosedPlaceOracle
from paired_mixed_place_pressure import ScaledPairedMixedPressure
from reflection_laurent_physical_frames import digest
from screen_restricted_closed_partitions import place_label
from search_fused_mixed_weighted import SearchDeadline, atomic_write, point_key


def selected_cases(catalogue, adaptive, weighted, limit):
    by_hash = {}
    for row in weighted['records']:
        if row.get('best') and (row['tensor_sha256'] not in by_hash or
                               row['best']['value'] < by_hash[row['tensor_sha256']]['best']['value']):
            by_hash[row['tensor_sha256']] = row
    eligible = [row for row in adaptive['records'] if any(point['factor'] != [0, 1]
                for point in catalogue['records'][row['catalogue_index']]['closed_points'])]
    ranked = sorted(eligible, key=lambda row: row['best']['rho'])
    by_weight = sorted(eligible, key=lambda row: by_hash[row['tensor_sha256']]['best']['value'])
    chosen, seen = [], set()
    for index in range(max(len(ranked), len(by_weight))):
        for rows in (ranked, by_weight):
            if index >= len(rows):
                continue
            row = rows[index]
            if row['tensor_sha256'] not in seen:
                chosen.append((row, catalogue['records'][row['catalogue_index']], by_hash[row['tensor_sha256']]))
                seen.add(row['tensor_sha256'])
                if len(chosen) >= limit:
                    return chosen
    return chosen


def partition_cases(adaptive, catalogue, old):
    qx, qy = len(catalogue['subset_x']), len(catalogue['subset_y'])
    n = qx*qy
    points = {tuple(row['factor']): row for row in catalogue['closed_points']}
    proposals = []
    # The complete adaptive partition and its active groups guide one
    # candidate. Removing a zero-capacity group is allowed explicitly.
    best = adaptive['best']
    cuts = []
    for cut, key, rate in ((0, 'row_places', best['row_rate']), (1, 'column_places', best['column_rate'])):
        labels = []
        for label in best[key]:
            factor = point_key(label)
            profile = catalogue['homogeneous_infinity_orders'][cut] if factor is None else points[factor]['orders'][cut]
            if log(profile.count(0)) <= rate+1e-6:
                labels.append(label)
        cuts.append(labels or best[key])
    proposals.append(('adaptive_partition', cuts[0], cuts[1]))
    # Natural zero/infinity orientation with one useful new factor at a
    # time avoids forcing a high-floor factor into every numerical trial.
    extras = []
    for point in catalogue['closed_points']:
        factor = tuple(point['factor'])
        if factor == (0, 1):
            continue
        for cut in (0, 1):
            profile = point['orders'][cut]
            if max(profile):
                score = (profile.count(0)/n, -point['residue_degree']*sum(profile)/n)
                extras.append((score, cut, factor, point['residue_degree']))
    extras.sort()
    best_each = {}
    for _, cut, factor, multiplicity in extras[:4]:
        labels = [place_label(factor, i) for i in range(multiplicity)]
        proposals.append((f'single_factor_cut{cut}_{factor}', [0]+(labels if cut == 0 else []),
                          [None]+(labels if cut == 1 else [])))
        best_each.setdefault(cut, (factor, labels))
    if len(best_each) == 2 and best_each[0][0] != best_each[1][0]:
        proposals.append(('two_new_factors', [0]+best_each[0][1], [None]+best_each[1][1]))
    result, seen = [], set()
    for label, row, column in proposals:
        key = tuple(map(str, row)), tuple(map(str, column))
        if key in seen or not set(key[0]).isdisjoint(key[1]):
            continue
        seen.add(key)
        allocations = []
        for places, preferred in ((row, 0), (column, None)):
            if preferred in places and len(places) > 1:
                allocations.append([.8 if place == preferred else .2/(len(places)-1) for place in places])
            else:
                allocations.append([1/len(places)]*len(places))
        result.append({'label': label, 'row_places': row, 'column_places': column,
                       'row_allocations': allocations[0], 'column_allocations': allocations[1],
                       'outer_log_scale': old['best']['parameters'][-1]})
    return result


def warm_start(pressure, baseline, old, rng, restart):
    values = pressure.initial(baseline, rng, restart)
    if restart:
        return values
    params = np.asarray(old['best']['parameters'])
    height_size = pressure.height_size
    assert len(params) == 2*height_size+1
    scale = params[-1]
    values[pressure.scale_index] = scale
    for cut, groups in enumerate(pressure.groups):
        for index, group in enumerate(groups):
            if cut == 0 and group['point'] == (0, 1):
                values[pressure.height_slices[cut][index]] = params[height_size:2*height_size]/(scale*group['allocation'])
            elif cut == 1 and group['point'] is None:
                values[pressure.height_slices[cut][index]] = params[:height_size]/(scale*group['allocation'])
    details = old['best']['details']
    for current, key in zip(pressure.type_slices, ('type_x', 'type_y')):
        probabilities = np.asarray(details[key])
        values[current] = np.log(probabilities[1:]/probabilities[0])
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/finite-fiber-restriction-catalogue.json'))
    parser.add_argument('--adaptive', type=Path, default=Path('evidence/restricted-closed-adaptive-screen.json'))
    parser.add_argument('--weighted', type=Path, default=Path('evidence/finite-fiber-restriction-search.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--trial-seconds', type=float, default=15)
    parser.add_argument('--candidates', type=int, default=24)
    parser.add_argument('--starts', type=int, default=3)
    parser.add_argument('--iterations', type=int, default=180)
    parser.add_argument('--seed', type=int, default=312729)
    args = parser.parse_args()
    source_bytes = {key: getattr(args, key).read_bytes() for key in ('catalogue', 'adaptive', 'weighted')}
    catalogue, adaptive, weighted = [json.loads(source_bytes[key]) for key in ('catalogue', 'adaptive', 'weighted')]
    assert catalogue['completed'] and adaptive['completed'] and weighted['completed']
    assert adaptive['input_sha256'] == sha256(source_bytes['catalogue']).hexdigest()
    selected = selected_cases(catalogue, adaptive, weighted, args.candidates)
    names = (Path(__file__).name, 'paired_mixed_place_pressure.py', 'f4_closed_place_oracle.py',
             'screen_restricted_closed_partitions.py', 'search_fused_mixed_weighted.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'scope': __doc__, 'arguments': vars(args) | {key: str(getattr(args, key)) for key in ('catalogue', 'adaptive', 'weighted', 'output')},
              'input_sha256': {key: sha256(value).hexdigest() for key, value in source_bytes.items()},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'selected_cases': [], 'records': [], 'winners': {}, 'completed_trials': 0,
              'attempted_trials': 0, 'completed': False, 'unresolved_objective_calls': 0}
    rng = np.random.default_rng(args.seed)
    oracles, old_rows = [], []
    for index, (ad, cat, old) in enumerate(selected):
        seed = old['seed']
        assert digest(seed['values']) == cat['tensor_sha256']
        oracle = F4ClosedPlaceOracle(seed)
        control = {'row_places': [0], 'column_places': [None], 'row_allocations': [1.], 'column_allocations': [1.],
                   'outer_log_scale': old['best']['parameters'][-1]}
        pressure = ScaledPairedMixedPressure(oracle, control)
        values = warm_start(pressure, control, old, rng, 0)
        pressure.evaluate(values)
        error = abs(pressure.best['rho']-old['best']['value'])
        assert error < 2e-8, (index, error)
        case = {key: cat[key] for key in ('source_index', 'subset_x', 'subset_y', 'degree', 'tensor_sha256')}
        case.update(seed=seed, adaptive_rho=ad['best']['rho'], old_weighted_rho=old['best']['value'],
                    two_place_control_error=error, partitions=partition_cases(ad, cat, old))
        report['selected_cases'].append(case)
        oracles.append(oracle)
        old_rows.append(old)
    tasks = [(i, j, start) for start in range(args.starts)
             for j in range(max(len(case['partitions']) for case in report['selected_cases']))
             for i, case in enumerate(report['selected_cases']) if j < len(case['partitions'])]
    began, stop = time.monotonic(), False
    deadline = began+args.seconds
    for case_index, partition_index, restart in tasks:
        if time.monotonic() >= deadline or stop:
            break
        case = report['selected_cases'][case_index]
        baseline = case['partitions'][partition_index]
        trial_start = time.monotonic()
        pressure = ScaledPairedMixedPressure(oracles[case_index], baseline, True,
                                             min(deadline, trial_start+args.trial_seconds))
        values = warm_start(pressure, baseline, old_rows[case_index], rng, restart)
        bounds = pressure.bounds()
        bounds[pressure.scale_index] = (1e-5, 4.)
        for index in range(pressure.type_slice.start, pressure.type_slice.stop):
            bounds[index] = (-20., 20.)
        for cut in pressure.height_slices:
            for current in cut:
                bounds[current] = [(-48., 48.)]*(current.stop-current.start)
        outcomes, completed, interruption, unresolved = [], True, None, 0
        def objective(value, smoothing):
            nonlocal unresolved
            try:
                return pressure.evaluate(value, smoothing)
            except ArithmeticError:
                unresolved += 1
                return 1e5+1e-6*float(value @ value), 2e-6*value
        try:
            objective(values, .025)
            for smoothing in (.025, .003, .0003):
                fit = minimize(lambda value: objective(value, smoothing), values, jac=True, method='L-BFGS-B',
                               bounds=bounds, options={'maxiter': args.iterations, 'maxfun': 1800,
                                                      'ftol': 2e-11, 'gtol': 3e-7, 'maxls': 30})
                outcomes.append({'smoothing': smoothing, 'success': bool(fit.success), 'iterations': int(fit.nit),
                                 'evaluations': int(fit.nfev), 'message': str(fit.message)})
                values = np.asarray(pressure.best['optimizer_parameters']) if pressure.best else fit.x
        except SearchDeadline:
            completed, interruption = False, 'Declared trial or overall wall-clock limit reached.'
        row = {'case_index': case_index, 'partition_index': partition_index, 'restart': restart,
               'tensor_sha256': case['tensor_sha256'], 'partition': baseline, 'completed': completed,
               'queries': pressure.queries, 'unresolved_objective_calls': unresolved,
               'seconds': time.monotonic()-trial_start, 'outcomes': outcomes, 'interruption': interruption,
               'best': {key: value for key, value in (pressure.best or {}).items() if key != 'place_records'}}
        report['records'].append(row)
        if pressure.best:
            key = str(case_index)
            if key not in report['winners'] or pressure.best['rho'] < report['winners'][key]['best']['rho']:
                report['winners'][key] = row | {'best': pressure.best}
            stop = pressure.best['rho'] < 2**-.5-1e-8
        report.update(completed_trials=sum(row['completed'] for row in report['records']),
                      attempted_trials=len(report['records']), elapsed_seconds=time.monotonic()-began,
                      unresolved_objective_calls=report['unresolved_objective_calls']+unresolved)
        atomic_write(args.output, report)
        print(json.dumps({'case': case_index, 'partition': partition_index, 'restart': restart,
                          'rho': pressure.best['rho'] if pressure.best else None,
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic()-began,
                  complete_trial_list=len(report['records']) == len(tasks), requested_trials=len(tasks),
                  stopping_reason='strict numerical threshold; exact certificate still required' if stop else
                  ('declared wall-clock cap' if time.monotonic() >= deadline else 'complete trial list'))
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
