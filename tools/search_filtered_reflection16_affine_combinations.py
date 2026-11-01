#!/usr/bin/env python3
"""Bounded exact affine sums across distinct filtered q16 frame directions.

Every component comes from a checked three-scalar affine span.  Different
directions and different pole orders may be combined.  Both complete
weighted cuts are recomputed; no infinity spectrum is assumed invariant.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar

from dense_coupled_reflection_frames import (
    rank_one_terms, transformed_from_terms, digest, save,
)
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor, physical_bit_energies
from reflection_laurent_physical_frames import (
    add, scale, degree, valuation, bell_normalize, coupled_admissible_poles,
)
from search_dense_coupled_reflection_frames import profile_rho
from search_filtered_reflection16_affine import aligned_directions, scaled_terms, affine_terms
from search_physical_code_gauges import start_params


def component_catalogue(base, poles):
    _, factors = bell_normalize(base)
    records, components, seen = [], [], set()
    for direction, (left, right) in enumerate(aligned_directions()):
        terms = rank_one_terms(base, left, right, factors)
        scalar_poles = [coupled_admissible_poles(scaled_terms(terms, scalar), base['degree'])
                        for scalar in (1, 2, 3)]
        allowed = set(scalar_poles[0]).intersection(*map(set, scalar_poles[1:]))
        for pole in poles:
            if pole not in allowed:
                continue
            for component in range(3):
                parameters = [int(i == component) for i in range(3)]
                changed = transformed_from_terms(base, affine_terms(terms, parameters), pole)
                delta = {index: add(first, second)
                         for index, (first, second) in enumerate(zip(base['values'], changed['values']))
                         if first != second}
                if not delta:
                    continue
                first = delta[min(delta)]
                least = ((first[0] | first[1]) & -(first[0] | first[1])).bit_length() - 1
                coefficient = ((first[0] >> least) & 1) | (((first[1] >> least) & 1) << 1)
                inverse = next(c for c in (1, 2, 3) if F4_PRODUCT[c][coefficient] == 1)
                delta = {index: scale(value, inverse) for index, value in delta.items()}
                key = sha256(json.dumps(sorted(delta.items()), separators=(',', ':')).encode()).hexdigest()
                if key in seen:
                    continue
                seen.add(key)
                assert max(degree(value) for value in delta.values()) < base['degree']
                assert all(not (a.bit_count() % 2 or b.bit_count() % 2) for a, b in delta.values())
                components.append(delta)
                records.append({'direction': direction, 'u': left, 'v': right, 'pole': pole,
                                'component': component, 'normalizing_scalar': inverse,
                                'nonzero_tensor_entries': len(delta), 'delta_sha256': key,
                                'all_three_scalar_frames_polynomial': True})
    return records, components


def combined_seed(base, components, word):
    values = list(base['values'])
    for index, coefficient in word:
        for position, value in components[index].items():
            values[position] = add(values[position], scale(value, coefficient))
    for index, (a, b) in enumerate(values):
        x, y, u, v = (index // 16 ** leg % 16 for leg in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('An affine sum changed a strict endpoint entry.')
    actual = max(degree(value) for value in values)
    if actual != base['degree']:
        raise ArithmeticError('An affine sum changed the common leading coefficient.')
    return dict(base, values=values, degree=actual,
                monomial_content=min(valuation(value) for value in values if value != (0, 0)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=2400)
    parser.add_argument('--trials', type=int, default=12000)
    parser.add_argument('--seed', type=int, default=312401)
    parser.add_argument('--poles', default='1,2,3,4,6,8,12,16')
    args = parser.parse_args()
    if args.seconds <= 0 or args.trials <= 0:
        parser.error('Use positive time and trial limits.')
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    base = tensor(8)
    weights = physical_bit_energies(8)
    energy = np.asarray([sum(value for bit, value in enumerate(weights) if state >> bit & 1)
                         for state in range(16)], float)
    natural = start_params(energy, 1)
    raw = RectangularOracle(base['values'], 16, 16).evaluate(
        natural, [1 / 16] * 16, [1 / 16] * 16, want_spectra=True)[2]
    baseline_orders = tuple(round(value) for value in raw)
    fit = minimize_scalar(lambda s: profile_rho(baseline_orders, energy, s), bounds=(.005, .4),
                          method='bounded', options={'xatol': 1e-9})
    names = (Path(__file__).name, 'search_filtered_reflection16_affine.py',
             'dense_coupled_reflection_frames.py', 'reflection_laurent_physical_frames.py',
             'reflection_dual_rail_gaussian.py', 'search_dense_coupled_reflection_frames.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_physical_code_gauges.py')
    report = {'scope': 'bounded cross-direction affine sums, exact complete tensors and full natural profiles; no optimality certificate',
              'started_utc': started_utc,
              'requested_seconds': args.seconds, 'requested_trials': args.trials, 'seed': args.seed,
              'requested_poles': [int(x) for x in args.poles.split(',')],
              'baseline_rho': float(fit.fun), 'baseline_scale': float(fit.x),
              'profiles': [{'orders': list(baseline_orders), 'rho': float(fit.fun),
                            'scale': float(fit.x), 'is_original_baseline': True}],
              'records': [], 'counters': {'proposed': 0, 'duplicate_tensors': 0,
                                        'unresolved_profiles': 0},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    records, components = component_catalogue(base, report['requested_poles'])
    report['components'] = records
    report['preparation_seconds'] = time.monotonic() - started
    profiles = {baseline_orders: 0}
    seen = {digest(base['values'])}
    rng = np.random.default_rng(args.seed)
    pool = []
    save(args.output, report)
    for trial in range(args.trials):
        if time.monotonic() - started >= args.seconds:
            break
        report['counters']['proposed'] += 1
        if pool and trial % 5 < 2:
            parent = pool[int(rng.integers(len(pool)))]
            word = dict(parent['word'])
            changes = 1 if trial % 2 else 2
            for _ in range(changes):
                index = int(rng.integers(len(components)))
                word[index] = word.get(index, 0) ^ int(rng.integers(1, 4))
            word = sorted((index, coefficient) for index, coefficient in word.items() if coefficient)
        else:
            length = 2 if trial % 4 else int(rng.integers(3, 6))
            word = sorted((int(index), int(rng.integers(1, 4)))
                          for index in rng.choice(len(components), length, replace=False))
        changed = combined_seed(base, components, word)
        exact_hash = digest(changed['values'])
        if exact_hash in seen:
            report['counters']['duplicate_tensors'] += 1
            continue
        seen.add(exact_hash)
        record = {'trial': trial, 'word': word, 'seed_sha256': exact_hash,
                  'degree': changed['degree'], 'monomial_content': changed['monomial_content']}
        try:
            raw = RectangularOracle(changed['values'], 16, 16).evaluate(
                natural, [1 / 16] * 16, [1 / 16] * 16, want_spectra=True)[2]
            key = tuple(round(value) for value in raw)
            assert max(abs(a - b) for a, b in zip(raw, key)) < 1e-8
            if key not in profiles:
                fit = minimize_scalar(lambda s: profile_rho(key, energy, s), bounds=(.005, .4),
                                      method='bounded', options={'xatol': 2e-7, 'maxiter': 30})
                profiles[key] = len(report['profiles'])
                report['profiles'].append({'orders': list(key), 'rho': float(fit.fun),
                                           'scale': float(fit.x)})
            record['profile_index'] = profiles[key]
            record['natural_rho'] = report['profiles'][record['profile_index']]['rho']
            pool.append(record)
            pool.sort(key=lambda row: row['natural_rho'])
            retained, used = [], set()
            for row in pool:
                if row['profile_index'] not in used:
                    retained.append(row)
                    used.add(row['profile_index'])
                if len(retained) == 24:
                    break
            pool = retained
        except (ArithmeticError, ValueError) as error:
            report['counters']['unresolved_profiles'] += 1
            record['status'] = 'No spectrum or pressure assigned: ' + str(error)
        report['records'].append(record)
        if len(report['records']) % 20 == 0:
            report.update(elapsed_seconds=time.monotonic() - started,
                          best_natural_rho=min((row['natural_rho'] for row in report['records']
                                                if 'natural_rho' in row), default=None),
                          completed=False)
            save(args.output, report)
            print(json.dumps({'records': len(report['records']), 'profiles': len(profiles),
                              'best_natural_rho': report['best_natural_rho'],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        if record.get('natural_rho', math.inf) < 1 / math.sqrt(2) - 1e-6:
            report['stopping_reason'] = 'Numerical threshold candidate requires exact independent replay.'
            save(args.output, report)
            break
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started,
                  best_natural_rho=min((row['natural_rho'] for row in report['records']
                                        if 'natural_rho' in row), default=None),
                  stopping_reason=report.get('stopping_reason', 'Declared time or trial limit reached.'))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in
                      ('completed', 'elapsed_seconds', 'best_natural_rho', 'stopping_reason')}))


if __name__ == '__main__':
    main()
