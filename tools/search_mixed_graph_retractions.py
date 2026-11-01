#!/usr/bin/env python3
"""Registered finite screen and bounded numerical walks on literal split graphs.

All pressure evaluations use complete profiles. Numerical best points are
search candidates only; an improving point requires a separate exact proof.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import random
import subprocess
import time

import numpy as np
from scipy.optimize import minimize

from mixed_graph_retractions import (FixedProfileSupport, canonical, literal_digest,
    natural, sparse_neighborhood, tensor)
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure
from search_f4_recentered_crossings import root_for_profiles


class BudgetReached(Exception):
    pass


def write_json(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    temporary.replace(path)


def source_record():
    root = Path(__file__).resolve().parents[1]
    names = ['mixed_graph_retractions.py', 'search_mixed_graph_retractions.py',
             'physical_split_crossings.py', 'reflection_dual_rail_gaussian.py',
             'scaled_asymmetric_crossing_pressure.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_f4_recentered_crossings.py']
    return {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
            'dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True)),
            'source_sha256': {n: sha256((root/'tools'/n).read_bytes()).hexdigest() for n in names}}


def fixed_fits(profiles):
    fixed = FixedProfileSupport(np.asarray(profiles).reshape(2, 25),
                                (0, 10, 6, 2, 8), (0, 10, 6, 2, 8))
    fits = []
    for angle in (.35, math.pi/4, math.pi/2-.35):
        start = np.asarray((math.log(.3), math.log(math.tan(angle))))
        result = minimize(fixed.evaluate, start, jac=True, method='SLSQP',
                          bounds=((-9., math.log(4.)), (-2., 2.)),
                          options={'maxiter': 120, 'ftol': 1e-11})
        value = fixed.evaluate(result.x, False)
        fits.append({'rate': value, 'variables': result.x.tolist(), 'success': bool(result.success),
                     'message': str(result.message), 'evaluations': int(result.nfev)})
    balanced = root_for_profiles(np.asarray(profiles).reshape(2, 25),
                                 (0, 10, 6, 2, 8), (0, 10, 6, 2, 8))
    if balanced['rho'] < 5/7-1e-8:
        raise ArithmeticError('The global natural-shape floor control failed.')
    return {'interior': min(fits, key=lambda f: f['rate']), 'fits': fits, 'balanced': balanced}


def free_fit(seed, start, args, rejected, deadline):
    support = ScaledAsymmetricPressure(seed)
    trial_deadline = min(deadline, time.monotonic()+args.trial_seconds)
    best, last, queries = None, None, 0
    started = time.monotonic()

    def objective(point):
        nonlocal best, last, queries
        if time.monotonic() >= trial_deadline:
            raise BudgetReached('trial or global wall-clock budget reached between queries')
        queries += 1
        try:
            value, derivative, details = support.evaluate(point)
        except (ArithmeticError, ValueError, OverflowError) as error:
            rejected.append({'parameters': point.tolist(), 'error': str(error), 'query': queries})
            raise
        last = {'rate': float(value), 'parameters': point.tolist(), 'details': details}
        if best is None or value < best['rate']:
            best = last
        return value, derivative

    result = None
    try:
        result = minimize(objective, start, jac=True, method='SLSQP',
                          bounds=[(-64., 64.)]*32 + [(-9., math.log(4.)), (-2., 2.)],
                          options={'maxiter': args.maxiter, 'ftol': 1e-10})
        status, message = 'optimizer_returned', str(result.message)
    except BudgetReached as error:
        status, message = 'budget_reached_partial', str(error)
    except (ArithmeticError, ValueError, OverflowError) as error:
        status, message = 'query_rejected_partial', str(error)
    return {'best': best, 'last': last, 'queries': queries, 'seconds': time.monotonic()-started,
            'status': status, 'message': message,
            'optimizer_success': bool(result.success) if result is not None else False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('screen', 'free', 'walk'), required=True)
    parser.add_argument('--output', type=Path, required=True, help='Summary JSON; .jsonl contains every trial.')
    parser.add_argument('--input', type=Path)
    parser.add_argument('--seconds', type=float, default=7200.)
    parser.add_argument('--trial-seconds', type=float, default=5.)
    parser.add_argument('--maxiter', type=int, default=100)
    parser.add_argument('--seed', type=int)
    args = parser.parse_args()
    if args.seconds <= 0 or args.trial_seconds <= 0 or args.maxiter < 1:
        parser.error('Budgets must be positive.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    trace_path = args.output.with_suffix('.jsonl')
    if args.output.exists() or trace_path.exists():
        parser.error('Refusing to overwrite a prior summary or trace.')
    started = time.monotonic()
    deadline = started + args.seconds
    seed_number = args.seed if args.seed is not None else int.from_bytes(os.urandom(16), 'big')
    rng = random.Random(seed_number)
    summary = {'status': 'running_incomplete', 'stage': args.stage, 'source': source_record(),
               'seed': seed_number, 'budget_seconds': args.seconds, 'trial_seconds': args.trial_seconds,
               'maxiter': args.maxiter, 'completed_trials': 0, 'successful_queries': 0,
               'rejected_queries': 0, 'best': None, 'trace': trace_path.name,
               'scope': 'Numerical finite search. Every returned profile is complete; unresolved queries '
                        'are rejected. No global optimum, exact strict support or score improvement.'}
    write_json(args.output, summary)
    with trace_path.open('x', buffering=1) as trace:
        def record(value):
            trace.write(json.dumps(value, separators=(',', ':'))+'\n')

        def update_best(coefficients, digest, best):
            if best is None or (summary['best'] is not None and best['rate'] >= summary['best']['rate']):
                return
            summary['best'] = dict(best, coefficients=list(coefficients), values_sha256=digest)
            literal = tensor(coefficients)
            write_json(args.output.with_suffix('.best-tensor.json'),
                       {'best': summary['best'], 'seed': literal})
            print(json.dumps({'event': 'best', 'stage': args.stage, 'trial': summary['completed_trials'],
                              'rate': best['rate'], 'seconds': time.monotonic()-started}), flush=True)

        if args.stage == 'screen':
            cache, selected, tensor_hashes = {}, [], set()
            representatives = sparse_neighborhood()
            summary['registered_representatives'] = len(representatives)
            for i, coefficients in enumerate(representatives):
                if time.monotonic() >= deadline:
                    break
                literal = tensor(coefficients)
                digest = literal_digest(literal)
                pressure = ScaledAsymmetricPressure(literal)
                rejected = []
                try:
                    profiles = pressure.base.oracle.prefixes(natural(), 1)['orders']
                    key = tuple(profiles)
                    if key not in cache:
                        cache[key] = fixed_fits(profiles)
                    fits = cache[key]
                    fit = fits['interior']
                    best = {'rate': fit['rate'], 'parameters': natural()[:-1].tolist()+fit['variables']}
                    update_best(coefficients, digest, best)
                    selected.append(dict(best, coefficients=coefficients, values_sha256=digest))
                    summary['successful_queries'] += 1
                    tensor_hashes.add(digest)
                    row = {'index': i, 'coefficients': coefficients, 'values_sha256': digest,
                           'unit_profiles': profiles, 'fits': fits, 'status': 'complete'}
                except (ArithmeticError, ValueError, OverflowError) as error:
                    summary['rejected_queries'] += 1
                    row = {'index': i, 'coefficients': coefficients, 'status': 'rejected',
                           'parameters': natural().tolist(), 'error': str(error)}
                record(row)
                summary['completed_trials'] += 1
                if i % 10 == 0:
                    summary['seconds'] = time.monotonic()-started
                    write_json(args.output, summary)
            chosen, seen = [], set()
            for point in sorted(selected, key=lambda p: p['rate']):
                if point['values_sha256'] not in seen and any(point['coefficients']):
                    chosen.append(point)
                    seen.add(point['values_sha256'])
                    if len(chosen) == 16:
                        break
            reference = next((p for p in selected if not any(p['coefficients'])), None)
            if reference is not None:
                chosen.insert(0, reference)
            summary.update(selected=chosen, distinct_literal_tensors=len(tensor_hashes),
                           distinct_unit_profile_pairs=len(cache),
                           status='complete' if summary['completed_trials'] == len(representatives)
                           else 'budget_reached_incomplete')
        elif args.stage == 'free':
            if args.input is None:
                parser.error('The free stage needs the completed screen summary.')
            previous = json.loads(args.input.read_text())
            if previous['stage'] != 'screen' or previous['status'] != 'complete':
                raise ValueError('Only a completed registered screen can select the free trials.')
            starts = previous['selected']
            summary['registered_trials'] = 3*len(starts)
            summary['input_sha256'] = sha256(args.input.read_bytes()).hexdigest()
            for item in starts:
                literal = tensor(item['coefficients'])
                for angle in (.35, math.pi/4, math.pi/2-.35):
                    if time.monotonic() >= deadline:
                        break
                    initial = np.asarray(item['parameters'])
                    initial[-1] = math.log(math.tan(angle))
                    rejected = []
                    fit = free_fit(literal, initial, args, rejected, deadline)
                    update_best(item['coefficients'], item['values_sha256'], fit['best'])
                    record({'index': summary['completed_trials'], 'coefficients': item['coefficients'],
                            'values_sha256': item['values_sha256'], 'initial_parameters': initial.tolist(),
                            'fit': fit, 'rejected_queries': rejected})
                    summary['completed_trials'] += 1
                    summary['successful_queries'] += fit['queries']-len(rejected)
                    summary['rejected_queries'] += len(rejected)
                    summary['seconds'] = time.monotonic()-started
                    write_json(args.output, summary)
            summary['status'] = ('complete' if summary['completed_trials'] == summary['registered_trials']
                                 else 'budget_reached_incomplete')
        else:
            if args.input is None:
                parser.error('The walk needs a source best point.')
            previous = json.loads(args.input.read_text())
            current = previous['best']
            if current is None:
                raise ValueError('The input has no successful point.')
            summary['input_sha256'] = sha256(args.input.read_bytes()).hexdigest()
            update_best(current['coefficients'], current['values_sha256'], current)
            accepted = 0
            while time.monotonic() < deadline:
                i = summary['completed_trials']
                if i and i % 200 == 0:
                    current = summary['best']
                coefficients = list(current['coefficients'])
                coordinate, increment = rng.randrange(20), rng.randrange(1, 4)
                coefficients[coordinate] ^= increment
                # The automorphism preserves profiles at every height, not just the initial gauge.
                coefficients = canonical(coefficients)
                literal = tensor(coefficients)
                digest = literal_digest(literal)
                initial = np.asarray(current['parameters'])
                rejected = []
                fit = free_fit(literal, initial, args, rejected, deadline)
                temperature = (.001, .005, .02, .1)[(i//50) % 4]
                best = fit['best']
                take = best is not None and (best['rate'] < current['rate'] or
                        rng.random() < math.exp(min(0., (current['rate']-best['rate'])/temperature)))
                update_best(coefficients, digest, best)
                record({'index': i, 'coefficients': coefficients, 'values_sha256': digest,
                        'mutation': [coordinate, increment], 'initial_parameters': initial.tolist(),
                        'fit': fit, 'temperature': temperature, 'accepted': take,
                        'rejected_queries': rejected})
                if take:
                    current = dict(best, coefficients=list(coefficients), values_sha256=digest)
                    accepted += 1
                summary['completed_trials'] += 1
                summary['successful_queries'] += fit['queries']-len(rejected)
                summary['rejected_queries'] += len(rejected)
                summary['accepted'] = accepted
                summary['seconds'] = time.monotonic()-started
                if i % 10 == 0:
                    write_json(args.output, summary)
            summary['status'] = 'complete_timed_budget'
            summary['last_state'] = current
    summary['seconds'] = time.monotonic()-started
    summary['trace_sha256'] = sha256(trace_path.read_bytes()).hexdigest()
    write_json(args.output, summary)
    print(json.dumps({k: summary[k] for k in ('status', 'stage', 'completed_trials',
                                            'successful_queries', 'rejected_queries', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
