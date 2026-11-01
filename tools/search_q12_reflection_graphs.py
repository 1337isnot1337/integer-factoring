#!/usr/bin/env python3
"""Timed numerical walks in the full 192-coordinate twelve-letter split chart."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import subprocess
import time

import numpy as np
from scipy.optimize import minimize

from mixed_graph_retractions import FixedProfileSupport, literal_digest
from normalized_support_gap import normalized_gap
from q12_reflection_graphs import HEIGHTS, PARAMETERS, Q, WORDS, canonical, natural, tensor
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure
from search_f4_recentered_crossings import root_for_profiles


class BudgetReached(Exception):
    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.output.exists() or args.output.with_suffix('.jsonl').exists():
        parser.error('Refusing to overwrite an earlier run.')
    root = Path(__file__).resolve().parents[1]
    for name, expected in manifest.get('source_sha256', {}).items():
        if sha256((root/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Registered source hash differs: {name}. Use its frozen source commit.')
    source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True):
        raise ValueError('Run from a clean frozen source.')
    begun = time.monotonic()
    report = {'status': 'running_incomplete', 'source_commit': source,
              'manifest': manifest, 'manifest_sha256': sha256(args.manifest.read_bytes()).hexdigest(),
              'tracks': [], 'best': None,
              'scope': 'Timed numerical actual-tensor search, not an exhaustive chart census, '
                       'an exact pressure bound, a global optimum or a stronger factoring score.'}
    objective_name = manifest.get('objective', 'support')
    if objective_name not in ('support', 'normalized_gap'):
        raise ValueError('Unknown search objective.')
    report['objective'] = objective_name

    def save():
        report['seconds'] = time.monotonic()-begun
        p = args.output.with_suffix('.tmp')
        p.write_text(json.dumps(report, indent=2)+'\n')
        p.replace(args.output)

    def improve(c, literal, best):
        if best is not None and (report['best'] is None or best['objective'] < report['best']['objective']):
            report['best'] = dict(best, coefficients=list(c), values_sha256=literal_digest(literal))
            out = {'best': report['best'], 'tensor': literal}
            p = args.output.with_suffix('.best-tensor.tmp')
            p.write_text(json.dumps(out, separators=(',', ':'))+'\n')
            p.replace(args.output.with_suffix('.best-tensor.json'))
            print(json.dumps({'event': 'best', 'rate': best['rate'], 'objective': best['objective'],
                              'seconds': time.monotonic()-begun}), flush=True)

    save()
    with args.output.with_suffix('.jsonl').open('x', buffering=1) as trace:
        def record(row):
            trace.write(json.dumps(row, separators=(',', ':'))+'\n')

        def fit(c, literal, initial, deadline, track_label, balanced=False):
            model = ScaledAsymmetricPressure(literal)
            best, last, queries, rejected = None, None, 0, []
            trial_end = min(deadline, time.monotonic()+manifest['trial_seconds'])
            started = time.monotonic()

            def objective(point):
                nonlocal best, last, queries
                if time.monotonic() >= trial_end:
                    raise BudgetReached
                queries += 1
                try:
                    value, gradient, details = model.evaluate(point)
                except (ArithmeticError, ValueError, OverflowError) as error:
                    rejected.append({'parameters': point.tolist(), 'error': str(error), 'query': queries})
                    raise
                objective_value, objective_gradient = ((value, gradient) if objective_name == 'support'
                    else normalized_gap(value, point, gradient))
                last = {'rate': value, 'objective': float(objective_value),
                        'parameters': point.tolist(), 'details': details}
                if best is None or objective_value < best['objective']:
                    best = last
                if value < report.get('minimum_support_query', {}).get('rate', math.inf):
                    report['minimum_support_query'] = dict(last, coefficients=list(c))
                return objective_value, objective_gradient

            opt = None
            try:
                bound = manifest['height_bound']
                initial = np.asarray(initial).copy()
                if balanced:
                    initial[-1] = 0.
                opt = minimize(objective, initial, jac=True, method='SLSQP',
                    bounds=[(-bound, bound)]*(model.parameters-2)+[(-9., math.log(4.)), (0., 0.) if balanced else (-2., 2.)],
                    options={'maxiter': manifest['maxiter'], 'ftol': 1e-10})
                status, message = 'optimizer_returned', str(opt.message)
            except BudgetReached:
                status, message = 'budget_reached_partial', 'trial/global time cap checked between queries'
            except (ArithmeticError, ValueError, OverflowError) as error:
                status, message = 'query_rejected_partial', str(error)
            improve(c, literal, best)
            row = {'track': track_label, 'coefficients': list(c), 'values_sha256': literal_digest(literal),
                   'initial_parameters': initial.tolist(), 'best': best, 'last': last,
                   'queries': queries, 'rejected_queries': rejected, 'status': status, 'message': message,
                   'optimizer_success': bool(opt.success) if opt is not None else False,
                   'balanced_angle_fixed': balanced,
                   'seconds': time.monotonic()-started}
            return best, row

        # One complete unit query supplies all initial scalar/angle evaluations.
        coefficients = (0,)*PARAMETERS
        base = tensor(coefficients)
        model = ScaledAsymmetricPressure(base)
        unit_orders = model.base.oracle.prefixes(natural(), 1)['orders']
        energy = [2*HEIGHTS[w] for w in WORDS[:Q]]
        fixed = FixedProfileSupport(np.asarray(unit_orders).reshape(2, Q*Q), energy, energy)
        report['coordinate_control'] = {'unit_orders': unit_orders, 'energy': energy,
            'balanced_fit': root_for_profiles(np.asarray(unit_orders).reshape(2, Q*Q), energy, energy)}
        start_deadline = time.monotonic()+manifest['initial_total_seconds']
        for angle in (.35, math.pi/4, math.pi/2-.35):
            def fixed_objective(point):
                value, gradient = fixed.evaluate(point)
                return ((value, gradient) if objective_name == 'support' else normalized_gap(value, point, gradient))
            fixed_fit = minimize(fixed_objective, (math.log(.07), math.log(math.tan(angle))), jac=True,
                method='SLSQP', bounds=((-9., math.log(4.)), (-2., 2.)),
                options={'maxiter': 120, 'ftol': 1e-11})
            initial = np.r_[natural()[:-1], fixed_fit.x]
            best, row = fit(coefficients, base, initial, start_deadline, 'coordinate-start')
            record(row)
        if report['best'] is None:
            raise ArithmeticError('No successful coordinate starting query.')
        frozen_start = report['best']
        for registered in manifest['tracks']:
            started = time.monotonic()
            deadline = started+registered['seconds']
            rng = random.Random(registered['seed'])
            current = frozen_start
            balanced = registered.get('balanced', False)
            if balanced:
                literal = tensor(current['coefficients'])
                initial = np.asarray(current['parameters'])
                best, row = fit(current['coefficients'], literal, initial, deadline,
                                registered['label']+'-fixed-angle-start', balanced=True)
                record(row)
                if best is None:
                    report['tracks'].append({'label': registered['label'], 'status': 'initial_query_rejected'})
                    save()
                    continue
                current = dict(best, coefficients=current['coefficients'], values_sha256=literal_digest(literal))
            state = {'label': registered['label'], 'status': 'running_incomplete', 'trials': 0,
                     'successful_queries': 0, 'rejected_queries': 0, 'accepted': 0, 'best': current}
            report['tracks'].append(state)
            while time.monotonic() < deadline:
                i = state['trials']
                if i and i % manifest['restart_period'] == 0:
                    current = state['best']
                c = list(current['coefficients'])
                changes = 1 if registered['mutation'] == 'single' else rng.choice((2, 4, 8))
                mutations = []
                for coordinate in rng.sample(range(PARAMETERS), changes):
                    increment = rng.randrange(1, 4)
                    c[coordinate] ^= increment
                    mutations.append([coordinate, increment])
                c = canonical(c)
                literal = tensor(c)
                initial = np.asarray(current['parameters'])
                best, row = fit(c, literal, initial, deadline, registered['label'], balanced=balanced)
                temperature = manifest['temperatures'][(i//manifest['temperature_block']) % len(manifest['temperatures'])]
                take = best is not None and (best['objective'] < current['objective'] or
                       rng.random() < math.exp(min(0., (current['objective']-best['objective'])/temperature)))
                row.update(mutations=mutations, accepted=take, temperature=temperature, index=i)
                record(row)
                if best is not None and best['objective'] < state['best']['objective']:
                    state['best'] = dict(best, coefficients=list(c), values_sha256=literal_digest(literal))
                if take:
                    current = dict(best, coefficients=list(c), values_sha256=literal_digest(literal))
                    state['accepted'] += 1
                state['trials'] += 1
                state['successful_queries'] += row['queries']-len(row['rejected_queries'])
                state['rejected_queries'] += len(row['rejected_queries'])
                state['seconds'] = time.monotonic()-started
                save()
            state.update(status='complete_timed_budget', last_state=current, seconds=time.monotonic()-started)
            save()
    report['status'] = 'complete_registered_timed_tracks'
    report['trace_sha256'] = sha256(args.output.with_suffix('.jsonl').read_bytes()).hexdigest()
    save()
    print(json.dumps({'status': report['status'], 'seconds': report['seconds'],
                      'best_support': report['best']['rate']}), flush=True)


if __name__ == '__main__':
    main()
