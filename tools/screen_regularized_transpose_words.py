#!/usr/bin/env python3
"""Registered literal transpose-word pilot, complete profiles and bounded fits."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from random import Random
import subprocess
import time

import numpy as np
from scipy.optimize import minimize

from bounded_complete_support_fit import fit
from mixed_graph_retractions import FixedProfileSupport
from native_rectangular_smith_f4 import RectangularOracle
from normalized_support_gap import normalized_gap
from regularized_transpose_words import construct
from search_f4_recentered_crossings import root_for_profiles


def gates(rng, q, count, shears):
    result = []
    for _ in range(count):
        gate = []
        for _ in range(shears):
            i, j = rng.sample(range(q), 2)
            gate.append([i, j, rng.randrange(1, 4)])
        result.append(gate)
    return result


def pair(rng, q, power, shears):
    a = gates(rng, q, 1, shears)[0]
    i, j = rng.sample(range(q), 2)
    return [{'shears': a, 'power': power},
            {'shears': a+[[i, j, rng.randrange(1, 4)]], 'power': power}]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refusing to overwrite an earlier run.')
    manifest = json.loads(args.manifest.read_text())
    root = Path(__file__).resolve().parents[1]
    for name, expected in manifest['source_sha256'].items():
        if sha256((root/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Registered source changed: {name}.')
    source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True):
        raise ValueError('Run from a clean frozen source.')
    begun = time.monotonic()
    report = {'status': 'running_incomplete', 'source_commit': source,
              'manifest': manifest, 'cases': [],
              'scope': 'Finite exploratory complete-profile pilot with bounded local fits. '
                       'No global minimum, exact support certificate or stronger factoring score.'}

    def save():
        report['seconds'] = time.monotonic()-begun
        temp = args.output.with_suffix('.tmp')
        temp.write_text(json.dumps(report, indent=2)+'\n')
        temp.replace(args.output)

    save()
    def remember(literal, metadata, candidate):
        if candidate is None or candidate['rate'] >= report.get('best_support', {}).get('rate', math.inf):
            return
        report['best_support'] = dict(candidate, values_sha256=metadata['values_sha256'], metadata=metadata)
        path = args.output.with_suffix('.best-tensor.tmp')
        path.write_text(json.dumps({'best': report['best_support'], 'tensor': literal}, separators=(',', ':'))+'\n')
        path.replace(args.output.with_suffix('.best-tensor.json'))

    for case in manifest['cases']:
        if time.monotonic()-begun > manifest['total_seconds']:
            break
        rng = Random(case['seed'])
        q, kind = case['q'], case['kind']
        if kind == 'diagonal':
            pre, post = [], []
        elif kind == 'one-pair':
            pre, post = [], pair(rng, q, case['append_power'], 48)
        elif kind == 'two-pairs':
            pre, post = pair(rng, q, case['prepend_power'], 48), pair(rng, q, case['append_power'], 48)
        elif kind == 'generic':
            pre = gates(rng, q, case['gates_per_side'], 48)
            post = gates(rng, q, case['gates_per_side'], 48)
        elif kind == 'fixture':
            pre, post = case['prepend'], case['append']
        else:
            raise ValueError('Unknown registered word kind.')
        literal, metadata = construct(q, pre, post, tuple(case['base_powers']))
        if 'expected_values_sha256' in case:
            assert metadata['values_sha256'] == case['expected_values_sha256']
        row = {'label': case['label'], 'metadata': metadata, 'status': 'running_incomplete', 'fits': []}
        report['cases'].append(row)
        params = [0.]*(8*(q-1))+[1.]
        try:
            profiles = RectangularOracle(literal['values'], q, q).prefixes(params, 1)['orders']
        except (ArithmeticError, ValueError) as error:
            row.update(status='unit_query_rejected', rejected_parameters=params, error=str(error))
            save()
            continue
        row['unit_profiles'] = profiles
        raw = np.asarray(profiles).reshape(2, q*q)
        row['bare_balanced'] = root_for_profiles(raw, [0]*q, [0]*q)
        fixed = FixedProfileSupport(raw, [0]*q, [0]*q)
        closed = []
        for angle in (.35, math.pi/4, math.pi/2-.35):
            def objective(point):
                value, gradient = fixed.evaluate(point)
                return normalized_gap(value, point, gradient)
            result = minimize(objective, (math.log(.1), math.log(math.tan(angle))), jac=True,
                method='SLSQP', bounds=((-9., math.log(4.)), (-2., 2.)),
                options={'maxiter': 120, 'ftol': 1e-11})
            closed.append({'variables': result.x.tolist(), 'objective': float(result.fun),
                           'rate': fixed.evaluate(result.x, False), 'success': bool(result.success)})
        row['bare_asymmetric'] = closed
        best = min(closed, key=lambda item: item['objective'])
        for value in closed:
            remember(literal, metadata, {'rate': value['rate'], 'parameters':
                     [0.]*(8*(q-1))+value['variables'], 'from': 'fixed-unit-profiles'})
        balanced_scale = row['bare_balanced']['scale']/math.sqrt(2)
        for balanced in (True, False):
            if time.monotonic()-begun >= manifest['total_seconds']:
                break
            point = np.r_[np.zeros(8*(q-1)),
                math.log(balanced_scale) if balanced else best['variables'][0],
                0. if balanced else best['variables'][1]]
            # Fixed recorded jitter avoids starting every free fit on a large tie face.
            local = np.random.default_rng(rng.getrandbits(128))
            point[:-2] += local.normal(0, .05, len(point)-2)
            record = fit(literal, point, seconds=min(manifest['trial_seconds'],
                manifest['total_seconds']-(time.monotonic()-begun)),
                maxiter=manifest['maxiter'], balanced=balanced, height_bound=128.)
            row['fits'].append(record)
            remember(literal, metadata, record['best'])
            save()
        row['status'] = 'complete_attempts' if len(row['fits']) == 2 else 'budget_reached_incomplete'
        save()
        print(json.dumps({'label': case['label'], 'degree': metadata['degree'],
                          'leading_vector_ranks': metadata['zero_R_vector_matrix_ranks'],
                          'bare_rho': row['bare_balanced']['rho'],
                          'best_fit_supports': [r['best']['rate'] if r['best'] else None for r in row['fits']]}), flush=True)
    report['status'] = ('complete_registered_attempts' if len(report['cases']) == len(manifest['cases'])
                        and all(c['status'] in ('complete_attempts', 'unit_query_rejected') for c in report['cases'])
                        else 'budget_reached_incomplete')
    save()
    print(json.dumps({'status': report['status'], 'seconds': report['seconds']}), flush=True)


if __name__ == '__main__':
    main()
