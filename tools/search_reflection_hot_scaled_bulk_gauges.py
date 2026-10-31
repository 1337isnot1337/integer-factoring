#!/usr/bin/env python3
"""Explore bulk heights comparable with the interval on selected hot defects."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from reflection_hot_all_bulk_gauges import AllBulkHotStatePressure
from reflection_hot_core_profiles import profile_data, optimize_scale
from search_reflection_hot_all_bulk_gauges import parameterized_query
from search_reflection_hot_core_filters import factor_from_row
from search_reflection_hot_state_gauges import EXCEPTIONAL


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=128)
    parser.add_argument('--maxfun', type=int, default=100)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    M, floor = args.interval, 1e-6
    rows = [(0,)*8]+list(EXCEPTIONAL[:2])
    baseline = optimize_scale(profile_data(M, 1, {0: {0: 1}}))
    result = {'status': 'Exact full profiles with exploratory interval-scaled bulk-height boxes.',
              'interval': M, 'core': 1, 'ramification': 1, 'complete': False,
              'planned_candidates': len(rows), 'planned_cones': 16*len(rows),
              'gaussian_natural': baseline, 'low_bulk_bounds': [0, M/2],
              'excess_bounds': [floor, 2*M], 'slack_bounds': [0, 2*M],
              'scale_bounds': [baseline['log_s']/4, baseline['log_s']*3],
              'starts_low_excess': [[0, floor], [M/16, M/8], [M/4, M/2]],
              'nominal_maxfun': args.maxfun,
              'optimizer_parameters_projected_to_stated_box_before_evaluation': True,
              'records': []}
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/search_reflection_hot_scaled_bulk_gauges.py',
               'tools/search_reflection_hot_all_bulk_gauges.py',
               'tools/reflection_hot_all_bulk_gauges.py']
    hashes = {name: sha256((root/name).read_bytes()).hexdigest() for name in sources}

    def save():
        result['seconds'] = time.monotonic()-start
        result['completed_candidates'] = sum(len(row['cones']) == 16 for row in result['records'])
        result['completed_cones'] = sum(len(row['cones']) for row in result['records'])
        result['source_sha256'] = hashes
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        temporary.replace(args.output)

    bounds = [(math.log(baseline['log_s']/4), math.log(baseline['log_s']*3))]
    bounds += [(0., M/2)]*4+[(floor, 2*M)]*4+[(0., 2*M)]*8
    lower_bounds, upper_bounds = np.array(bounds).T
    for row in rows:
        factor = factor_from_row(row, 1)
        model = AllBulkHotStatePressure(M, factor, profile_precision=max(512, 8*M))
        natural = np.zeros(17)
        natural[0] = math.log(baseline['log_s'])
        natural_point = parameterized_query(model, natural, (0, 0, 0, 0))[2]
        record = {'row': row, 'factor': factor, 'natural': natural_point,
                  'best': natural_point, 'cones': [],
                  'box_projection_queries': 0, 'maximum_box_projection': 0.}
        result['records'].append(record)
        for larger in product(range(2), repeat=4):
            best, runs = None, []

            def objective(parameters):
                nonlocal best
                # L-BFGS-B can return a negative slack of one floating ulp.
                # Project to its stated box; never weaken the regularity gate.
                projected = np.clip(parameters, lower_bounds, upper_bounds)
                difference = float(np.max(np.abs(projected-parameters)))
                if difference:
                    record['box_projection_queries'] += 1
                    record['maximum_box_projection'] = max(record['maximum_box_projection'], difference)
                value, gradient, point = parameterized_query(model, projected, larger)
                if best is None or value < best['rho']:
                    best = point
                return value, gradient

            for low, excess in result['starts_low_excess']:
                initial = natural.copy()
                initial[1:5] = low
                initial[5:9] = excess
                run = minimize(objective, initial, jac=True, method='L-BFGS-B', bounds=bounds,
                               options={'maxfun': args.maxfun, 'maxiter': args.maxfun,
                                        'maxls': 20, 'ftol': 1e-13, 'gtol': 2e-9})
                runs.append({'rho': float(run.fun), 'evaluations': int(run.nfev),
                             'success': bool(run.success), 'message': str(run.message)})
            record['cones'].append({'larger_leg_choices': larger, 'best': best, 'starts': runs})
            if best['rho'] < record['best']['rho']:
                record['best'] = best
            record['queries'] = model.queries
            record['source_branches'] = [
                {'place': key[0], 'pivot': key[1], 'swapped': key[2],
                 'source_metadata': value[0][0]['metadata'], 'sector_metadata': value[0][1]}
                for key, value in sorted(model.branch_cache.items())]
            save()
        print(json.dumps({'completed_candidates': result['completed_candidates'],
                          'completed_cones': result['completed_cones'], 'rho': record['best']['rho'],
                          'queries': model.queries, 'seconds': result['seconds']}), flush=True)
    result['complete'] = True
    result['best'] = min(result['records'], key=lambda row: row['best']['rho'])
    save()
    print(json.dumps({'complete': True, 'candidates': len(rows),
                      'best_rho': result['best']['best']['rho'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    main()
