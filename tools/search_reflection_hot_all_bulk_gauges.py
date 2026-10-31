#!/usr/bin/env python3
"""Bounded physical-height searches in all sixteen proved bulk-order cones.

This driver uses the frozen complete-profile evaluator. It explores
finite coefficient rows and numerical boxes; it is not an obstruction
theorem or a pressure certificate.
"""
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
from reflection_hot_leg_bulk_gauges import BOUNDARIES
from reflection_hot_core_profiles import profile_data, optimize_scale
from search_reflection_hot_core_filters import factor_from_row
from search_reflection_hot_state_gauges import EXCEPTIONAL


def bulk_and_offsets(low, excess, slack, larger):
    """Choose the larger physical leg independently at all four boundaries."""
    bulk = np.zeros((2, 4))
    offsets = np.zeros((2, 4, 2))
    for place, sides in enumerate(BOUNDARIES):
        for side, legs in enumerate(sides):
            big = legs[larger[2*place+side]]
            bulk[place, list(legs)] = low[place, side]
            bulk[place, big] += excess[place, side]
            for leg in legs:
                offsets[place, leg, 1] = slack[place, leg]-bulk[place, leg]-low[place, side]
    return bulk, offsets


def parameterized_query(model, parameters, larger):
    ell = math.exp(parameters[0])
    low = parameters[1:5].reshape(2, 2)
    excess = parameters[5:9].reshape(2, 2)
    slack = parameters[9:].reshape(2, 4)
    bulk, offsets = bulk_and_offsets(low, excess, slack, larger)
    value, (dh, db, ds), detail = model.evaluate(ell, offsets, bulk, True)
    dl, de = np.zeros((2, 2)), np.zeros((2, 2))
    for place, sides in enumerate(BOUNDARIES):
        for side, legs in enumerate(sides):
            big = legs[larger[2*place+side]]
            dl[place, side] = db[place, list(legs)].sum()-2*dh[place, list(legs), 1].sum()
            de[place, side] = db[place, big]-dh[place, big, 1]
    gradient = np.r_[ell*ds, dl.reshape(-1), de.reshape(-1), dh[:, :, 1].reshape(-1)]
    point = {'rho': value, 'log_s': ell, 'low_bulks': low.tolist(), 'bulk_excesses': excess.tolist(),
             'slacks': slack.tolist(), 'physical_leg_bulks': bulk.tolist(),
             'relative_hot_heights': offsets[:, :, 1].tolist(),
             'total_added_hot_heights': (bulk+offsets[:, :, 1]).tolist(),
             'active_orders': detail['active_orders'], 'spectator_orders': detail['spectator_orders'],
             'type_x_core': detail['type_x_core'], 'type_y_core': detail['type_y_core'],
             'cold_x_energies': detail['cold_x_energies'], 'cold_y_energies': detail['cold_y_energies'],
             'actual_cut_branches': [model.records[place][0]['data']['bulk_branch'] for place in range(2)]}
    return value, gradient, point


def derivative_controls():
    """Finite differences check the new 17-parameter chain rule in all cones."""
    model = AllBulkHotStatePressure(4, factor_from_row((0, 1, 1, 1, 1, 0, 0, 1), 1))
    point = np.r_[math.log(.317), [.73, .91, 1.17, 1.31], [.83, 1.09, 1.27, 1.43],
                  [.23, .37, .41, .59, .67, .71, .89, 1.03]]
    errors = []
    epsilon = 1e-6
    for larger in product(range(2), repeat=4):
        _, gradient, _ = parameterized_query(model, point, larger)
        for i in range(17):
            delta = np.zeros(17)
            delta[i] = epsilon
            forward = parameterized_query(model, point+delta, larger)[0]
            backward = parameterized_query(model, point-delta, larger)[0]
            errors.append(abs((forward-backward)/(2*epsilon)-gradient[i]))
    if max(errors) >= 5e-8:
        raise ArithmeticError('The bulk/slack chain-rule control failed.')
    return {'comparisons': len(errors), 'maximum_absolute_error': max(errors)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=64)
    parser.add_argument('--random', type=int, default=4)
    parser.add_argument('--seed', type=int, default=133749)
    parser.add_argument('--maxfun', type=int, default=70)
    parser.add_argument('--excess-floor', type=float, default=1e-6)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.excess_floor < 1:
        raise ValueError('Use a small positive excess floor to remain inside each fixed cone.')
    start = time.monotonic()
    rng = np.random.default_rng(args.seed)
    rows = [(0,)*8, (1, 1, 0, 0, 1, 1, 1, 1)]+list(EXCEPTIONAL)
    while len(rows) < 9+args.random:
        row = tuple(map(int, rng.integers(0, 4, size=8)))
        if row not in rows:
            rows.append(row)
    baseline = optimize_scale(profile_data(args.interval, 1, {0: {0: 1}}))
    result = {'status': 'Exact full source profiles with bounded numerical optimization across all sixteen bulk-order cones.',
              'interval': args.interval, 'core': 1, 'ramification': 1, 'seed': args.seed,
              'planned_candidates': len(rows), 'planned_cones': 16*len(rows), 'complete': False,
              'gaussian_natural': baseline, 'low_bulk_bounds': [0, 8],
              'excess_bounds': [args.excess_floor, 16], 'slack_bounds': [0, 24],
              'scale_bounds': [baseline['log_s']/4, baseline['log_s']*3],
              'starts_per_cone': 2, 'nominal_maxfun': args.maxfun,
              'chain_rule_control': derivative_controls(), 'records': []}
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/search_reflection_hot_all_bulk_gauges.py',
               'tools/reflection_hot_all_bulk_gauges.py', 'tools/reflection_hot_leg_bulk_gauges.py']

    def save():
        result['seconds'] = time.monotonic()-start
        result['completed_candidates'] = sum(len(row['cones']) == 16 for row in result['records'])
        result['completed_cones'] = sum(len(row['cones']) for row in result['records'])
        result['source_sha256'] = {name: sha256((root/name).read_bytes()).hexdigest() for name in sources}
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        temporary.replace(args.output)

    bounds = [(math.log(baseline['log_s']/4), math.log(baseline['log_s']*3))]
    bounds += [(0., 8.)]*4+[(args.excess_floor, 16.)]*4+[(0., 24.)]*8
    for row in rows:
        factor = factor_from_row(row, 1)
        model = AllBulkHotStatePressure(args.interval, factor, profile_precision=max(512, 8*args.interval))
        natural = np.zeros(17)
        natural[0] = math.log(baseline['log_s'])
        natural_point = parameterized_query(model, natural, (0, 0, 0, 0))[2]
        record = {'row': row, 'factor': factor, 'natural': natural_point,
                  'best': natural_point, 'cones': []}
        result['records'].append(record)
        for larger in product(range(2), repeat=4):
            best, runs = None, []

            def objective(parameters):
                nonlocal best
                value, gradient, point = parameterized_query(model, parameters, larger)
                if best is None or value < best['rho']:
                    best = point
                return value, gradient

            for low, excess in ((0., args.excess_floor), (1., 1.)):
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
                          'completed_cones': result['completed_cones'],
                          'rho': record['best']['rho'], 'queries': model.queries,
                          'seconds': result['seconds']}), flush=True)
    result['complete'] = True
    result['best'] = min(result['records'], key=lambda row: row['best']['rho'])
    save()
    print(json.dumps({'complete': True, 'candidates': len(rows),
                      'best_rho': result['best']['best']['rho'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    main()
