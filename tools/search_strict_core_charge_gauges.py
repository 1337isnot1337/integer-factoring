#!/usr/bin/env python3
"""Screen every two-mode charge class with arbitrary physical core gauges."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from strict_core_charge_profiles import affine_charge_representatives
from strict_core_charge_gauges import CoreChargePressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=16)
    parser.add_argument('--initial', type=Path)
    parser.add_argument('--maxfun', type=int, default=180)
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard < args.shards:
        raise ValueError('The shard index must lie in the selected shard count.')
    initial = np.zeros(27)
    initial[-1] = 0.25
    if args.initial:
        initial = np.asarray(json.loads(args.initial.read_text())['best']['parameters'], dtype=float)
        if initial.shape != (27,):
            raise ValueError('Expected a two-mode physical-core parameter checkpoint.')
    records = []
    started = time.monotonic()

    def save():
        args.output.write_text(json.dumps({'status': 'exploratory full core-gauge search; no exact rate claim',
            'interval': args.interval, 'shard': args.shard, 'shards': args.shards,
            'maxfun': args.maxfun, 'seconds': time.monotonic() - started,
            'records': records}, indent=2) + '\n')

    for index, charges in enumerate(affine_charge_representatives(2)):
        if index % args.shards != args.shard:
            continue
        # A constant charge shift leaves the projection identical. Choosing
        # a most frequent value as zero makes the exact phase program shorter.
        common = max(range(3), key=lambda value: (charges.count(value), -value))
        effective = tuple((value - common) % 3 for value in charges)
        case_started = time.monotonic()
        model = CoreChargePressure(args.interval, effective)
        record = {'class_index': index, 'charges': list(charges),
                  'equivalent_program_charges': list(effective), 'best': None, 'runs': []}
        records.append(record)

        def objective(parameters):
            offsets = np.zeros((2, 4, 4))
            offsets[:, :, 1:] = parameters[:-3].reshape((2, 4, 3))
            try:
                value, (do, db, ds), detail = model.evaluate(parameters[-1], offsets,
                                                            parameters[-3:-1], True)
            except ArithmeticError:
                return 100.0, np.zeros(27)
            if record['best'] is None or value < record['best']['rho'] - 1e-10:
                record['best'] = {**detail, 'parameters': parameters.tolist()}
                record['queries'] = model.queries
                record['seconds'] = time.monotonic() - case_started
                save()
            return value, np.r_[do[:, :, 1:].reshape(-1), db, ds]

        starts = [initial.copy(), np.r_[np.zeros(26), initial[-1]]]
        for seed in starts:
            run = minimize(objective, seed, jac=True, method='L-BFGS-B',
                           bounds=[(-32, 32)] * 24 + [(-4 * args.interval, 4 * args.interval)] * 2
                           + [(0.002, 3.0)],
                           options={'maxfun': args.maxfun, 'maxiter': args.maxfun,
                                    'maxls': 20, 'ftol': 1e-12, 'gtol': 1e-8})
            record['runs'].append({'rho': float(run.fun), 'success': bool(run.success),
                                   'message': str(run.message), 'evaluations': int(run.nfev)})
        record['queries'] = model.queries
        record['seconds'] = time.monotonic() - case_started
        save()
        print(json.dumps({'class_index': index, 'charges': charges,
                          'rho': record['best']['rho'], 'queries': model.queries,
                          'seconds': record['seconds']}), flush=True)


if __name__ == '__main__':
    main()
