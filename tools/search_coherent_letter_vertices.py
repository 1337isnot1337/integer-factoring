#!/usr/bin/env python3
"""Sample complete profiles at vertices of the shared letter-gauge cone.

Linear objectives use positive, randomly tilted physical weights. Every
selected cone vertex gets both complete profiles, followed by a separate
scale grid and local refinements. This is a bounded exploratory sample,
not a proof of global optimality or an exact pressure witness.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import linprog

from coherent_weighted_gauges import SharedGaugeOracle
from search_coherent_core_products import pressure_grid
from search_coherent_hot_retractions import _atomic_write, q32_catalogue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials', type=int, default=60)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--seed', type=int, default=104729)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    description, embedding, retraction, modes = list(q32_catalogue())[55]
    names = (Path(__file__).name, 'coherent_weighted_gauges.py', 'coherent_hot_mixed.py',
             'coherent_hot_mixed.cpp', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds_fast.py', 'coherent_hot_compounds.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'search_coherent_core_products.py')
    report = {'scope': __doc__, 'completed': False, 'description': description,
              'embedding': embedding, 'retraction': retraction, 'modes': modes,
              'seed': args.seed, 'trial_budget': args.trials, 'seconds_budget': args.seconds,
              'records': [], 'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'exact_pressure_certificate': False}
    _atomic_write(args.output, report)
    started = time.monotonic()
    oracle = SharedGaugeOracle(embedding, retraction, modes)
    random = np.random.default_rng(args.seed)
    report.update(deletion_edges=oracle.edges, cone=oracle.cone.astype(int).tolist(),
                  free_labels=oracle.free)
    bounds = [(-float(oracle.heights[j]) / 2., 32.) for j in oracle.free]
    report['normalized_potential_bounds'] = bounds
    seen = set()
    best = float('inf')
    attempts = 0
    while len(report['records']) < args.trials and time.monotonic() - started < args.seconds:
        index = len(report['records'])
        if attempts == 0:
            theta = np.zeros(oracle.q)
            lp = None
        else:
            tilt = random.choice((.3, .7, 1.5, 3.))
            scale = random.uniform(.1, .8)
            weights = np.exp(-scale * oracle.heights[oracle.free]
                             + random.normal(0., tilt, oracle.q - 1))
            fit = linprog(weights, A_ub=-oracle.cone[:, 1:], b_ub=oracle.cone[:, 0],
                          bounds=bounds, method='highs')
            if not fit.success:
                raise ArithmeticError(f'Cone linear program failed: {fit.message}')
            theta = np.zeros(oracle.q)
            theta[oracle.free] = fit.x
            lp = {'objective': weights.tolist(), 'value': float(fit.fun),
                  'message': str(fit.message)}
        attempts += 1
        rational = [Fraction(float(value)).limit_denominator(4096) for value in theta]
        key = tuple((value.numerator, value.denominator) for value in rational)
        if key in seen:
            if attempts > 20 * args.trials:
                break
            continue
        seen.add(key)
        theta = np.asarray([float(value) for value in rational])
        parameters = np.asarray([1.] + theta[oracle.free].tolist())
        slacks = oracle.cone @ parameters
        if slacks.min() < -1e-8:
            raise ArithmeticError('Rational reconstruction left the sufficient cone.')
        report['in_progress'] = {'index': index, 'potentials': theta.tolist(),
                                 'rational_potentials': key, 'linear_program': lp}
        _atomic_write(args.output, report)
        cuts = [oracle._cut(cut, 1., theta) for cut in (0, 1)]
        profile = {'zero_orders_t': sorted(cuts[0]['orders']),
                   'infinity_orders_t': sorted(cuts[1]['orders']),
                   'heights': (oracle.heights + 2 * theta).tolist()}
        pressure = pressure_grid(profile)
        record = {'index': index, 'potentials': theta.tolist(), 'rational_potentials': key,
                  'linear_program': lp, 'minimum_cone_slack': float(slacks.min()),
                  'profile': profile, 'pressure': pressure,
                  'elapsed_seconds': time.monotonic() - started}
        report['records'].append(record)
        report['in_progress'] = None
        report['attempts'] = attempts
        report['elapsed_seconds'] = record['elapsed_seconds']
        if pressure['best']['rho'] < best:
            best = pressure['best']['rho']
            report['best_index'] = index
            print(json.dumps({'index': index, 'rho': best, 'scale_t': pressure['best']['scale_t'],
                              'seconds': report['elapsed_seconds']}), flush=True)
        _atomic_write(args.output, report)
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started,
                  requested_sample_complete=len(report['records']) == args.trials,
                  global_optimality_claimed=False, factoring_exponent_improved=False)
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
