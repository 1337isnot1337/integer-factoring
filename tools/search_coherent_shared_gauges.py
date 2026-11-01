#!/usr/bin/env python3
"""Exploratory complete-profile search of shared coherent letter gauges.

The recorded cone is a sufficient condition for both cut conjugations to
remain regular. Numerical pivot weights and pressure roots need exact
rational replay before any asymptotic gain is claimed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import LinearConstraint, minimize

from coherent_weighted_gauges import SharedGaugeOracle
from search_coherent_hot_retractions import _atomic_write, q32_catalogue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coordinate', action='store_true')
    parser.add_argument('--maxiter', type=int, default=100)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    index = 56 if args.coordinate else 55
    description, p, l, modes = list(q32_catalogue())[index]
    names = (Path(__file__).name, 'coherent_weighted_gauges.py', 'coherent_hot_mixed.py',
             'coherent_hot_mixed.cpp', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds_fast.py', 'coherent_hot_compounds.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp', 'selected_gaussian_coefficients.cpp')
    report = {'scope': __doc__, 'completed': False, 'description': description,
              'embedding': p, 'retraction': l, 'modes': modes, 'records': [],
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'seconds_budget': args.seconds, 'iteration_budget': args.maxiter,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'exact_pressure_certificate': False, 'global_optimality_claimed': False}
    _atomic_write(args.output, report)
    started = time.monotonic()
    oracle = SharedGaugeOracle(p, l, modes)
    report.update(deletion_edges=oracle.edges, cone=oracle.cone.astype(int).tolist(),
                  natural=oracle.natural, free_labels=oracle.free)
    last = None
    queries = 0
    best = float('inf')

    def objective(parameters):
        nonlocal last, queries, best
        if time.monotonic() - started > args.seconds:
            raise TimeoutError('The requested exploratory time budget ended.')
        if last is not None and np.array_equal(parameters, last[0]):
            return last[1]['rho'], last[1]['gradient']
        value = oracle.evaluate(parameters, enforce=False, details=True)
        queries += 1
        last = (np.array(parameters), value)
        if value['minimum_cone_slack'] >= -1e-8 and value['rho'] < best:
            best = value['rho']
            report['best'] = {key: item for key, item in value.items() if key != 'gradient'}
            report['best']['parameters'] = list(map(float, parameters))
            report['queries'] = queries
            report['elapsed_seconds'] = time.monotonic() - started
            _atomic_write(args.output, report)
            print(json.dumps({'rho': best, 'queries': queries,
                              'slack': value['minimum_cone_slack'],
                              'seconds': report['elapsed_seconds']}), flush=True)
        return value['rho'], value['gradient']

    start = np.zeros(oracle.q)
    start[0] = .3295
    random = np.random.default_rng(104729)
    perturbed = start.copy()
    perturbed[1:] = .001 * oracle.heights[oracle.free] + random.uniform(-.0002, .0002, oracle.q - 1)
    starts = [start, perturbed]
    for initial in starts:
        try:
            fit = minimize(objective, initial, method='SLSQP', jac=True,
                           constraints=LinearConstraint(oracle.cone, 0., np.inf),
                           bounds=[(.01, 2.)] + [(-10., 10.)] * (oracle.q - 1),
                           options={'maxiter': args.maxiter, 'ftol': 1e-11})
            report['records'].append({'rho': float(fit.fun), 'parameters': fit.x.tolist(),
                                       'success': bool(fit.success), 'message': str(fit.message),
                                       'iterations': int(fit.nit), 'queries': int(fit.nfev)})
        except (TimeoutError, ArithmeticError, ValueError) as error:
            report['records'].append({'error': str(error), 'type': type(error).__name__})
            if isinstance(error, TimeoutError):
                break
        _atomic_write(args.output, report)
    report.update(completed=True, queries=queries, elapsed_seconds=time.monotonic() - started,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
