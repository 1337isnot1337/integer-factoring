#!/usr/bin/env python3
"""Bounded search with independent complete-cut gauges and a common scale."""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import LinearConstraint, minimize

from coherent_asymmetric_gauges import AsymmetricGaugeOracle
from search_coherent_hot_retractions import _atomic_write, q32_catalogue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=int, default=55)
    parser.add_argument('--maxiter', type=int, default=100)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    description, embedding, retraction, modes = list(q32_catalogue())[args.index]
    started = time.monotonic()
    oracle = AsymmetricGaugeOracle(embedding, retraction, modes)
    report = {'scope': __doc__, 'completed': False, 'index': args.index,
              'description': description, 'embedding': embedding, 'retraction': retraction,
              'modes': modes, 'free_labels': oracle.free, 'records': [],
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'seconds_budget': args.seconds, 'iteration_budget': args.maxiter,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'coherent_asymmetric_gauges.py',
                                             'coherent_weighted_gauges.py')},
              'exact_pressure_certificate': False, 'global_optimality_claimed': False}
    _atomic_write(args.output, report)
    last, queries, best = None, 0, float('inf')

    def objective(parameters):
        nonlocal last, queries, best
        if time.monotonic() - started > args.seconds:
            raise TimeoutError('The exploratory time budget ended.')
        if last is not None and np.array_equal(parameters, last[0]):
            return last[1]['rho'], last[1]['gradient']
        value = oracle.evaluate(parameters, enforce=False, details=True)
        last, queries = (np.array(parameters), value), queries + 1
        if value['minimum_cone_slack'] >= -1e-8 and value['rho'] < best:
            best = value['rho']
            report['best'] = {key: item for key, item in value.items() if key != 'gradient'}
            report['best']['parameters'] = list(map(float, parameters))
            report.update(queries=queries, elapsed_seconds=time.monotonic() - started)
            _atomic_write(args.output, report)
            print(json.dumps({'rho': best, 'queries': queries,
                              'seconds': report['elapsed_seconds']}), flush=True)
        return value['rho'], value['gradient']

    initial = np.zeros(2 * oracle.q - 1)
    initial[0] = .3295
    starts = [initial]
    rng = np.random.default_rng(202611014)
    for amount in (.001, .015):
        other = initial.copy()
        for place in (1, oracle.q):
            other[place:place + oracle.q - 1] = (amount * oracle.heights[oracle.free]
                                                + rng.uniform(-amount / 10, amount / 10, oracle.q - 1))
        starts.append(other)
    for initial in starts:
        try:
            fit = minimize(objective, initial, method='SLSQP', jac=True,
                           constraints=LinearConstraint(oracle.asymmetric_cone, 0., np.inf),
                           bounds=[(.01, 2.)] + [(-10., 10.)] * (2 * oracle.q - 2),
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
