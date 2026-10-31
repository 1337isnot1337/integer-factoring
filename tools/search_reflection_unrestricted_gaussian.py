#!/usr/bin/env python3
"""Explore independent reflection Gaussian gauges without a regularity cone.

The Clifford profile reduction is proved separately. This floating search
optimizes the actual Gibbs root with both physical type costs. It cannot
certify an improvement or a universal obstruction on its own.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from reflection_unrestricted_gaussian import ReflectionGaussianPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=16)
    parser.add_argument('--restarts', type=int, default=20)
    parser.add_argument('--iterations', type=int, default=900)
    parser.add_argument('--seconds', type=float, default=3600.)
    parser.add_argument('--seed', type=int, default=611703)
    parser.add_argument('--precision', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.restarts < 1 or args.iterations < 1 or args.seconds <= 0:
        raise ValueError('Use positive restart, iteration, and time budgets.')
    pressure = ReflectionGaussianPressure(args.interval, args.precision)
    rng = np.random.default_rng(args.seed)
    fit = minimize_scalar(lambda x: pressure.evaluate(pressure.natural(math.exp(x))),
                          bounds=(-8., math.log(2.)), method='bounded')
    scale = math.exp(fit.x)
    natural = pressure.natural(scale)
    pressure.evaluate(natural)
    baseline = pressure.last.copy()
    best = baseline.copy()
    # A mode perturbation makes the control less likely to lie on a
    # weighted-Smith wall. The scale step accounts for M**2 base orders.
    control = natural+rng.normal(size=len(natural))*.0007
    control[-1] = scale
    _, analytic = pressure.evaluate(control, True)
    gradient_error = 0.
    coordinates = list(range(min(17, len(control)-1)))+[len(control)-1]
    for coordinate in coordinates:
        step = 1e-6 if coordinate+1 < len(control) else 1e-6/(args.interval*args.interval)
        delta = np.zeros_like(control)
        delta[coordinate] = step
        finite = (pressure.evaluate(control+delta)-pressure.evaluate(control-delta))/(2*step)
        gradient_error = max(gradient_error, abs(float(finite-analytic[coordinate])))
    if gradient_error > 2e-5:
        raise ArithmeticError(f'The local gradient control is unresolved: {gradient_error}.')
    start = time.monotonic()
    deadline = start+args.seconds
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/reflection_unrestricted_gaussian.py',
               'tools/search_reflection_unrestricted_gaussian.py',
               'research/clifford-gaussian-smith-profiles.md']
    result = {'status': 'exploratory floating unrestricted Gaussian gauge search',
              'interval': args.interval, 'independent_heights': 4*args.interval,
              'baseline': baseline, 'best': best, 'runs': [],
              'gradient_max_error': gradient_error, 'random_seed': args.seed,
              'rejected_queries': 0, 'budget_exhausted': False,
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in sources},
              'scope': 'The full Clifford/Fock profile theorem is exact; this finite numerical optimization proves no strict improvement or universal exclusion.'}

    def save():
        result['best'] = best
        result['queries'] = pressure.queries
        result['seconds'] = time.monotonic()-start
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        temporary.replace(args.output)

    save()
    print(json.dumps({'interval': args.interval, 'baseline_rho': baseline['rho'],
                      'scale': scale, 'gradient_max_error': gradient_error}), flush=True)
    for trial in range(args.restarts):
        if time.monotonic() >= deadline:
            result['budget_exhausted'] = True
            break
        initial = np.r_[natural[:-1]/(scale*args.interval), math.log(scale)]
        if trial:
            initial[:-1] += rng.normal(size=len(initial)-1)*(.1 if trial < args.restarts//2 else .75)
        initial_scale = math.exp(initial[-1])
        initial_parameters = np.r_[initial_scale*args.interval*initial[:-1], initial_scale]
        initial_margin = pressure.internal_margin(initial_parameters)

        def objective(parameters):
            nonlocal best
            if time.monotonic() >= deadline:
                raise TimeoutError('Requested search budget exhausted.')
            current_scale = math.exp(parameters[-1])
            current = np.r_[current_scale*args.interval*parameters[:-1], current_scale]
            try:
                value, derivative = pressure.evaluate(current, True)
            except ArithmeticError:
                result['rejected_queries'] += 1
                return 2., np.zeros_like(parameters)
            if value < best['rho']-1e-13:
                best = pressure.last.copy()
                save()
            transformed = np.r_[derivative[:-1]*current_scale*args.interval, derivative@current]
            return value, transformed

        try:
            fitted = minimize(objective, initial, jac=True, method='SLSQP',
                              bounds=[(-8., 8.)]*(len(initial)-1)+[(math.log(.0001), math.log(2.))],
                              options={'maxiter': args.iterations, 'ftol': 1e-11})
            result['runs'].append({'trial': trial, 'rho': float(fitted.fun),
                                   'success': bool(fitted.success), 'iterations': int(fitted.nit),
                                   'initial_internal_regularity_margin': initial_margin,
                                   'message': str(fitted.message)})
        except TimeoutError:
            result['budget_exhausted'] = True
            break
        save()
        print(json.dumps({'trial': trial, 'best_rho': best['rho'],
                          'queries': pressure.queries, 'seconds': result['seconds']}), flush=True)
    save()


if __name__ == '__main__':
    main()
