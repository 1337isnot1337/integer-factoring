#!/usr/bin/env python3
"""Numerical whole-symbol search on coordinate Pascal crossings over F4."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from f4_pascal_coordinate_sectors import F4SectorPascalLetterPressure
from f4_pascal_flags import change_flag
from permuted_pascal_kernels import permuted_pascal_kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--images', type=int, nargs='+', required=True)
    parser.add_argument('--alphabet', type=int, nargs='+', required=True)
    parser.add_argument('--flag', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--iterations', type=int, default=500)
    parser.add_argument('--restarts', type=int, default=12)
    parser.add_argument('--precision', type=int, default=512)
    parser.add_argument('--cached-profiles', action='store_true',
                        help='Check adaptive regular LDU cones before rebuilding a profile.')
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.seconds, args.iterations, args.restarts) <= 0:
        parser.error('Use positive search budgets.')
    root, start = Path(__file__).resolve().parents[1], time.monotonic()
    names = ['tools/search_f4_pascal_coordinate_flags.py', 'tools/f4_pascal_coordinate_sectors.py',
             'tools/f4_pascal_flags.py', 'tools/pascal_coordinate_sectors.py',
             'tools/pascal_letter_gauges.py', 'tools/permuted_pascal_kernels.py',
             'tools/cube_f4_bivariate_smith.py', 'tools/native_matrix_smith_f4.py',
             'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp',
             'tools/native_weighted_smith_f4.py', 'tools/weighted_pascal_flags.py',
             'tools/tensor_pascal_box_crossings.py']
    if args.cached_profiles:
        names += ['tools/native_smith_chambers_f4.py','tools/native_smith_chambers_f4.cpp',
                  'tools/native_matrix_smith_f4.cpp']
    model = None
    flag = json.loads(args.flag.read_text())
    report = {'status':'registered_incomplete', 'images':args.images, 'alphabet':args.alphabet,
              'constant_f4_flag':flag, 'seed':args.seed, 'seconds_budget':args.seconds,
              'cached_profiles':args.cached_profiles,
              'restarts_ceiling':args.restarts, 'iterations_ceiling':args.iterations,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names}, 'runs':[],
              'rejected_queries':0, 'scope':'Numerical complete original-polynomial '
              'profiles with both physical costs. No exact contraction certificate.'}

    def save():
        report['total_seconds'] = time.monotonic()-start
        if model is not None:
            report['queries'] = model.queries
            if args.cached_profiles:
                from native_smith_chambers_f4 import sector_cache_statistics
                report['profile_cache_statistics'] = sector_cache_statistics(model)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.output.with_suffix('.tmp')
        tmp.write_text(json.dumps(report,indent=2)+'\n')
        tmp.replace(args.output)
    save()
    seed = change_flag(permuted_pascal_kernel(args.images), flag)
    model = F4SectorPascalLetterPressure(seed, args.alphabet, precision=args.precision)
    if args.cached_profiles:
        from native_smith_chambers_f4 import cache_sector_pressure
        cache_sector_pressure(model)
    report['construction_seconds'] = time.monotonic()-start
    report['construction'] = model.construction
    natural = model.fitted_natural_start()
    # Remove the eight constant-height gauge directions. This preserves the
    # restored pressure and both physical costs exactly up to scalar shifts.
    heights = natural[:-1].reshape(8,model.q)
    heights -= heights[:, :1]
    model.evaluate(natural)
    best = model.last.copy()
    report.update(baseline=best, best=best)
    rng = np.random.default_rng(args.seed)
    control = natural.copy()
    control[:-1] += rng.normal(size=len(control)-1)*.0001
    _, derivative = model.evaluate(control, True)
    gradients = []
    for coordinate in (0,1,model.q,4*model.q,8*model.q):
        step = 1e-7/(seed['outer']*seed['n']) if coordinate == 8*model.q else 1e-7
        delta = np.zeros_like(control)
        delta[coordinate] = step
        finite = (model.evaluate(control+delta)-model.evaluate(control-delta))/(2*step)
        error = abs(float(finite-derivative[coordinate]))
        if error > 1e-4:
            raise ArithmeticError('The independent field-valued symbol derivative control failed.')
        gradients.append({'coordinate':coordinate, 'analytic':float(derivative[coordinate]),
                          'finite_difference':float(finite), 'error':error})
    report['gradient_controls'] = gradients
    search_start = time.monotonic()
    report['presearch_seconds'] = search_start-start
    bounds = []
    for i in range(8*model.q):
        bounds.append((0.,0.) if i % model.q == 0 else (-12.,12.))
    bounds.append((math.log(.0001), math.log(2.)))
    report['status'] = 'running_incomplete'
    save()
    print(json.dumps({'baseline':best['rho'], 'construction':report['construction'],
                      'presearch_seconds':report['presearch_seconds']}),flush=True)
    for trial in range(args.restarts):
        if time.monotonic()-search_start >= args.seconds:
            report['status'] = 'complete_time_budget'
            break
        scale = natural[-1]
        initial = np.r_[natural[:-1]/(scale*seed['outer']), math.log(scale)]
        if trial:
            initial[:-1] += rng.normal(size=len(initial)-1)*(.02 if trial % 2 else .15)
            initial[np.arange(8)*model.q] = 0.

        def objective(parameters):
            nonlocal best
            if time.monotonic()-search_start >= args.seconds:
                raise TimeoutError('The registered search budget is exhausted.')
            scale = math.exp(parameters[-1])
            current = np.r_[scale*seed['outer']*parameters[:-1], scale]
            try:
                value, derivative = model.evaluate(current, True)
            except ArithmeticError:
                report['rejected_queries'] += 1
                return 2., np.zeros_like(parameters)
            if value < best['rho']-1e-12:
                best = model.last.copy()
                report['best'] = best
                save()
            elif model.queries % 25 == 0:
                report['search_seconds'] = time.monotonic()-search_start
                save()
            transformed = np.r_[derivative[:-1]*scale*seed['outer'], derivative @ current]
            return value, transformed
        try:
            fitted = minimize(objective, initial, jac=True, method='SLSQP', bounds=bounds,
                              options={'maxiter':args.iterations, 'ftol':1e-11})
            report['runs'].append({'trial':trial, 'rho':float(fitted.fun),
                                  'iterations':int(fitted.nit), 'success':bool(fitted.success),
                                  'message':str(fitted.message)})
        except TimeoutError:
            report['status'] = 'complete_time_budget'
            break
        report.update(queries=model.queries, search_seconds=time.monotonic()-search_start)
        save()
        print(json.dumps({'trial':trial,'best':best['rho'],'queries':model.queries,
                          'search_seconds':report['search_seconds']}),flush=True)
    else:
        report['status'] = 'complete_restart_ceiling'
    report.update(queries=model.queries, search_seconds=time.monotonic()-search_start, best=best)
    save()
    print(json.dumps({'status':report['status'],'best':best['rho'],
                      'queries':model.queries,'search_seconds':report['search_seconds']}),flush=True)


if __name__ == '__main__':
    main()
