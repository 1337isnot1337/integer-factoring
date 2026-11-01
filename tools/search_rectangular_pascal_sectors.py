#!/usr/bin/env python3
"""Registered whole-symbol follow-on for a completed rectangular screen."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import minimize

from permuted_pascal_kernels import permuted_pascal_kernel
from rectangular_pascal_sectors import RectangularPascalPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, started = Path(__file__).resolve().parents[1], time.monotonic()
    screen = json.loads(args.screen.read_text())
    if screen['status'] != 'complete_registered_screen':
        raise ValueError('The registered screen must finish before selection.')
    eligible = [r for r in screen['completed'] if r['case']['alphabet_b'] != r['case']['alphabet_x']]
    selected = min(eligible,key=lambda r:(r['balanced']['rho'],r['case']['label']))
    case, policy = selected['case'], screen['registration']['follow_on_policy']
    names = list(screen['source_sha256'])+['tools/search_rectangular_pascal_sectors.py']
    seed = permuted_pascal_kernel(screen['registration']['images'])
    model = None
    report = {'status':'registered_incomplete','policy':policy,'selected_case':case,
              'screen_sha256':sha256(args.screen.read_bytes()).hexdigest(),
              'screen_source_commit':screen['git_head'],
              'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names},
              'runs':[],'rejected_queries':0,'scope':'Numerical arbitrary-symbol heights on the '
                    'registered selected rectangular crossing, complete original-polynomial '
                    'profiles and both physical costs. No exact contraction certificate.'}
    def save():
        report['total_seconds'] = time.monotonic()-started
        if model is not None:
            report['queries'] = model.queries
        args.output.parent.mkdir(parents=True,exist_ok=True)
        tmp = args.output.with_suffix('.tmp')
        tmp.write_text(json.dumps(report,indent=2)+'\n')
        tmp.replace(args.output)
    save()
    model = RectangularPascalPressure(seed,case['alphabet_b'],case['alphabet_x'],
                                     precision=selected['attempts'][-1]['precision'])
    natural = np.asarray(selected['balanced']['parameters'])
    constant_coordinates = [p*model.size+o for p in range(2) for o in model.offsets]
    for p in range(2):
        for leg,o in enumerate(model.offsets):
            segment = natural[p*model.size+o:p*model.size+o+len(model.leg_alphabets[leg])]
            segment -= segment[0]
    model.evaluate(natural)
    if abs(model.last['rho']-selected['balanced']['rho']) > 1e-10:
        raise ArithmeticError('The frozen selected baseline did not replay.')
    report['construction'] = model.construction
    best = model.last.copy()
    report.update(baseline=best,best=best)
    rng = np.random.default_rng(policy['seed'])
    control = natural.copy()
    control[:-1] += rng.normal(size=len(control)-1)*.0001
    _, derivative = model.evaluate(control,True)
    gradients = []
    for coordinate in (0,model.a,model.a+model.b,model.size,2*model.size):
        step = 1e-7/(seed['outer']*seed['n']) if coordinate == 2*model.size else 1e-7
        delta = np.zeros_like(control)
        delta[coordinate] = step
        finite = (model.evaluate(control+delta)-model.evaluate(control-delta))/(2*step)
        error = abs(float(finite-derivative[coordinate]))
        if error > 1e-4:
            raise ArithmeticError('A selected rectangular derivative control failed.')
        gradients.append({'coordinate':coordinate,'analytic':float(derivative[coordinate]),
                          'finite_difference':float(finite),'error':error})
    report['gradient_controls'] = gradients
    bounds = [(0.,0.) if i in constant_coordinates else (-12.,12.) for i in range(2*model.size)]
    bounds.append((math.log(.0001),math.log(2.)))
    search_start = time.monotonic()
    report.update(status='running_incomplete',presearch_seconds=search_start-started)
    save()
    print(json.dumps({'case':case['label'],'baseline':best['rho'],
                      'presearch_seconds':report['presearch_seconds']}),flush=True)
    for trial in range(policy['restarts_ceiling']):
        if time.monotonic()-search_start >= policy['whole_symbol_search_seconds']:
            report['status'] = 'complete_time_budget'
            break
        scale = natural[-1]
        initial = np.r_[natural[:-1]/(scale*seed['outer']),math.log(scale)]
        if trial:
            initial[:-1] += rng.normal(size=len(initial)-1)*(.02 if trial % 2 else .15)
            initial[constant_coordinates] = 0.
        def objective(parameters):
            nonlocal best
            if time.monotonic()-search_start >= policy['whole_symbol_search_seconds']:
                raise TimeoutError('The registered whole-symbol budget is exhausted.')
            scale = math.exp(parameters[-1])
            current = np.r_[scale*seed['outer']*parameters[:-1],scale]
            try:
                value, derivative = model.evaluate(current,True)
            except ArithmeticError:
                report['rejected_queries'] += 1
                return 2.,np.zeros_like(parameters)
            if value < best['rho']-1e-12:
                best = model.last.copy()
                report['best'] = best
                save()
            elif model.queries % 25 == 0:
                report['search_seconds'] = time.monotonic()-search_start
                save()
            return value,np.r_[derivative[:-1]*scale*seed['outer'],derivative @ current]
        try:
            fit = minimize(objective,initial,jac=True,method='SLSQP',bounds=bounds,
                           options={'maxiter':policy['iterations_ceiling'],'ftol':1e-11})
            report['runs'].append({'trial':trial,'rho':float(fit.fun),'iterations':int(fit.nit),
                                  'success':bool(fit.success),'message':str(fit.message)})
        except TimeoutError:
            report['status'] = 'complete_time_budget'
            break
        report['search_seconds'] = time.monotonic()-search_start
        save()
        print(json.dumps({'trial':trial,'best':best['rho'],'queries':model.queries}),flush=True)
    else:
        report['status'] = 'complete_restart_ceiling'
    report.update(best=best,search_seconds=time.monotonic()-search_start)
    save()
    print(json.dumps({'status':report['status'],'best':best['rho'],'queries':model.queries,
                      'search_seconds':report['search_seconds']}),flush=True)


if __name__ == '__main__':
    main()
