#!/usr/bin/env python3
"""Registered free-height asymmetric search on retained F4 pencil classes."""
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

from check_scaled_asymmetric_pressure import inputs
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest,cases = inputs(root)
    names = ['tools/search_two_place_free_asymmetric.py','tools/check_scaled_asymmetric_pressure.py',
             'tools/scaled_asymmetric_crossing_pressure.py','tools/asymmetric_crossing_pressure.py',
             'tools/native_rectangular_smith_f4.py','tools/native_rectangular_smith_f4.cpp',
             'tools/native_weighted_smith_f4.py','tools/two_place_f4_census.py',
             'tools/joint_energy_conjugations.py','tools/search_f4_recentered_crossings.py']
    report = {'status':'registered_incomplete','manifest':manifest,'trials':[],
              'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_tree':subprocess.check_output(['git','rev-parse','HEAD^{tree}'],cwd=root,text=True).strip(),
              'input_profiles_verified':[],
              'successful_queries':0,'rejected_queries':0,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names},
              'best':None,'exact_contraction_certificate':False,
              'scope':manifest['scope']}
    started = time.monotonic()
    deadline = started+manifest['global_seconds']
    rng = np.random.default_rng(manifest['seed'])
    natural = np.asarray(manifest['natural_heights_and_scale'],float)
    previous = {}

    def save():
        report['seconds'] = time.monotonic()-started
        args.output.parent.mkdir(parents=True,exist_ok=True)
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(report,indent=2)+'\n')
        temporary.replace(args.output)
    save()
    bounds = [tuple(manifest['bounds']['H'])]*32
    bounds += [tuple(manifest['bounds']['log_t']),tuple(manifest['bounds']['log_tan'])]
    angles = manifest['restart_angles']
    exhausted = False
    for restart,angle in enumerate(angles):
        for case,values in cases:
            if time.monotonic()>=deadline:
                exhausted = True
                break
            index = case['index']
            model = ScaledAsymmetricPressure({'q_x':5,'q_u':5,'values':values})
            if restart == 0:
                assert model.base.oracle.prefixes(natural,1)['orders'] == case['natural_orders']
                report['input_profiles_verified'].append(index)
            if restart and index in previous:
                initial = np.array(previous[index]['optimizer_variables'],float)
                initial[:-2] += rng.normal(size=32)*(.03 if restart%2 else .25)
            else:
                initial = np.r_[natural[:-1],math.log(.2),0.]
            initial[-1] = math.log(math.tan(angle))
            initial = np.maximum([b[0] for b in bounds],np.minimum([b[1] for b in bounds],initial))
            trial_started = time.monotonic()
            trial_deadline = min(deadline,trial_started+manifest['per_trial_seconds'])
            local = {'best':None,'successful_queries':0,'rejected_queries':0}

            def objective(point):
                if time.monotonic()>=trial_deadline:
                    raise TimeoutError('The registered trial or global deadline is exhausted.')
                try:
                    value,gradient,data = model.evaluate(point)
                except (ArithmeticError,ValueError):
                    local['rejected_queries'] += 1
                    report['rejected_queries'] += 1
                    return 1e6,np.zeros_like(point)
                local['successful_queries'] += 1
                report['successful_queries'] += 1
                if local['best'] is None or value<local['best']['support_value']:
                    local['best'] = data
                if report['best'] is None or value<report['best']['support_value']-1e-10:
                    report['best'] = dict(data,case=case,restart=restart,literal_values=values)
                    save()
                return value,gradient
            completed,error = True,None
            try:
                fitted = minimize(objective,initial,jac=True,method='SLSQP',bounds=bounds,
                                  options={'maxiter':manifest['iteration_ceiling'],'ftol':1e-11})
                outcome = {'success':bool(fitted.success),'message':str(fitted.message),
                           'iterations':int(fitted.nit),'evaluations':int(fitted.nfev)}
            except TimeoutError as exception:
                completed,error,outcome = False,str(exception),None
            if local['best'] is not None:
                previous[index] = local['best']
            record = {'case_index':index,'representative':case['representative'],'restart':restart,
                      'initial_variables':initial.tolist(),'completed_optimizer_call':completed,
                      'interruption':error,'outcome':outcome,**local,
                      'seconds':time.monotonic()-trial_started}
            report['trials'].append(record)
            report['status'] = 'running_incomplete'
            save()
            print(json.dumps({'case':index,'restart':restart,'best':None if local['best'] is None else
                local['best']['support_value'],'global_best':report['best']['support_value'] if report['best'] else None,
                'queries':report['successful_queries'],'seconds':report['seconds']}),flush=True)
        if exhausted:
            break
    report['status'] = 'complete_global_time_budget' if exhausted else 'complete_registered_trial_list'
    report['registered_trials'] = len(cases)*len(angles)
    report['completed_trial_records'] = len(report['trials'])
    report['completed_optimizer_calls'] = sum(t['completed_optimizer_call'] for t in report['trials'])
    save()
    print(json.dumps({'status':report['status'],'attempted_trials':len(report['trials']),
                      'best':report['best']['support_value'] if report['best'] else None,
                      'seconds':report['seconds']}),flush=True)


if __name__ == '__main__':
    main()
