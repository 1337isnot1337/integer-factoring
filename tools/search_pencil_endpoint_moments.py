#!/usr/bin/env python3
"""Registered endpoint second-moment search, with complete profiles throughout."""
from hashlib import sha256
import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import minimize

from pencil_endpoint_moments import EndpointMoments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root,started = Path(__file__).resolve().parents[1],time.monotonic()
    manifest = json.loads(args.manifest.read_text())
    names = ['tools/search_pencil_endpoint_moments.py','tools/pencil_endpoint_moments.py',
             'tools/pencil_principal_alphabets.py','tools/search_f4_recentered_crossings.py',
             'tools/native_rectangular_smith_f4.py','tools/native_rectangular_smith_f4.cpp',
             'tools/native_weighted_smith_f4.py']
    report = {'status':'registered_incomplete','registration':manifest,'trials':[],
              'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names},'best_by_axis':{},
              'successful_queries':0,'rejected_queries':0,
              'scope':'Numerical variance search on a saturated simple-minimum face. '
                      'A nonpositive numerical value needs independent exact replay; '
                      'positive values give no global lower bound or score improvement.'}
    model = EndpointMoments()
    def save():
        report['seconds'] = time.monotonic()-started
        args.output.parent.mkdir(parents=True,exist_ok=True)
        temp = args.output.with_suffix('.tmp')
        temp.write_text(json.dumps(report,indent=2)+'\n')
        temp.replace(args.output)
    save()
    try:
        for strong in manifest['strong_places']:
            rng = np.random.default_rng(manifest['seed']+strong)
            natural = np.asarray(model.natural[12*strong:12*(strong+1)])
            axis_start = time.monotonic()
            axis_deadline = axis_start+manifest['seconds_per_axis']
            for trial in range(manifest['restarts_ceiling']):
                if time.monotonic()>=axis_deadline:
                    report['trials'].append({'strong_place':strong,'status':'complete_axis_time_budget'})
                    break
                sigma = 0 if trial==0 else manifest['perturbation_sigmas'][(trial-1)%len(manifest['perturbation_sigmas'])]
                initial = np.clip(natural+rng.normal(size=12)*sigma,-manifest['height_bound'],manifest['height_bound'])
                local = {'best':None,'queries':0,'rejections':0}
                def objective(point):
                    if time.monotonic()>=axis_deadline:
                        raise TimeoutError
                    try:
                        value,gradient,data = model.evaluate(point,strong)
                    except ArithmeticError as error:
                        if 'locus changed' in str(error):
                            raise
                        report['rejected_queries'] += 1
                        local['rejections'] += 1
                        return 1e6,np.zeros_like(point)
                    report['successful_queries'] += 1
                    local['queries'] += 1
                    if local['best'] is None or value<local['best']['variance_difference']:
                        local['best'] = data
                    key = str(strong)
                    if key not in report['best_by_axis'] or value<report['best_by_axis'][key]['variance_difference']-1e-10:
                        report['best_by_axis'][key] = data
                        save()
                    if value<1e-8:
                        save()
                        raise StopIteration
                    return value,gradient
                completed,outcome = True,None
                try:
                    fit = minimize(objective,initial,jac=True,method='SLSQP',
                        bounds=[(-manifest['height_bound'],manifest['height_bound'])]*12,
                        options={'maxiter':manifest['iterations_ceiling'],'ftol':1e-11})
                    outcome = {'success':bool(fit.success),'message':str(fit.message),
                               'iterations':int(fit.nit),'evaluations':int(fit.nfev)}
                except TimeoutError:
                    completed = False
                report['trials'].append(dict(strong_place=strong,restart=trial,initial=initial.tolist(),
                    completed_optimizer_call=completed,outcome=outcome,**local))
                report['status'] = 'running_incomplete'
                save()
                print(json.dumps({'strong_place':strong,'restart':trial,
                    'best':report['best_by_axis'][str(strong)]['variance_difference'],
                    'queries':report['successful_queries'],'seconds':report['seconds']}),flush=True)
        report['status'] = 'complete_registered_limits'
    except StopIteration:
        report['status'] = 'candidate_requires_exact_replay'
    save()


if __name__=='__main__':
    main()
