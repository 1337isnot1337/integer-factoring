#!/usr/bin/env python3
"""Execute the registered asymmetric follow-on on proper pencil restrictions."""
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
from pencil_principal_alphabets import restrict
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    summary = json.loads((args.screen/'summary.json').read_text())
    policy = json.loads(args.manifest.read_text())
    assert summary['status']=='complete_registered_screen'
    assert summary['completed_cases']==summary['registration']['registered_cases']==45292
    assert summary['resolved_orders']==753750 and not summary['unresolved_profiles']
    assert sha256((args.screen/'profiles.jsonl').read_bytes()).hexdigest()==summary['profile_ledger_sha256']
    assert sha256((root/'evidence/pencil-principal-alphabets-manifest.json').read_bytes()).hexdigest()==policy['screen_manifest_sha256']
    fits = json.loads((args.screen/'scalar-fits.json').read_text())
    records = [json.loads(line) for line in (args.screen/'profiles.jsonl').read_text().splitlines()]
    records = sorted((r for r in records if len(r['alphabet_x'])<5 or len(r['alphabet_y'])<5),
        key=lambda r:(fits[r['scalar_fit_index']]['rho'],r['source_case_index'],r['alphabet_x'],r['alphabet_y']))
    selected,seen = [],set()
    for record in records:
        key = len(record['alphabet_x']),len(record['alphabet_y']),record['literal_values_sha256']
        if key not in seen:
            selected.append(record)
            seen.add(key)
            if len(selected)==policy['selected_ceiling']:
                break
    source,cases = inputs(root)
    values_by_index = {c['index']:v for c,v in cases}
    names = ['tools/search_pencil_principal_asymmetric.py','tools/pencil_principal_alphabets.py',
             'tools/check_scaled_asymmetric_pressure.py','tools/scaled_asymmetric_crossing_pressure.py',
             'tools/asymmetric_crossing_pressure.py','tools/native_rectangular_smith_f4.py',
             'tools/native_rectangular_smith_f4.cpp','tools/native_weighted_smith_f4.py',
             'tools/two_place_f4_census.py','tools/joint_energy_conjugations.py',
             'tools/search_f4_recentered_crossings.py']
    report = {'status':'registered_incomplete','registration':policy,'selected':selected,
              'screen_source_commit':summary['source_commit'],'screen_ledger_sha256':summary['profile_ledger_sha256'],
              'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names},'trials':[],
              'successful_queries':0,'rejected_queries':0,'best':None,'verified_initial_profiles':[],
              'scope':'Registered proper-subalphabet local free-height asymmetric search. '
                      'Exact native field arithmetic, floating weights and pressure, no global exclusion.'}
    started,previous = time.monotonic(),{}
    rng = np.random.default_rng(policy['seed'])
    deadline = started+policy['global_seconds']
    def save():
        report['seconds'] = time.monotonic()-started
        args.output.parent.mkdir(parents=True,exist_ok=True)
        temp = args.output.with_suffix('.tmp')
        temp.write_text(json.dumps(report,indent=2)+'\n')
        temp.replace(args.output)
    save()
    exhausted = False
    for restart,angle in enumerate(policy['restart_angles']):
        for selection,record in enumerate(selected):
            if time.monotonic()>=deadline:
                exhausted = True
                break
            seed,natural,_,_ = restrict(values_by_index[record['source_case_index']],
                source['natural_heights_and_scale'],record['alphabet_x'],record['alphabet_y'])
            assert sha256(json.dumps(seed['values'],separators=(',',':')).encode()).hexdigest()==record['literal_values_sha256']
            model = ScaledAsymmetricPressure(seed)
            if restart==0:
                assert model.base.oracle.prefixes(natural,1)['orders']==sum(record['orders'],[])
                report['verified_initial_profiles'].append(selection)
            height_count = len(natural)-1
            if restart and selection in previous:
                initial = np.asarray(previous[selection]['optimizer_variables'])
                initial[:-2] += rng.normal(size=height_count)*(.03 if restart%2 else .25)
            else:
                initial = np.r_[natural[:-1],math.log(.2),0.]
            initial[-1] = math.log(math.tan(angle))
            bounds = [tuple(policy['bounds']['H'])]*height_count
            bounds += [tuple(policy['bounds']['log_t']),tuple(policy['bounds']['log_tan'])]
            initial = np.clip(initial,[b[0] for b in bounds],[b[1] for b in bounds])
            trial_start = time.monotonic()
            trial_deadline = min(deadline,trial_start+policy['per_trial_seconds'])
            local = {'best':None,'successful_queries':0,'rejected_queries':0}
            def objective(point):
                if time.monotonic()>=trial_deadline:
                    raise TimeoutError
                try:
                    value,gradient,data = model.evaluate(point)
                except (ArithmeticError,ValueError):
                    report['rejected_queries'] += 1
                    local['rejected_queries'] += 1
                    return 1e6,np.zeros_like(point)
                report['successful_queries'] += 1
                local['successful_queries'] += 1
                if local['best'] is None or value<local['best']['support_value']:
                    local['best'] = data
                if report['best'] is None or value<report['best']['support_value']-1e-10:
                    report['best'] = dict(data,selection=selection,restricted_record=record,literal_values=seed['values'])
                    save()
                return value,gradient
            completed = True
            try:
                fit = minimize(objective,initial,jac=True,method='SLSQP',bounds=bounds,
                    options={'maxiter':policy['iteration_ceiling'],'ftol':1e-11})
                outcome = {'success':bool(fit.success),'message':str(fit.message),
                           'iterations':int(fit.nit),'evaluations':int(fit.nfev)}
            except TimeoutError:
                completed,outcome = False,None
            if local['best'] is not None:
                previous[selection] = local['best']
            report['trials'].append(dict(selection=selection,restart=restart,
                completed_optimizer_call=completed,outcome=outcome,seconds=time.monotonic()-trial_start,**local))
            report['status'] = 'running_incomplete'
            save()
            print(json.dumps({'selection':selection,'restart':restart,
                'best':local['best']['support_value'] if local['best'] else None,
                'global_best':report['best']['support_value'] if report['best'] else None,
                'queries':report['successful_queries'],'seconds':report['seconds']}),flush=True)
        if exhausted:
            break
    report.update(status='complete_global_time_budget' if exhausted else 'complete_registered_trial_list',
                  registered_trials=len(selected)*len(policy['restart_angles']))
    save()


if __name__ == '__main__':
    main()
