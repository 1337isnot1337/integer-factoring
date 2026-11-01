#!/usr/bin/env python3
"""Adapt the existing fixed-type proximal bundle method to Fock coordinates.

The source method is search_rectangular_core_bundle.py. Every new supporting
plane is checked against every earlier queried point in both directions.
Numerical master solutions and stopping conditions certify no global bound.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize
from scipy.special import logsumexp

from permuted_pascal_kernels import permuted_pascal_kernel
from rectangular_pascal_sectors import RectangularPascalPressure


def physical_vector(model, probabilities):
    return np.r_[[v for _ in range(2) for leg in range(4)
                  for v in probabilities[0 if leg in (0,3) else 1]],0.]


def log_profile_query(model, parameters):
    rho, derivative = model.evaluate(parameters,True)
    detail = model.last.copy()
    physical = physical_vector(model,detail['pair_probabilities'])
    value = sum(float(logsumexp(-np.asarray(c['orders']))) for c in detail['complete_original_profiles'])
    gradient = 2*detail['entropy']*derivative-physical
    return rho,value,gradient,detail,physical


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root, begun = Path(__file__).resolve().parents[1], time.monotonic()
    manifest = json.loads(args.manifest.read_text())
    names = ['tools/search_pascal_coordinate_bundle.py','tools/search_rectangular_core_bundle.py',
             'tools/rectangular_pascal_sectors.py','tools/permuted_pascal_kernels.py',
             'tools/pascal_coordinate_sectors.py','tools/native_matrix_smith_f4.py',
             'tools/native_matrix_smith_f4.cpp','tools/native_rectangular_smith_f4.cpp',
             'tools/native_weighted_smith_f4.py','tools/cube_f4_bivariate_smith.py',
             'tools/weighted_pascal_flags.py','tools/tensor_pascal_box_crossings.py']
    report = {'status':'running_incomplete','registration':manifest,
              'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names},'cases':[],
              'scope':'Timed numerical fixed-type convex bundle steps with Gibbs updates. '
                      'Complete original-polynomial profiles, no exact contraction/global optimum.'}
    def save():
        report['seconds'] = time.monotonic()-begun
        args.output.parent.mkdir(parents=True,exist_ok=True)
        temp = args.output.with_suffix('.tmp')
        temp.write_text(json.dumps(report,indent=2)+'\n')
        temp.replace(args.output)
    save()
    for case in manifest['cases']:
        start = time.monotonic()
        seed = permuted_pascal_kernel(case['images'])
        model = RectangularPascalPressure(seed,case['alphabet_b'],case['alphabet_x'],precision=case['precision'])
        dimension, cuts, best, minimum_slack = 2*model.size+1, [], None, 0.
        record = {'label':case['label'],'status':'running_incomplete','cuts':[],
                  'steps':[],'profile_rejections':0,'master_failures':0}
        report['cases'].append(record)
        search_start = None
        def query(theta):
            nonlocal best,minimum_slack
            if search_start is not None and time.monotonic()-search_start >= manifest['search_seconds_per_case']:
                raise TimeoutError
            result = log_profile_query(model,theta)
            rho,f,g,detail,physical = result
            if cuts:
                slacks = [f-(oldf+oldg@(theta-oldx)) for oldx,oldf,oldg in cuts]
                slacks += [oldf-(f+g@(oldx-theta)) for oldx,oldf,oldg in cuts]
                minimum_slack = min(minimum_slack,min(slacks))
                if min(slacks) < -2e-7:
                    raise ArithmeticError(f'Invalid convex supporting plane: {min(slacks)}.')
            cuts.append((theta.copy(),f,g.copy()))
            record['cuts'].append({'parameters':theta.tolist(),'log_profile_product':f,'subgradient':g.tolist()})
            if best is None or rho < best['rho']-1e-12:
                best = detail
                record['best'] = best
                print(json.dumps({'case':case['label'],'query':model.queries,'rho':rho,
                                  'seconds':time.monotonic()-start}),flush=True)
            record.update(queries=model.queries,minimum_heldout_convexity_slack=minimum_slack)
            return result
        theta = np.asarray(case['parameters'],dtype=float)
        constant = [p*model.size+o for p in range(2) for o in model.offsets]
        for p in range(2):
            for leg,o in enumerate(model.offsets):
                block = theta[p*model.size+o:p*model.size+o+len(model.leg_alphabets[leg])]
                block -= block[0]
        initial = query(theta)
        if abs(initial[0]-case['baseline_rho']) > 1e-10:
            raise ArithmeticError('The registered starting rate did not replay.')
        record['baseline'] = initial[3]
        record['construction_seconds'] = time.monotonic()-start
        normal = np.zeros((2*(dimension-1),dimension+1))
        for j in range(dimension-1):
            normal[2*j,j],normal[2*j,dimension-1] = 1,manifest['normalized_height_bound']
            normal[2*j+1,j],normal[2*j+1,dimension-1] = -1,manifest['normalized_height_bound']
        low = np.r_[np.full(dimension-1,-60.),.0001,-1e7]
        high = np.r_[np.full(dimension-1,60.),2.,1e7]
        low[constant] = high[constant] = 0.
        bounds = Bounds(low,high)
        metric = np.r_[np.ones(dimension-1),seed['outer']]
        rng = np.random.default_rng(manifest['seed'])
        search_start = time.monotonic()
        save()
        try:
            for restart in range(manifest['restarts_ceiling']):
                if restart:
                    theta = np.asarray(best['parameters'])
                    theta[:-1] += rng.normal(0,.025,dimension-1)
                    theta[-1] *= np.exp(rng.normal(0,.1))
                    theta[:-1] = np.clip(theta[:-1],-.9*manifest['normalized_height_bound']*theta[-1],
                                        .9*manifest['normalized_height_bound']*theta[-1])
                    theta[constant] = 0.
                rho,f,g,detail,physical = query(theta)
                anchor,entropy,mu,serious = theta.copy(),detail['entropy'],.03,0
                current = (f+physical@anchor)/(2*entropy)
                for iteration in range(manifest['iterations_ceiling']):
                    if time.monotonic()-search_start >= manifest['search_seconds_per_case']:
                        raise TimeoutError
                    A = np.asarray([np.r_[-cg,1.] for _,_,cg in cuts])
                    lower = np.asarray([cf-cg@cp for cp,cf,cg in cuts])
                    def objective(x):
                        delta = (x[:dimension]-anchor)*metric
                        return (x[dimension]+physical@x[:dimension])/(2*entropy)+.5*mu*(delta@delta)
                    def derivative(x):
                        return np.r_[physical/(2*entropy)+mu*(metric**2)*(x[:dimension]-anchor),1/(2*entropy)]
                    initial_master = np.r_[anchor,max(cf+cg@(anchor-cp) for cp,cf,cg in cuts)]
                    opt = minimize(objective,initial_master,jac=derivative,method='SLSQP',
                        constraints=[LinearConstraint(A,lower,np.inf),LinearConstraint(normal,0,np.inf)],
                        bounds=bounds,options={'ftol':1e-11,'maxiter':300})
                    candidate = opt.x[:dimension]
                    if (not np.all(np.isfinite(candidate)) or np.min(A@opt.x-lower)<-1e-6
                            or np.min(normal@opt.x)<-1e-6):
                        record['master_failures'] += 1
                        record['steps'].append({'restart':restart,'iteration':iteration,
                            'status':'master_unresolved','message':str(opt.message)})
                        save()
                        break
                    predicted = current-(opt.x[dimension]+physical@candidate)/(2*entropy)
                    if predicted < 2e-9 or max(abs(candidate-anchor)) < 1e-9:
                        result = query(anchor)
                        change = max(abs(result[4]-physical))
                        physical,entropy = result[4],result[3]['entropy']
                        current = (result[1]+physical@anchor)/(2*entropy)
                        record['steps'].append({'restart':restart,'iteration':iteration,
                            'status':'type_update_at_bundle_stop','type_change':float(change)})
                        save()
                        if change < 1e-7:
                            break
                        mu = .03
                        continue
                    try:
                        result = query(candidate)
                    except ArithmeticError as error:
                        if 'Invalid convex' in str(error):
                            raise
                        record['profile_rejections'] += 1
                        mu = min(100.,mu*4)
                        record['steps'].append({'restart':restart,'iteration':iteration,
                            'status':'profile_unresolved','error':str(error)})
                        save()
                        continue
                    value = (result[1]+physical@candidate)/(2*entropy)
                    actual = current-value
                    accepted = actual >= .2*predicted
                    if accepted:
                        anchor,current,serious = candidate,value,serious+1
                        mu = max(.0001,mu*.9)
                        if serious % 6 == 0:
                            physical,entropy = result[4],result[3]['entropy']
                            current = (result[1]+physical@anchor)/(2*entropy)
                    else:
                        mu = min(10.,mu*1.15)
                    record['steps'].append({'restart':restart,'iteration':iteration,'rho':result[0],
                        'fixed_type_rho':float(value),'predicted':float(predicted),'actual':float(actual),
                        'serious_step':bool(accepted),'proximal_weight':mu,
                        'qp_success':bool(opt.success),'qp_iterations':int(opt.nit)})
                    save()
                    if best['rho'] < 1/np.sqrt(2)-1e-9:
                        raise StopIteration
            record['status'] = 'complete_iteration_restart_ceiling'
        except TimeoutError:
            record['status'] = 'complete_time_budget'
        except StopIteration:
            record['status'] = 'candidate_requires_exact_validation'
        record.update(best=best,queries=model.queries,search_seconds=time.monotonic()-search_start,
                      total_seconds=time.monotonic()-start,minimum_heldout_convexity_slack=minimum_slack)
        save()
        print(json.dumps({'case':case['label'],'status':record['status'],'best':best['rho'],
                          'queries':model.queries,'search_seconds':record['search_seconds']}),flush=True)
    report['status'] = 'complete_registered_campaign'
    save()


if __name__ == '__main__':
    main()
