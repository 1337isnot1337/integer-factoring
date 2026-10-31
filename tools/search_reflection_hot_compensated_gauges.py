#!/usr/bin/env python3
"""Search the bounded cone of cold-compensated negative hot-bit heights."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from reflection_hot_compensated_gauges import CompensatedHotStatePressure,boundary_bulks,physical_bulks
from reflection_hot_core_profiles import profile_data,optimize_scale
from search_reflection_hot_core_filters import factor_from_row
from search_reflection_hot_state_gauges import EXCEPTIONAL


def candidate(interval,row,baseline,starts,maxfun,rng,asymmetric=False):
    factor=factor_from_row(row,1);model=CompensatedHotStatePressure(interval,factor)
    count=4 if asymmetric else 2
    initial=np.zeros(9+count);initial[0]=math.log(baseline['log_s'])
    points=[initial]
    bulk_starts=((1,0,0,1),(0,1,1,0),(1,1,1,1)) if asymmetric else ((1,1),(1,0),(0,1))
    for bulk in bulk_starts:
        if len(points)>=starts:
            break
        point=initial.copy();point[1:1+count]=bulk;points.append(point)
    while len(points)<starts:
        point=initial.copy();point[1:1+count]=rng.uniform(0,2,size=count);point[1+count:]=rng.uniform(0,4,size=8);points.append(point)
    best=None;runs=[]
    def objective(parameters):
        ell=math.exp(parameters[0]);bulk=np.maximum(parameters[1:1+count],0.)
        if asymmetric:
            bulk=bulk.reshape(2,2)
        expanded,_=boundary_bulks(bulk);physical_bulk=physical_bulks(expanded)
        slack=np.maximum(parameters[1+count:],0.).reshape(2,4)
        offsets=np.zeros((2,4,2));offsets[:,:,1]=slack-2*physical_bulk
        value,(gradient,db,ds),detail=model.evaluate(ell,offsets,bulk,True)
        nonlocal best
        if best is None or value<best['rho']:
            best={'rho':value,'log_s':ell,'bulk_heights':bulk.tolist(),'relative_hot_heights':offsets[:,:,1].tolist(),
                  'total_added_hot_heights':(offsets[:,:,1]+physical_bulk).tolist(),'slacks':slack.tolist(),
                  'active_orders':detail['active_orders'],'core_type_x':detail['type_x_core'],
                  'core_type_y':detail['type_y_core'],'cold_physical_energies':detail['cold_physical_energies']}
        if asymmetric:
            db=db.copy()
            for place,sides in enumerate((((0,1),(2,3)),((2,0),(3,1)))):
                for side,legs in enumerate(sides):
                    db[place,side]-=2*gradient[place,list(legs),1].sum()
        else:
            db=db-2*gradient[:,:,1].sum(axis=1)
        derivative=np.r_[ell*ds,db.reshape(-1),gradient[:,:,1].reshape(-1)]
        return value,derivative
    for point in points:
        run=minimize(objective,point,jac=True,method='L-BFGS-B',
                     bounds=[(math.log(baseline['log_s']/4),math.log(baseline['log_s']*3))]+[(0.,8.)]*count+[(0.,24.)]*8,
                     options={'maxfun':maxfun,'maxiter':maxfun,'maxls':20,'ftol':1e-13,'gtol':2e-9})
        runs.append({'rho':float(run.fun),'evaluations':int(run.nfev),'success':bool(run.success),'message':str(run.message)})
    return {'row':list(row),'factor':factor,'best':best,'starts':runs,'queries':model.queries,
            'local_factors':model.local_factors,'spectator_orders_before_bulk':model.records[0][0]['data']['spectator_orders'],
            'coefficient_metadata':model.records[0][1]}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--interval',type=int,default=16)
    parser.add_argument('--representatives',type=Path);parser.add_argument('--starts',type=int,default=4)
    parser.add_argument('--maxfun',type=int,default=120);parser.add_argument('--random',type=int,default=0)
    parser.add_argument('--asymmetric',action='store_true')
    parser.add_argument('--seed',type=int,default=133746);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic();rng=np.random.default_rng(args.seed)
    rows=[(0,)*8,(1,1,0,0,1,1,1,1)]+list(EXCEPTIONAL)
    if args.representatives:
        source=json.loads(args.representatives.read_text())
        rows=list(dict.fromkeys(rows+[tuple(record['representative']) for record in source['profile_classes']]))
    target=len(rows)+args.random
    while len(rows)<target:
        row=tuple(map(int,rng.integers(0,4,size=8)))
        if row not in rows:
            rows.append(row)
    baseline=optimize_scale(profile_data(args.interval,1,{0:{0:1}}))
    result={'status':'exact coefficient profiles with exploratory floating compensated-height optimization',
            'interval':args.interval,'gaussian_natural':baseline,'planned_candidates':len(rows),
            'independent_boundary_bulks':args.asymmetric,
            'bulk_height_bounds':[0,8],'slack_bounds':[0,24],'seed':args.seed,'complete':False,'records':[]}
    root=Path(__file__).resolve().parent.parent
    files=['tools/search_reflection_hot_compensated_gauges.py','tools/reflection_hot_compensated_gauges.py']
    def save():
        result['seconds']=time.monotonic()-start;result['completed_candidates']=len(result['records'])
        result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in files}
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for row in rows:
        record=candidate(args.interval,row,baseline,args.starts,args.maxfun,rng,args.asymmetric)
        result['records'].append(record);save()
        print(json.dumps({'completed':len(result['records']),'planned':len(rows),'row':row,'rho':record['best']['rho'],
                          'bulk':record['best']['bulk_heights'],'queries':record['queries'],'seconds':result['seconds']}),flush=True)
    result['complete']=True;result['best']=min(result['records'],key=lambda record:record['best']['rho']);save()


if __name__=='__main__':
    main()
