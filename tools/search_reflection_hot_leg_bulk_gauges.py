#!/usr/bin/env python3
"""Explore separate physical-leg cold bulks in the proved ordering cone."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from reflection_hot_leg_bulk_gauges import LegBulkHotStatePressure,BOUNDARIES
from reflection_hot_core_profiles import profile_data,optimize_scale
from search_reflection_hot_core_filters import factor_from_row
from search_reflection_hot_state_gauges import EXCEPTIONAL


def bulk_and_offsets(low,excess,slack):
    bulk=np.zeros((2,4));offsets=np.zeros((2,4,2))
    for place,sides in enumerate(BOUNDARIES):
        for side,legs in enumerate(sides):
            big=legs[0]if side==0 else legs[1]
            bulk[place,list(legs)]=low[place,side]
            bulk[place,big]+=excess[place,side]
            offsets[place,list(legs),1]=slack[place,list(legs)]-2*low[place,side]
            offsets[place,big,1]=slack[place,big]-(bulk[place,big]+low[place,side])
    return bulk,offsets


def candidate(interval,row,baseline,starts,maxfun,rng):
    factor=factor_from_row(row,1);model=LegBulkHotStatePressure(interval,factor,profile_precision=max(512,8*interval))
    initial=np.zeros(17);initial[0]=math.log(baseline['log_s']);points=[initial]
    for low,excess in ((1,0),(0,1),(1,1)):
        if len(points)>=starts:break
        point=initial.copy();point[1:5]=low;point[5:9]=excess;points.append(point)
    while len(points)<starts:
        point=initial.copy();point[1:5]=rng.uniform(0,2,4);point[5:9]=rng.uniform(0,3,4);point[9:]=rng.uniform(0,4,8);points.append(point)
    best=None;runs=[]
    def objective(parameters):
        ell=math.exp(parameters[0]);low=np.maximum(parameters[1:5],0).reshape(2,2)
        excess=np.maximum(parameters[5:9],0).reshape(2,2);slack=np.maximum(parameters[9:],0).reshape(2,4)
        bulk,offsets=bulk_and_offsets(low,excess,slack)
        value,(dh,db,ds),detail=model.evaluate(ell,offsets,bulk,True)
        nonlocal best
        if best is None or value<best['rho']:
            best={'rho':value,'log_s':ell,'low_bulks':low.tolist(),'bulk_excesses':excess.tolist(),
                  'physical_leg_bulks':bulk.tolist(),'slacks':slack.tolist(),'relative_hot_heights':offsets[:,:,1].tolist(),
                  'total_added_hot_heights':(bulk+offsets[:,:,1]).tolist(),'active_orders':detail['active_orders'],
                  'spectator_orders':detail['spectator_orders'],'type_x_core':detail['type_x_core'],
                  'type_y_core':detail['type_y_core'],'cold_x_energies':detail['cold_x_energies'],'cold_y_energies':detail['cold_y_energies']}
        dl=np.zeros((2,2));de=np.zeros((2,2))
        for place,sides in enumerate(BOUNDARIES):
            for side,legs in enumerate(sides):
                big=legs[0]if side==0 else legs[1]
                dl[place,side]=db[place,list(legs)].sum()-2*dh[place,list(legs),1].sum()
                de[place,side]=db[place,big]-dh[place,big,1]
        return value,np.r_[ell*ds,dl.reshape(-1),de.reshape(-1),dh[:,:,1].reshape(-1)]
    for point in points:
        run=minimize(objective,point,jac=True,method='L-BFGS-B',
                     bounds=[(math.log(baseline['log_s']/4),math.log(baseline['log_s']*3))]+[(0.,8.)]*4+[(0.,16.)]*4+[(0.,24.)]*8,
                     options={'maxfun':maxfun,'maxiter':maxfun,'maxls':20,'ftol':1e-13,'gtol':2e-9})
        runs.append({'rho':float(run.fun),'evaluations':int(run.nfev),'success':bool(run.success),'message':str(run.message)})
    return {'row':list(row),'factor':factor,'best':best,'starts':runs,'queries':model.queries,
            'local_factors':model.local_factors,'coefficient_metadata':model.records[0][1]}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--interval',type=int,default=64)
    parser.add_argument('--representatives',type=Path);parser.add_argument('--starts',type=int,default=4)
    parser.add_argument('--maxfun',type=int,default=100);parser.add_argument('--random',type=int,default=16)
    parser.add_argument('--seed',type=int,default=133748);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic();rng=np.random.default_rng(args.seed)
    rows=[(0,)*8,(1,1,0,0,1,1,1,1)]+list(EXCEPTIONAL)
    if args.representatives:
        source=json.loads(args.representatives.read_text());rows=list(dict.fromkeys(rows+[tuple(r['representative'])for r in source['profile_classes']]))
    target=len(rows)+args.random
    while len(rows)<target:
        row=tuple(map(int,rng.integers(0,4,size=8)))
        if row not in rows:rows.append(row)
    baseline=optimize_scale(profile_data(args.interval,1,{0:{0:1}}))
    result={'status':'exact source profiles with exploratory floating physical-leg bulk optimization',
            'interval':args.interval,'gaussian_natural':baseline,'planned_candidates':len(rows),
            'low_bulk_bounds':[0,8],'bulk_excess_bounds':[0,16],'slack_bounds':[0,24],
            'seed':args.seed,'complete':False,'records':[]}
    root=Path(__file__).resolve().parent.parent
    files=['tools/search_reflection_hot_leg_bulk_gauges.py','tools/reflection_hot_leg_bulk_gauges.py']
    def save():
        result['seconds']=time.monotonic()-start;result['completed_candidates']=len(result['records'])
        result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest()for name in files}
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for row in rows:
        record=candidate(args.interval,row,baseline,args.starts,args.maxfun,rng);result['records'].append(record);save()
        print(json.dumps({'completed':len(result['records']),'planned':len(rows),'rho':record['best']['rho'],
                          'excess':record['best']['bulk_excesses'],'seconds':result['seconds']}),flush=True)
    result['complete']=True;result['best']=min(result['records'],key=lambda r:r['best']['rho']);save()


if __name__=='__main__':main()
