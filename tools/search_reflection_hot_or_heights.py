#!/usr/bin/env python3
"""Nonlinear compensated two-bit heights with a monotone OR-state slack."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.special import logsumexp

from check_reflection_hot_compensated_gauges import original_profiles
from reflection_hot_compensated_gauges import CompensatedHotStatePressure
from reflection_hot_block_profiles import profile_data
from reflection_hot_core_profiles import optimize_scale
from reflection_dual_rail_gaussian import physical_bit_energies
from search_reflection_hot_joint_heights import DEFAULT_SINGLE_DEFECTS


def heights(bulk,slack):
    result=np.zeros((2,4,4))
    for place in range(2):
        for state in range(1,4):
            result[place,:,state]=slack[place]-2*bulk[place]*state.bit_count()
    return result


def fixed_pressure(interval,spectators,profiles,bulk,slack,ell):
    original=np.asarray(physical_bit_energies(interval),dtype=float)
    total_b,total_u=sum(bulk),sum(slack)
    core=np.asarray([0.,original[-1],original[-2],original[-1]+original[-2]])
    core+=np.asarray([0.,-2*total_b+2*total_u,-2*total_b+2*total_u,-4*total_b+2*total_u])
    cold=original[:-2]+2*total_b
    log_b=sum(float(logsumexp(-ell*np.asarray(part))) for part in profiles)
    log_b+=sum(float(np.logaddexp(0,-ell*(np.asarray(spectators)+2*b)).sum()) for b in bulk)
    def equation(rate):
        return log_b-rate*(2*float(logsumexp(-ell*core/rate))+2*float(np.logaddexp(0,-ell*cold/rate).sum()))
    upper=4.
    while equation(upper)>0:
        upper*=2
    return brentq(equation,1e-9,upper,xtol=1e-13)/2


def control():
    factor={0:{0:1},255:{0:3,1:3}};bulk=[1.,2.];slack=[1.,2.];offsets=heights(bulk,slack)
    expected,rho,clearing,count=original_profiles(8,1,2,factor,offsets,bulk,.25)
    model=CompensatedHotStatePressure(8,factor,core=2);actual,detail=model.evaluate(.25,offsets,bulk)
    for place in range(2):
        data=model.records[place][0]['data']
        profile=np.asarray(detail['active_orders'][place])+.25*(data['scalar_order']+(1 if place==0 else -1)*clearing)
        assert np.array_equal(np.sort(profile),np.sort(expected[place]))
    separate=fixed_pressure(8,[],[np.asarray(part)/.25 for part in detail['active_orders']],bulk,slack,.25)
    assert abs(rho-actual)<1e-12 and abs(rho-separate)<1e-12
    return {'original_weighted_cut_profiles':2,'endpoint_coefficients':count,'pressure_error':abs(rho-separate)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--interval',type=int,default=64)
    parser.add_argument('--candidates',type=int,default=3);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic()
    factors=[{0:{0:1}}]+[{0:{0:1},mask:{first:coefficient,second:coefficient}}
                        for mask,coefficient,first,second in DEFAULT_SINGLE_DEFECTS[:args.candidates]]
    grid=list(dict.fromkeys((first,second,ratio*first,ratio*second)
                            for first,second,ratio in product((0,1,2,4),(0,1,2,4),(0,1,2))))
    baseline=optimize_scale(profile_data(args.interval,1,{0:{0:1}},2))
    result={'status':'exact compensated integer-height profiles; exploratory floating scale optimization',
            'interval':args.interval,'core':2,'gaussian_natural':baseline,'original_nonlinear_control':control(),
            'height_rule':'h(state)=u_cut*indicator(state!=0)-2*b_cut*popcount(state) on every physical leg',
            'height_grid':grid,'planned_candidates':len(factors),'complete':False,'records':[]}
    root=Path(__file__).resolve().parent.parent
    files=['tools/search_reflection_hot_or_heights.py','tools/reflection_hot_compensated_gauges.py']
    def save():
        result['seconds']=time.monotonic()-start;result['completed_candidates']=len(result['records'])
        result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in files}
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for factor in factors:
        model=CompensatedHotStatePressure(args.interval,factor,core=2,profile_precision=1024)
        spectators=model.records[0][0]['data']['spectator_orders'];records=[];best=None
        for first,second,u0,u1 in grid:
            bulk=[first,second];slack=[u0,u1];offsets=heights(bulk,slack)
            _,detail=model.evaluate(.125,offsets,bulk)
            profiles=[[int(order/.125) for order in part] for part in detail['active_orders']]
            assert all(int(order/.125)==order/.125 for part in detail['active_orders'] for order in part)
            objective=lambda ell:fixed_pressure(args.interval,spectators,profiles,bulk,slack,ell)
            optimized=minimize_scalar(objective,bounds=(baseline['log_s']/4,baseline['log_s']*4),method='bounded',options={'xatol':1e-12})
            record={'bulk':bulk,'slack':slack,'rho':float(optimized.fun),'log_s':float(optimized.x),
                    'scale_evaluations':int(optimized.nfev),'scale_success':bool(optimized.success)}
            records.append(record)
            if best is None or record['rho']<best['rho']:
                best={**record,'active_profiles':profiles}
        result['records'].append({'factor':factor,'best':best,'grid_records':records,
                                  'spectator_orders_before_bulk':spectators,'coefficient_metadata':model.records[0][1]})
        save();print(json.dumps({'completed':len(result['records']),'planned':len(factors),'rho':best['rho'],
                                'bulk':best['bulk'],'slack':best['slack'],'seconds':result['seconds']}),flush=True)
    result['complete']=True;save()


if __name__=='__main__':
    main()
