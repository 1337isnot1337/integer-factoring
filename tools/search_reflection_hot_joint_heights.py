#!/usr/bin/env python3
"""Two-hot-mode defects under nonlinear joint-occupation height penalties."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from check_reflection_hot_state_gauges import original_weighted
from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_hot_block_profiles import normalize_factor, profile_data
from reflection_hot_core_profiles import optimize_scale
from reflection_hot_state_gauges import HotStatePressure


DEFAULT_SINGLE_DEFECTS = (
    (255,3,0,1),(15,1,-16,-15),(255,1,-2,-1),(255,3,0,2),
    (255,2,-2,0),(255,1,2,3),(15,1,-16,-14),(255,2,-4,-3),
    (255,2,-2,2),(255,3,2,4),
)


def fixed_pressure(interval,ramification,spectators,profiles,first,second,ell):
    natural = ramification*np.asarray(physical_bit_energies(interval))
    energies = np.asarray([0.,natural[-1],natural[-2],natural[-1]+natural[-2]+2*(first+second)])
    log_b = sum(float(logsumexp(-ell*np.asarray(part))) for part in profiles)
    log_b += 2*float(np.logaddexp(0,-ell*np.asarray(spectators)).sum())
    def equation(rate):
        return log_b-rate*(2*float(logsumexp(-ell*energies/rate))
                          +2*float(np.logaddexp(0,-ell*natural[:-2]/rate).sum()))
    upper = 4.
    while equation(upper)>0:
        upper *= 2
    return brentq(equation,1e-9,upper,xtol=1e-13)/2


def nonlinear_control():
    records = []
    for factor in ({0:{0:1}},{0:{0:1},255:{0:3,1:3}}):
        offsets = np.zeros((2,4,4));offsets[0,:,3]=2;offsets[1,:,3]=3
        model = HotStatePressure(8,factor,core=2,internal_regular=True)
        expected,rho,clearing,count = original_weighted(8,1,factor,2,None,offsets,.375)
        actual,detail = model.evaluate(.375,offsets)
        for place in range(2):
            data = model.records[place][0]['data']
            profile = np.asarray(detail['active_orders'][place])+.375*(data['scalar_order']+(1 if place==0 else -1)*clearing)
            assert np.array_equal(np.sort(profile),np.sort(expected[place]))
        assert abs(actual-rho)<1e-12
        separate = fixed_pressure(8,1,[],[np.asarray(part)/.375 for part in detail['active_orders']],2,3,.375)
        assert abs(separate-rho)<1e-12
        records.append({'factor':factor,'original_complete_weighted_profiles':2,'endpoint_coefficients':count,
                        'independent_pressure_error':abs(separate-rho)})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval',type=int,default=16)
    parser.add_argument('--source',type=Path,help='Optional alternative single-defect search record; the saved ten-row catalogue is built in.')
    parser.add_argument('--candidates',type=int,default=10)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args();start=time.monotonic()
    if args.source:
        source=json.loads(args.source.read_text())
        selected=[normalize_factor(record['factor'],2) for record in source['best_classes'][:args.candidates]]
    else:
        selected=[{0:{0:1},mask:{first:coefficient,second:coefficient}}
                  for mask,coefficient,first,second in DEFAULT_SINGLE_DEFECTS[:args.candidates]]
    factors=[{0:{0:1}}]+selected
    baseline=optimize_scale(profile_data(args.interval,1,{0:{0:1}},2))
    result={'status':'exact integer-height profiles; exploratory floating scale optimization',
            'interval':args.interval,'core':2,'gaussian_natural':baseline,
            'height_states':'All physical legs: h00=h01=h10=0, h11=t_place.',
            'height_grid':[0,1,2,4,8,16],'planned_candidates':len(factors),'complete':False,
            'original_nonlinear_height_controls':nonlinear_control(),'records':[]}
    root=Path(__file__).resolve().parent.parent
    files=['tools/search_reflection_hot_joint_heights.py','tools/reflection_hot_state_gauges.py',
           'tools/reflection_hot_state_coordinates.py','tools/reflection_hot_state_support.py']
    def save():
        result['seconds']=time.monotonic()-start
        result['completed_candidates']=len(result['records'])
        result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in files}
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for factor in factors:
        model=HotStatePressure(args.interval,factor,core=2,internal_regular=True)
        spectators=model.records[0][0]['data']['spectator_orders']
        records=[];best=None
        for first,second in product(result['height_grid'],repeat=2):
            offsets=np.zeros((2,4,4));offsets[0,:,3]=first;offsets[1,:,3]=second
            _,detail=model.evaluate(1.,offsets)
            assert all(order==int(order) for part in detail['active_orders'] for order in part)
            profiles=[list(map(int,part)) for part in detail['active_orders']]
            objective=lambda ell:fixed_pressure(args.interval,1,spectators,profiles,first,second,ell)
            optimized=minimize_scalar(objective,bounds=(baseline['log_s']/8,baseline['log_s']*8),method='bounded',
                                      options={'xatol':1e-12})
            record={'zero_joint_height':first,'infinity_joint_height':second,'rho':float(optimized.fun),
                    'log_s':float(optimized.x),'scale_evaluations':int(optimized.nfev),'scale_success':bool(optimized.success)}
            records.append(record)
            if best is None or record['rho']<best['rho']:
                best={**record,'active_profiles':profiles}
        result['records'].append({'factor':factor,'best':best,'grid_records':records,
                                   'spectator_orders':spectators,'coefficient_metadata':model.records[0][1]})
        save();print(json.dumps({'completed':len(result['records']),'planned':len(factors),'rho':best['rho'],
                                'zero_joint_height':best['zero_joint_height'],'infinity_joint_height':best['infinity_joint_height'],
                                'seconds':result['seconds']}),flush=True)
    result['complete']=True;save()


if __name__=='__main__':
    main()
