#!/usr/bin/env python3
"""Literal full-tensor checks for cold-compensated negative hot gauges."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from check_reflection_hot_core_profiles import polynomial_seed
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import logical_heights,physical_bit_energies
from reflection_hot_compensated_gauges import CompensatedHotStatePressure,retained_internal_factor,boundary_bulks,physical_bulks
from reflection_hot_state_gauges import HotStatePressure
from relabelled_gaussian_sums import flatten


def original_profiles(interval,ramification,core,factor,offsets,bulk,ell):
    values,clearing=polynomial_seed(interval,ramification,factor,core)
    m,q=interval//2,1<<(interval//2)
    core_labels=[sum(((state>>(m-1-j))&1)<<j for j in range(core)) for state in range(q)]
    bulk,_=boundary_bulks(bulk);physical_bulk=physical_bulks(bulk)
    profiles=[]
    for place,name in enumerate(('zero','infinity')):
        h=logical_heights(interval,name)
        heights=[[ramification*sum(h[leg*m+j] for j in range(m) if state>>j&1)
                  +physical_bulk[place,leg]*state.bit_count()+offsets[place,leg,core_labels[state]] for state in range(q)]
                 for leg in range(4)]
        first,second=(((0,1),(2,3)) if place==0 else ((0,2),(1,3)))
        rw=[heights[first[0]][a]+heights[first[1]][b] for a,b in product(range(q),repeat=2)]
        cw=[heights[second[0]][a]+heights[second[1]][b] for a,b in product(range(q),repeat=2)]
        profiles.append(MatrixOracle(flatten(values,m,bool(place)),q*q,precision=1024).profile(
            ell*np.asarray(rw),ell*np.asarray(cw),scale=ell,infinity=bool(place))['orders'])
    natural=ramification*np.asarray(physical_bit_energies(interval))+bulk.sum()
    gx=[];gy=[]
    for state in range(q):
        energy=sum(natural[j] for j in range(m) if state>>j&1);label=core_labels[state]
        gx.append(energy+offsets[:,(0,3),label].sum());gy.append(energy+offsets[:,(1,2),label].sum())
    log_b=sum(float(logsumexp(-np.asarray(part))) for part in profiles)
    def equation(rate):
        return log_b-rate*(logsumexp(-ell*np.asarray(gx)/rate)+logsumexp(-ell*np.asarray(gy)/rate))
    upper=4.
    while equation(upper)>0:
        upper*=2
    return profiles,brentq(equation,1e-9,upper,xtol=1e-13)/2,clearing,len(values)


def literal_control(interval,ramification,core,factor,seed,asymmetric=False):
    rng=np.random.default_rng(seed);bulk=np.asarray([[1.25,2.5],[.75,3.25]] if asymmetric else [1.25,2.5])
    expanded,_=boundary_bulks(bulk);physical_bulk=physical_bulks(expanded)
    positive=rng.integers(0,9,size=(2,4,core))/4.
    offsets=np.zeros((2,4,1<<core))
    for state in range(1<<core):
        offsets[:,:,state]=-2*physical_bulk*state.bit_count()
        for j in range(core):
            if state>>j&1:
                offsets[:,:,state]+=positive[:,:,j]
        if state.bit_count()>1:
            offsets[:,:,state]+=.75
    ell=.5
    expected,rho,clearing,count=original_profiles(interval,ramification,core,factor,offsets,bulk,ell)
    model=CompensatedHotStatePressure(interval,factor,ramification,core)
    actual,detail=model.evaluate(ell,offsets,bulk)
    for place in range(2):
        data=model.records[place][0]['data']
        spectrum=np.asarray(detail['active_orders'][place])+ell*(data['scalar_order']+(1 if place==0 else -1)*clearing)
        for spectator in data['spectator_orders']:
            spectrum=np.r_[spectrum,spectrum+ell*(spectator+expanded[place].sum())]
        assert np.array_equal(np.sort(spectrum),np.sort(expected[place])),(interval,ramification,core,place)
    assert abs(actual-rho)<2e-12
    return {'interval':interval,'ramification':ramification,'core':core,'factor':factor,'bulk':bulk.tolist(),
            'minimum_total_hot_height':float(np.min(offsets[:,:,1:]+physical_bulk[:,:,None])),
            'complete_weighted_cut_profiles':2,'endpoint_coefficients':count,'pressure_error':abs(actual-rho)}


def derivative_control(asymmetric=False):
    rng=np.random.default_rng(96381)
    factor={0:{0:1},3:{-1:2,1:2},15:{0:3,2:3}}
    model=CompensatedHotStatePressure(4,factor)
    bulk=np.array([[1.713,2.279],[1.389,2.625]] if asymmetric else [1.713,2.279]);offsets=np.zeros((2,4,2))
    expanded,_=boundary_bulks(bulk)
    offsets[:,:,1]=rng.uniform(.2,1.8,size=(2,4))-2*physical_bulks(expanded)
    ell,epsilon=.317,1e-6
    _,(gradient,db,ds),_=model.evaluate(ell,offsets,bulk,True)
    errors=[]
    for place,leg in product(range(2),range(4)):
        delta=np.zeros_like(offsets);delta[place,leg,1]=epsilon
        value=(model.evaluate(ell,offsets+delta,bulk)[0]-model.evaluate(ell,offsets-delta,bulk)[0])/(2*epsilon)
        errors.append(abs(value-gradient[place,leg,1]))
    for index in np.ndindex(bulk.shape):
        delta=np.zeros_like(bulk);delta[index]=epsilon
        value=(model.evaluate(ell,offsets,bulk+delta)[0]-model.evaluate(ell,offsets,bulk-delta)[0])/(2*epsilon)
        errors.append(abs(value-db[index]))
    value=(model.evaluate(ell+epsilon,offsets,bulk)[0]-model.evaluate(ell-epsilon,offsets,bulk)[0])/(2*epsilon)
    errors.append(abs(value-ds))
    assert max(errors)<3e-8
    return {'derivatives':len(errors),'maximum_error':max(errors)}


def verify(extended=False):
    assert retained_internal_factor(16,'zero')=={0:{0:1},3:{-4:1,-2:1},12:{0:1,2:1},15:{-4:1,0:1}}
    assert retained_internal_factor(16,'infinity')=={0:{0:1},5:{0:1,2:1},10:{0:1,2:1},15:{0:1,4:1}}
    gaussian={0:{0:1}};mixed={0:{0:1},3:{-1:3,1:3},10:{0:2,2:2}}
    records=[literal_control(M,P,1,f,300+i) for i,(M,P,f) in enumerate(
        ((4,1,gaussian),(4,1,mixed),(4,3,mixed),(8,1,gaussian),(8,1,mixed)))]
    if extended:
        records.append(literal_control(8,1,2,{0:{0:1},85:{-1:2,1:2},255:{0:3,2:3}},306))
    records += [literal_control(4,1,2,{0:{0:1},85:{-1:2,1:2},255:{0:3,2:3}},307,True),
                literal_control(8,1,1,mixed,308,True)]
    positive=[]
    for M in (4,8,16,32):
        offsets=np.zeros((2,4,2));offsets[:,:,1]=np.arange(8).reshape(2,4)/4
        a=CompensatedHotStatePressure(M,mixed).evaluate(.25,offsets)[0]
        b=HotStatePressure(M,mixed,internal_regular=True).evaluate(.25,offsets)[0]
        assert abs(a-b)<1e-12
        positive.append({'interval':M,'rho_difference':abs(a-b)})
    invalid=np.zeros((2,4,2));invalid[0,0,1]=-3
    try:
        CompensatedHotStatePressure(4,gaussian).evaluate(.3,invalid,[1,1])
    except ValueError:
        pass
    else:
        raise AssertionError('The compensated cone accepted an unsupported hot height.')
    root=Path(__file__).resolve().parent.parent
    files=['tools/reflection_hot_compensated_gauges.py','tools/check_reflection_hot_compensated_gauges.py',
           'research/reflection-hot-compensated-gauges.md']
    return {'status':'pass','extended':extended,'original_tensor_controls':records,
            'positive_cone_regressions':positive,'gradient_control':derivative_control(),
            'asymmetric_gradient_control':derivative_control(True),
            'complete_weighted_cut_profiles':sum(r['complete_weighted_cut_profiles'] for r in records),
            'endpoint_coefficient_equalities':sum(r['endpoint_coefficients'] for r in records),
            'source_sha256':{name:sha256((root/name).read_bytes()).hexdigest() for name in files},
            'scope':'Exact original tensor and weighted-cut controls for a compensated negative-height cone; no improved exponent claim.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--extended',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=verify(args.extended)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('status','complete_weighted_cut_profiles','endpoint_coefficient_equalities')}))


if __name__=='__main__':
    main()
