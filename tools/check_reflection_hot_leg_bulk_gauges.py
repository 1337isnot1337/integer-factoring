#!/usr/bin/env python3
"""Original physical cuts and derivatives for separate physical-leg bulks."""
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
from check_rectangular_defect_support import Rat,Z
from check_rectangular_defect_core_basis import residue_operator,left_car,independent_rows
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import logical_heights,physical_bit_energies
from reflection_hot_compensated_gauges import CompensatedHotStatePressure,physical_bulks
from reflection_hot_leg_bulk_gauges import LegBulkHotStatePressure,mode_legs
from reflection_hot_state_coordinates import coordinate_basis
from reflection_hot_state_support import core_data
from relabelled_gaussian_sums import flatten


def weighted_basis_control(interval,core,ramification=1):
    data=core_data(interval,'zero',ramification,core,True);m=interval//2
    checked=0
    for place in range(2):
        bulk=np.asarray([[3,1,2,5],[2,5,3,1]])[place]
        for side in ('left','right'):
            basis=coordinate_basis(interval,('zero','infinity')[place],side,ramification,core,True)
            assert basis['vacuum_rows']==basis['vacuum_columns']
            low,high=mode_legs(place,side)
            mode=[int(bulk[low if i<m else high]) for i in data['active_indices']]
            images={}
            for name,pair in basis['regular_images'].items():
                leg='xyuv'.index(name.split('_')[1][0]);creator=name.endswith('_create')
                a=[Rat(*x) for x in pair['creation']];b=[Rat(*x) for x in pair['annihilation']]
                sign=1 if side=='left' else -1
                own=-int(bulk[leg]) if creator else int(bulk[leg])
                a=[x*Z**(sign*(own+weight)) for x,weight in zip(a,mode)]
                b=[x*Z**(sign*(own-weight)) for x,weight in zip(b,mode)]
                assert all(x.order()>=0 for x in a+b if x)
                images[name]=residue_operator((a,b));checked+=1
            size=1<<len(data['active_indices']);projector=[1<<i for i in range(size)]
            for name in basis['names']:
                occupied=left_car(left_car(projector,images[name+'_annihilate']),images[name+'_create'])
                projector=[a^b for a,b in zip(projector,occupied)]
            selected=basis['vacuum_rows']
            minor=[sum(((projector[i]>>j)&1)<<k for k,j in enumerate(selected))for i in selected]
            assert len(independent_rows(minor)[0])==len(selected)
    return {'interval':interval,'core':core,'ramification':ramification,'weighted_generators':checked,'principal_unit_minors':4}


def original_profiles(interval,ramification,core,factor,offsets,bulk,ell):
    values,clearing=polynomial_seed(interval,ramification,factor,core)
    m,q=interval//2,1<<(interval//2)
    labels=[sum(((state>>(m-1-j))&1)<<j for j in range(core)) for state in range(q)]
    profiles=[]
    for place,name in enumerate(('zero','infinity')):
        h=logical_heights(interval,name)
        heights=[[ramification*sum(h[leg*m+j]for j in range(m)if state>>j&1)
                  +bulk[place,leg]*state.bit_count()+offsets[place,leg,labels[state]]for state in range(q)]for leg in range(4)]
        first,second=(((0,1),(2,3))if place==0 else((0,2),(1,3)))
        rows=[heights[first[0]][a]+heights[first[1]][b]for a,b in product(range(q),repeat=2)]
        columns=[heights[second[0]][a]+heights[second[1]][b]for a,b in product(range(q),repeat=2)]
        profiles.append(MatrixOracle(flatten(values,m,bool(place)),q*q,precision=1024).profile(
            ell*np.asarray(rows),ell*np.asarray(columns),scale=ell,infinity=bool(place))['orders'])
    natural=ramification*np.asarray(physical_bit_energies(interval));gx=[];gy=[]
    for state in range(q):
        energy=sum(natural[j]for j in range(m)if state>>j&1);label=labels[state]
        gx.append(energy+bulk[:,(0,3)].sum()*state.bit_count()+offsets[:,(0,3),label].sum())
        gy.append(energy+bulk[:,(1,2)].sum()*state.bit_count()+offsets[:,(1,2),label].sum())
    log_b=sum(float(logsumexp(-np.asarray(part)))for part in profiles)
    def equation(rate):
        return log_b-rate*(logsumexp(-ell*np.asarray(gx)/rate)+logsumexp(-ell*np.asarray(gy)/rate))
    upper=4.
    while equation(upper)>0:upper*=2
    return profiles,brentq(equation,1e-9,upper,xtol=1e-13)/2,clearing,len(values)


def literal_control(interval,ramification,core,factor,seed):
    rng=np.random.default_rng(seed);bulk=np.asarray([[3.,1.,2.,4.],[2.,4.,3.,1.]])
    positive=rng.integers(0,5,size=(2,4,core))/4;offsets=np.zeros((2,4,1<<core))
    for state in range(1<<core):
        offsets[:,:,state]=-bulk*state.bit_count()
        for bit in range(core):
            if state>>bit&1:offsets[:,:,state]+=positive[:,:,bit]
        if state.bit_count()>1:offsets[:,:,state]+=.5
    ell=.5;expected,rho,clearing,count=original_profiles(interval,ramification,core,factor,offsets,bulk,ell)
    model=LegBulkHotStatePressure(interval,factor,ramification,core);actual,detail=model.evaluate(ell,offsets,bulk)
    for place in range(2):
        data=model.records[place][0]['data']
        spectrum=np.asarray(detail['active_orders'][place])+ell*(data['scalar_order']+(1 if place==0 else -1)*clearing)
        for spectator in detail['spectator_orders'][place]:spectrum=np.r_[spectrum,spectrum+ell*spectator]
        assert np.array_equal(np.sort(spectrum),np.sort(expected[place])),(interval,ramification,core,place)
    assert abs(actual-rho)<2e-12
    return {'interval':interval,'ramification':ramification,'core':core,'factor':factor,'bulk':bulk.tolist(),
            'complete_weighted_cut_profiles':2,'endpoint_coefficients':count,'pressure_error':abs(actual-rho)}


def derivative_control():
    rng=np.random.default_rng(133792);model=LegBulkHotStatePressure(4,{0:{0:1},3:{-1:2,1:2},15:{0:3,2:3}})
    bulk=np.asarray([[3.713,1.279,2.389,4.625],[2.713,4.279,3.389,1.625]])
    offsets=np.zeros((2,4,2));offsets[:,:,1]=-.5*bulk+rng.uniform(.2,.8,size=(2,4))
    ell,epsilon=.317,1e-6;_,(dh,db,ds),_=model.evaluate(ell,offsets,bulk,True);errors=[]
    for index in np.ndindex(bulk.shape):
        delta=np.zeros_like(bulk);delta[index]=epsilon
        value=(model.evaluate(ell,offsets,bulk+delta)[0]-model.evaluate(ell,offsets,bulk-delta)[0])/(2*epsilon)
        errors.append(abs(value-db[index]))
        d=np.zeros_like(offsets);d[index+(1,)]=epsilon
        value=(model.evaluate(ell,offsets+d,bulk)[0]-model.evaluate(ell,offsets-d,bulk)[0])/(2*epsilon)
        errors.append(abs(value-dh[index+(1,)]))
    value=(model.evaluate(ell+epsilon,offsets,bulk)[0]-model.evaluate(ell-epsilon,offsets,bulk)[0])/(2*epsilon)
    errors.append(abs(value-ds));assert max(errors)<3e-8
    return {'derivatives':len(errors),'maximum_error':max(errors)}


def verify(extended=False):
    bases=[weighted_basis_control(M,k,P)for M,k,P in ((4,1,1),(4,2,1),(8,1,3),(16,2,1),(64,1,1),(64,2,3))]
    gaussian={0:{0:1}};mixed={0:{0:1},3:{-1:2,1:2},15:{0:3,2:3}}
    controls=[literal_control(M,P,k,f,610+i)for i,(M,P,k,f)in enumerate(
        ((4,1,1,gaussian),(4,1,1,mixed),(4,3,1,mixed),(8,1,1,mixed),
         (4,1,2,{0:{0:1},85:{-1:2,1:2},255:{0:3,2:3}})))]
    if extended:controls.append(literal_control(8,1,2,{0:{0:1},85:{-1:2,1:2},255:{0:3,2:3}},616))
    common=[]
    for M in (4,8,16,64):
        bulk=np.asarray([[1.25,2.5],[.75,3.25]]);physical=physical_bulks(bulk)
        h=np.zeros((2,4,2));h[:,:,1]=-physical+.25*np.arange(8).reshape(2,4)
        a=LegBulkHotStatePressure(M,mixed).evaluate(.25,h,physical)[0]
        b=CompensatedHotStatePressure(M,mixed).evaluate(.25,h,bulk)[0]
        assert abs(a-b)<2e-12;common.append({'interval':M,'pressure_error':abs(a-b)})
    root=Path(__file__).resolve().parent.parent
    files=['tools/reflection_hot_leg_bulk_gauges.py','tools/check_reflection_hot_leg_bulk_gauges.py','research/reflection-hot-leg-bulk-gauges.md']
    return {'status':'pass','extended':extended,'weighted_basis_controls':bases,'original_tensor_controls':controls,
            'common_boundary_regressions':common,'derivative_control':derivative_control(),
            'complete_weighted_cut_profiles':sum(r['complete_weighted_cut_profiles']for r in controls),
            'endpoint_coefficient_equalities':sum(r['endpoint_coefficients']for r in controls),
            'source_sha256':{f:sha256((root/f).read_bytes()).hexdigest()for f in files},
            'scope':'Separate nonnegative physical-leg bulks in an explicit regularity cone; no improved exponent claim.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--extended',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=verify(args.extended)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ('status','complete_weighted_cut_profiles','endpoint_coefficient_equalities')}))


if __name__=='__main__':main()
