#!/usr/bin/env python3
"""Exact original-tensor and weighted-basis controls for all bulk orderings."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_rectangular_defect_support import ZERO,ONE,Z,Rat
from check_rectangular_defect_core_basis import residue_operator,left_car,independent_rows
from check_reflection_hot_leg_bulk_gauges import original_profiles
from reflection_hot_alternate_coordinates import full_images,check_support
from reflection_hot_all_bulk_gauges import AllBulkHotStatePressure,branch_coordinates
from reflection_hot_leg_bulk_gauges import BOUNDARIES
from reflection_hot_state_support import full_images as original_images


def bulk_for_branch(place,branch):
    left,right=BOUNDARIES[place];out=np.zeros(4);swap=branch>=2;pivot='T'if branch%2 else'I'
    actual_left=left[::-1]if swap else left;actual_right=right[::-1]if swap else right
    out[list(actual_left)]=(3,1)if pivot=='I'else(1,3);out[list(actual_right)]=(2,4)
    return out


def generator_change_control(interval,ramification,core):
    old=original_images(interval,ramification,core,True);new=full_images(interval,ramification,core)
    m=interval//2;c=ONE+Z**(2*ramification);s=Z**(4*ramification);b=ONE+Z**(2*interval*ramification)
    matrix=[[c*s**(i-(i&j))for j in range(m)]for i in range(m)]
    def multiply(vector,transpose=False):
        return [sum(((matrix[j][i]if transpose else matrix[i][j])*vector[j]for j in range(m)),ZERO)for i in range(m)]
    for name,(creation,annihilation)in old.items():
        first=multiply(creation[:m]);second=multiply(creation[m:])
        actual_creation=[a+d for a,d in zip(first,creation[m:])]+[a/b for a in second]
        first=multiply(annihilation[:m],True);second=multiply(annihilation[m:],True)
        actual_annihilation=[a/b for a in first]+[a+d for a,d in zip(second,annihilation[:m])]
        assert (actual_creation,actual_annihilation)==new[name]
    return {'interval':interval,'ramification':ramification,'core':core,'generator_identities':len(old)}


def weighted_minor_control(interval,core,place,branch):
    pivot='T'if branch%2 else'I';swap=branch>=2
    co=branch_coordinates(interval,('zero','infinity')[place],1,core,pivot,swap,512)
    b=bulk_for_branch(place,branch);data=co['data'];checked=0
    for side in ('left','right'):
        basis=co[side];low,high=data['bulk_block_legs'][side]
        mode=[int(b[low if i<interval//2 else high])for i in data['active_indices']];ops={}
        for name,pair in basis['regular_images'].items():
            leg='xyuv'.index(name.split('_')[1][0]);creator=name.endswith('_create')
            own=-int(b[leg])if creator else int(b[leg]);sign=1 if side=='left'else-1
            creation=[Rat(*x)*Z**(sign*(own+w))for x,w in zip(pair['creation'],mode)]
            annihilation=[Rat(*x)*Z**(sign*(own-w))for x,w in zip(pair['annihilation'],mode)]
            assert all(x.order()>=0 for x in creation+annihilation if x)
            ops[name]=residue_operator((creation,annihilation));checked+=1
        size=1<<len(data['active_indices']);projector=[1<<i for i in range(size)]
        for name in basis['names']:
            occupied=left_car(left_car(projector,ops[name+'_annihilate']),ops[name+'_create'])
            projector=[a^b for a,b in zip(projector,occupied)]
        rows,columns=basis['vacuum_rows'],basis['vacuum_columns']
        minor=[sum(((projector[i]>>j)&1)<<k for k,j in enumerate(columns))for i in rows]
        assert len(independent_rows(minor)[0])==len(rows)
    return {'interval':interval,'core':core,'place':place,'branch':branch,'weighted_generators':checked,'unit_minors':2}


def literal_control(interval,ramification,core,branches,seed):
    factor=({0:{0:1},3:{-1:2,1:2},15:{0:3,2:3}}if core==1 else{0:{0:1},85:{-1:2,1:2},255:{0:3,2:3}})
    bulk=np.asarray([bulk_for_branch(p,branches[p])for p in range(2)]);rng=np.random.default_rng(seed)
    positive=rng.integers(0,5,size=(2,4,core))/4.;offsets=np.zeros((2,4,1<<core))
    for state in range(1<<core):
        offsets[:,:,state]=-bulk*state.bit_count()
        for bit in range(core):
            if state>>bit&1:offsets[:,:,state]+=positive[:,:,bit]
        if state.bit_count()>1:offsets[:,:,state]+=.5
    ell=.5;expected,rho,clearing,count=original_profiles(interval,ramification,core,factor,offsets,bulk,ell)
    model=AllBulkHotStatePressure(interval,factor,ramification,core);actual,detail=model.evaluate(ell,offsets,bulk)
    for place in range(2):
        data=model.records[place][0]['data'];spectrum=np.asarray(detail['active_orders'][place])+ell*(data['scalar_order']+(1 if place==0 else-1)*clearing)
        for spectator in detail['spectator_orders'][place]:spectrum=np.r_[spectrum,spectrum+ell*spectator]
        assert np.array_equal(np.sort(spectrum),np.sort(expected[place])),(interval,core,branches,place)
    assert abs(actual-rho)<2e-12
    return {'interval':interval,'ramification':ramification,'core':core,'branches':branches,
            'complete_weighted_cut_profiles':2,'endpoint_coefficients':count,'pressure_error':abs(actual-rho)}


def derivative_control(branch):
    model=AllBulkHotStatePressure(4,{0:{0:1},3:{-1:2,1:2},15:{0:3,2:3}})
    bulk=np.asarray([bulk_for_branch(place,branch)for place in range(2)])+np.asarray([[.13,.27,.31,.43],[.21,.19,.37,.29]])
    offsets=np.zeros((2,4,2));offsets[:,:,1]=-.5*bulk+np.asarray([[.21,.29,.41,.53],[.31,.27,.43,.51]])
    ell,eps=.317,1e-6;_,(dh,db,ds),_=model.evaluate(ell,offsets,bulk,True);errors=[]
    for index in np.ndindex(bulk.shape):
        d=np.zeros_like(bulk);d[index]=eps
        value=(model.evaluate(ell,offsets,bulk+d)[0]-model.evaluate(ell,offsets,bulk-d)[0])/(2*eps);errors.append(abs(value-db[index]))
        h=np.zeros_like(offsets);h[index+(1,)]=eps
        value=(model.evaluate(ell,offsets+h,bulk)[0]-model.evaluate(ell,offsets-h,bulk)[0])/(2*eps);errors.append(abs(value-dh[index+(1,)]))
    value=(model.evaluate(ell+eps,offsets,bulk)[0]-model.evaluate(ell-eps,offsets,bulk)[0])/(2*eps);errors.append(abs(value-ds))
    assert max(errors)<3e-8
    return {'branch':branch,'derivatives':len(errors),'maximum_error':max(errors)}


def verify(extended=False):
    generators=[generator_change_control(M,P,k)for M,P,k in ((4,1,1),(8,1,1),(4,3,2),(16,1,2))]
    support=[check_support(M,P,k)for M,P,k in ((4,1,1),(8,3,1),(16,1,2),(64,1,1),(64,3,2))]
    weighted=[weighted_minor_control(M,k,place,branch)for M,k in ((16,1),(16,2),(64,1))for place,branch in product(range(2),range(4))]
    controls=[literal_control(4,1,1,branches,800+i)for i,branches in enumerate(product(range(4),repeat=2))]
    controls += [literal_control(8,1,1,(branch,(branch+1)%4),830+branch)for branch in range(4)]
    controls += [literal_control(4,3,1,(1,3),835),literal_control(4,1,2,(1,2),836),literal_control(4,1,2,(3,0),837)]
    if extended:controls += [literal_control(8,1,2,(1,3),838),literal_control(8,1,2,(2,1),839)]
    expanded=branch_coordinates(128,'zero',1,1,'T',False,128)['metadata']
    assert expanded['exact_polynomial_matrix']and expanded['profile_precision']>128
    root=Path(__file__).resolve().parent.parent
    files=['tools/reflection_hot_alternate_coordinates.py','tools/reflection_hot_all_bulk_gauges.py',
           'tools/check_reflection_hot_all_bulk_gauges.py','tools/reflection_hot_leg_bulk_gauges.py',
           'tools/reflection_hot_state_support.py','tools/reflection_hot_state_coordinates.py','research/reflection-hot-all-bulk-gauges.md']
    return {'status':'pass','extended':extended,'generator_change_controls':generators,'support_controls':support,
            'weighted_basis_controls':weighted,'original_tensor_controls':controls,
            'derivative_controls':[derivative_control(branch)for branch in range(4)],
            'source_expansion_control':expanded,
            'complete_weighted_cut_profiles':sum(r['complete_weighted_cut_profiles']for r in controls),
            'endpoint_coefficient_equalities':sum(r['endpoint_coefficients']for r in controls),
            'source_sha256':{f:sha256((root/f).read_bytes()).hexdigest()for f in files},
            'scope':'All nonnegative physical-leg bulk orderings, with exact hot-cold regularity constraints and complete source matrices; no improved exponent claim.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--extended',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=verify(args.extended)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ('status','complete_weighted_cut_profiles','endpoint_coefficient_equalities')}))


if __name__=='__main__':main()
