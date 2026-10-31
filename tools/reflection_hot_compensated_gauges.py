#!/usr/bin/env python3
"""Bounded reflection cores with cold-bit compensation for negative hot heights."""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from check_rectangular_f4_repair import fmul
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import laurent_kernel,physical_bit_energies
from reflection_hot_block_profiles import normalize_factor
from reflection_hot_state_coordinates import gaussian_coordinates
from reflection_hot_state_gauges import defect_block


def exterior_product(first,second):
    result={}
    for a,pa in first.items():
        for b,pb in second.items():
            if a&b:
                continue
            target=result.setdefault(a|b,{})
            for i,x in pa.items():
                for j,y in pb.items():
                    target[i+j]=target.get(i+j,0)^fmul(x,y)
                    if not target[i+j]:
                        del target[i+j]
    return {mask:poly for mask,poly in result.items() if poly}


@lru_cache(maxsize=32)
def retained_internal_factor(interval,place,ramification=1,core=1):
    kernel=laurent_kernel(interval);m=interval//2
    boundaries=((0,1),(2,3)) if place=='zero' else ((2,0),(3,1))
    factor={0:{0:1}}
    for legs in boundaries:
        selected=[(leg*core+j,leg*m+m-1-j) for leg in legs for j in range(core)]
        for a,(ma,ia) in enumerate(selected):
            for mb,ib in selected[a+1:]:
                scalar={ramification*degree:value for degree,value in kernel[ia][ib].items()}
                if scalar:
                    factor=exterior_product(factor,{0:{0:1},(1<<ma)|(1<<mb):scalar})
    return normalize_factor(factor,core)


def occupations(basis,side):
    vacuum=basis['vacuum_rows' if side=='left' else 'vacuum_columns']
    return np.asarray([state.bit_count()+position.bit_count()
                       for state in range(len(next(iter(basis['labels'].values())))//len(vacuum))
                       for position in vacuum],dtype=float)


def boundary_bulks(bulk):
    """Rows are cuts; columns are the output and input boundary bulks."""
    source=np.zeros(2) if bulk is None else np.asarray(bulk,dtype=float)
    if source.shape==(2,):
        return np.repeat(source[:,None],2,axis=1),True
    if source.shape==(2,2):
        return source,False
    raise ValueError('Bulk heights have shape (2,) for uniform cuts or (2,2) for separate boundaries.')


def physical_bulks(bulk):
    return np.asarray([[bulk[0,0],bulk[0,0],bulk[0,1],bulk[0,1]],
                       [bulk[1,0],bulk[1,1],bulk[1,0],bulk[1,1]]])


class CompensatedHotStatePressure:
    def __init__(self,interval,factor,ramification=1,core=1,profile_precision=512):
        self.interval,self.ramification,self.core,self.q=interval,ramification,core,1<<core
        self.factor=normalize_factor(factor,core)
        self.local_factors=[];self.records=[];self.oracles=[];self.queries=0
        for place in ('zero','infinity'):
            coordinates=gaussian_coordinates(interval,place,ramification,core,profile_precision,True)
            factor=exterior_product(self.factor,retained_internal_factor(interval,place,ramification,core))
            self.local_factors.append(factor)
            sides=[];metadata=[]
            for sector in range(2):
                matrix,shifts,md=defect_block(coordinates,factor,sector)
                sides.append(MatrixOracle(matrix,len(matrix),shifts,md['profile_precision']))
                metadata.append(md)
            self.records.append((coordinates,metadata));self.oracles.append(sides)

    def evaluate(self,log_s,offsets=None,bulk=None,gradient=False):
        if not math.isfinite(log_s) or log_s<=0:
            raise ValueError('Use a positive finite logarithmic scale.')
        offsets=np.zeros((2,4,self.q)) if offsets is None else np.asarray(offsets,dtype=float)
        bulk,uniform=boundary_bulks(bulk)
        physical_bulk=physical_bulks(bulk)
        if offsets.shape!=(2,4,self.q) or np.any(bulk<0):
            raise ValueError('Use nonnegative two-place bulk heights and core state heights of shape (2,4,q).')
        for place in range(2):
            for state in range(self.q):
                for bit in range(self.core):
                    if not state>>bit&1 and np.any(offsets[place,:,state|(1<<bit)]-offsets[place,:,state]<-2*physical_bulk[place]):
                        raise ValueError('Every hot occupation increment must be at least minus twice its cut bulk height.')
        letters=dict(zip('xyuv',range(4)))
        log_b=0.;ds=0.;derivative=np.zeros_like(offsets);db=np.zeros((2,2));profiles=[]
        for place,((coordinates,metadata),oracles) in enumerate(zip(self.records,self.oracles)):
            left,right,data=coordinates['left'],coordinates['right'],coordinates['data']
            rn,cn=occupations(left,'left'),occupations(right,'right')
            rows=bulk[place,0]*rn;columns=bulk[place,1]*cn
            for letter,labels in left['labels'].items():
                rows+=offsets[place,letters[letter]][labels]
            for letter,labels in right['labels'].items():
                columns+=offsets[place,letters[letter]][labels]
            orders=[];ri=[];ci=[]
            for sector,(oracle,md) in enumerate(zip(oracles,metadata)):
                rg,cg=coordinates['row_groups'][sector],coordinates['column_groups'][sector]
                while True:
                    try:
                        result=oracle.profile(log_s*rows[rg],log_s*columns[cg],scale=log_s)
                        break
                    except ArithmeticError as error:
                        if '(status 3)' not in str(error) or oracle.precision>=8192:
                            raise
                        oracle.precision*=2
                md['smith_precision']=oracle.precision
                if not md['exact_polynomial_matrix']:
                    known=log_s*(md['known_below']+float(rows[rg].min())+float(columns[cg].min()))
                    if max(result['orders'])>=known:
                        raise ArithmeticError('The compensated profile exceeds its retained coefficient threshold.')
                orders+=result['orders'];ri += [rg[i] for i in result['rows']];ci += [cg[i] for i in result['cols']]
            orders=np.asarray(orders);lz=float(logsumexp(-orders));probabilities=np.exp(-orders-lz)
            log_b+=lz;ds-=float(probabilities@orders)/log_s
            db[place,0]-=log_s*float(probabilities@rn[ri]);db[place,1]-=log_s*float(probabilities@cn[ci])
            for basis,indices in ((left,ri),(right,ci)):
                for letter,labels in basis['labels'].items():
                    derivative[place,letters[letter]]-=log_s*np.bincount(np.asarray(labels)[indices],weights=probabilities,minlength=self.q)
            spectators=np.asarray(data['spectator_orders'],dtype=float)+bulk[place].sum()
            occupancy=np.exp(-np.logaddexp(0,log_s*spectators))
            log_b+=float(np.logaddexp(0,-log_s*spectators).sum())
            ds-=float(spectators@occupancy);db[place]-=log_s*float(occupancy.sum())
            profiles.append(orders.tolist())
        natural=self.ramification*np.asarray(physical_bit_energies(self.interval),dtype=float)+bulk.sum()
        hot=[self.interval//2-1-j for j in range(self.core)]
        cold=[i for i in range(len(natural)) if i not in hot]
        core_energy=np.asarray([sum(natural[hot[j]] for j in range(self.core) if state>>j&1) for state in range(self.q)])
        gx=core_energy+offsets[:,(0,3),:].sum(axis=(0,1));gy=core_energy+offsets[:,(1,2),:].sum(axis=(0,1))
        cold_energy=natural[cold]
        def equation(rate):
            physical=float(logsumexp(-log_s*gx/rate)+logsumexp(-log_s*gy/rate)
                           +2*np.logaddexp(0,-log_s*cold_energy/rate).sum())
            return log_b-rate*physical
        if equation(1e-9)<=0:
            raise ArithmeticError('A nonpositive pressure root requires separate review.')
        upper=4.
        while equation(upper)>0:
            upper*=2
        rate=brentq(equation,1e-9,upper,xtol=1e-13)
        px=np.exp(-log_s*gx/rate-logsumexp(-log_s*gx/rate));py=np.exp(-log_s*gy/rate-logsumexp(-log_s*gy/rate))
        cold_occupancies=np.exp(-np.logaddexp(0,log_s*cold_energy/rate))
        entropy=-sum(float(p[p>0]@np.log(p[p>0])) for p in (px,py))
        entropy+=2*float(np.logaddexp(0,-log_s*cold_energy/rate).sum()+cold_occupancies@(log_s*cold_energy/rate))
        self.queries+=1
        detail={'rho':rate/2,'log_s':log_s,'bulk_heights':bulk.tolist(),'active_orders':profiles,
                'log_profile_product':log_b,'type_x_core':px.tolist(),'type_y_core':py.tolist(),'entropy':entropy,
                'cold_physical_energies':cold_energy.tolist()}
        if not gradient:
            return rate/2,detail
        if entropy<1e-10:
            raise ArithmeticError('The physical entropy is too small for stable derivatives.')
        for place in range(2):
            derivative[place]+=log_s*np.stack((px,py,py,px))
        particles=np.asarray([state.bit_count() for state in range(self.q)])
        db+=log_s*float((px+py)@particles+2*cold_occupancies.sum())
        ds+=float(px@gx+py@gy+2*cold_occupancies@cold_energy)
        if uniform:
            db=db.sum(axis=1)
        return rate/2,(derivative/(2*entropy),db/(2*entropy),ds/(2*entropy)),detail
