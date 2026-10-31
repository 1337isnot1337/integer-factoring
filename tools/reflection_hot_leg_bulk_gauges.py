#!/usr/bin/env python3
"""Reflection hot cores with a separate cold bulk height on each physical leg.

The retained cross-only coordinate matrices are unchanged. Their vacuum
coordinates receive the low/high diagonal-block bulk height, while their
physical hot labels receive their own leg height. The regularity cone is
checked explicitly; complete source coefficients come from the existing
compensated engine.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from reflection_dual_rail_gaussian import gauged_kernel,physical_bit_energies
from reflection_hot_compensated_gauges import CompensatedHotStatePressure


BOUNDARIES=(((0,1),(2,3)),((2,0),(3,1)))


def mode_legs(place,side):
    """Physical leg whose bulk is assigned to diagonal blocks A and B."""
    left,right=BOUNDARIES[place]
    return (left[1],left[0]) if side=='left' else right


@lru_cache(maxsize=32)
def removed_internal_partners(interval,core,place):
    m=interval//2;kernel=gauged_kernel(interval,('zero','infinity')[place])
    hot={m-1-j for j in range(core)};answer={}
    for boundary in BOUNDARIES[place]:
        for leg in boundary:
            for j in range(core):
                i=m-1-j;entries=[]
                for partner in boundary:
                    for cold in range(m):
                        if cold in hot:
                            continue
                        coefficient=kernel[leg*m+i][partner*m+cold]
                        if coefficient:
                            entries.append((partner,(coefficient&-coefficient).bit_length()-1))
                answer[leg,j]=entries
    return answer


def coordinate_leg_counts(coordinates,side,place):
    basis=coordinates[side];data=coordinates['data'];m=data['M']//2
    if basis['vacuum_rows']!=basis['vacuum_columns']:
        raise ArithmeticError('Separate leg bulks require the proved principal vacuum minor.')
    vacuum=basis['vacuum_rows' if side=='left' else 'vacuum_columns']
    labels=basis['labels'];size=len(next(iter(labels.values())))
    result=np.zeros((size,4),dtype=float)
    for letter,values in labels.items():
        result[:, 'xyuv'.index(letter)]=[int(value).bit_count() for value in values]
    low,high=mode_legs(place,side)
    for index in range(size):
        state=vacuum[index%len(vacuum)]
        for bit,original in enumerate(data['active_indices']):
            if state>>bit&1:
                result[index,low if original<m else high]+=1
    return result


class LegBulkHotStatePressure(CompensatedHotStatePressure):
    def validate_bulk_order(self,place,bulk):
        left,right=BOUNDARIES[place]
        if bulk[left[0]]<bulk[left[1]] or bulk[right[1]]<bulk[right[0]]:
            raise ValueError('The natural cross basis requires left-first >= left-second and right-second >= right-first.')

    def coordinate_counts(self,coordinates,side,place):
        return coordinate_leg_counts(coordinates,side,place)

    def block_legs(self,place,side):
        return mode_legs(place,side)

    def evaluate(self,log_s,offsets=None,bulk=None,gradient=False):
        if not math.isfinite(log_s) or log_s<=0:
            raise ValueError('Use a positive finite logarithmic scale.')
        offsets=np.zeros((2,4,self.q)) if offsets is None else np.asarray(offsets,dtype=float)
        bulk=np.zeros((2,4)) if bulk is None else np.asarray(bulk,dtype=float)
        if offsets.shape!=(2,4,self.q) or bulk.shape!=(2,4) or np.any(bulk<0):
            raise ValueError('Use nonnegative bulk heights of shape (2,4) and core heights of shape (2,4,q).')
        for place in range(2):
            self.validate_bulk_order(place,bulk[place])
            partners=removed_internal_partners(self.interval,self.core,place)
            for leg in range(4):
                for bit in range(self.core):
                    lower=min((self.ramification*order+bulk[place,leg]+bulk[place,other]
                               for other,order in partners[leg,bit]),default=math.inf)
                    for state in range(self.q):
                        if not state>>bit&1 and offsets[place,leg,state|(1<<bit)]-offsets[place,leg,state]<-lower:
                            raise ValueError('A removed hot-cold Gaussian factor is not regular at the requested heights.')
        letters=dict(zip('xyuv',range(4)))
        log_b=0.;ds=0.;derivative=np.zeros_like(offsets);db=np.zeros_like(bulk);profiles=[];spectator_profiles=[]
        for place,((coordinates,metadata),oracles) in enumerate(zip(self.records,self.oracles)):
            left,right,data=coordinates['left'],coordinates['right'],coordinates['data']
            rn=self.coordinate_counts(coordinates,'left',place);cn=self.coordinate_counts(coordinates,'right',place)
            rows=rn@bulk[place];columns=cn@bulk[place]
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
                        raise ArithmeticError('The profile exceeds its retained coefficient threshold.')
                orders+=result['orders'];ri += [rg[i] for i in result['rows']];ci += [cg[i] for i in result['cols']]
            orders=np.asarray(orders);lz=float(logsumexp(-orders));probabilities=np.exp(-orders-lz)
            log_b+=lz;ds-=float(probabilities@orders)/log_s
            db[place]-=log_s*(probabilities@rn[ri]+probabilities@cn[ci])
            for basis,indices in ((left,ri),(right,ci)):
                for letter,labels in basis['labels'].items():
                    derivative[place,letters[letter]]-=log_s*np.bincount(np.asarray(labels)[indices],weights=probabilities,minlength=self.q)
            active=set(data['active_indices']);spectators=[];features=[]
            low_left,high_left=self.block_legs(place,'left');low_right,high_right=self.block_legs(place,'right')
            for index in range(self.interval):
                if index in active:
                    continue
                feature=np.zeros(4)
                for leg in ((low_left,low_right) if index<self.interval//2 else (high_left,high_right)):
                    feature[leg]+=1
                spectators.append(self.ramification*(2+4*index)+float(feature@bulk[place]));features.append(feature)
            spectators=np.asarray(spectators);features=np.asarray(features).reshape((-1,4))
            occupancy=np.exp(-np.logaddexp(0,log_s*spectators))
            log_b+=float(np.logaddexp(0,-log_s*spectators).sum())
            ds-=float(spectators@occupancy);db[place]-=log_s*(occupancy@features)
            profiles.append(orders.tolist());spectator_profiles.append(spectators.tolist())
        natural=self.ramification*np.asarray(physical_bit_energies(self.interval),dtype=float)
        shift_x=float(bulk[:,(0,3)].sum());shift_y=float(bulk[:,(1,2)].sum())
        hot=[self.interval//2-1-j for j in range(self.core)];cold=[i for i in range(len(natural)) if i not in hot]
        particles=np.asarray([state.bit_count() for state in range(self.q)])
        core_energy=np.asarray([sum(natural[hot[j]] for j in range(self.core) if state>>j&1) for state in range(self.q)])
        gx=core_energy+shift_x*particles+offsets[:,(0,3),:].sum(axis=(0,1))
        gy=core_energy+shift_y*particles+offsets[:,(1,2),:].sum(axis=(0,1))
        cold_x=natural[cold]+shift_x;cold_y=natural[cold]+shift_y
        def equation(rate):
            physical=float(logsumexp(-log_s*gx/rate)+logsumexp(-log_s*gy/rate)
                           +np.logaddexp(0,-log_s*cold_x/rate).sum()+np.logaddexp(0,-log_s*cold_y/rate).sum())
            return log_b-rate*physical
        if equation(1e-9)<=0:
            raise ArithmeticError('A nonpositive pressure root requires separate review.')
        upper=4.
        while equation(upper)>0:
            upper*=2
        rate=brentq(equation,1e-9,upper,xtol=1e-13)
        px=np.exp(-log_s*gx/rate-logsumexp(-log_s*gx/rate));py=np.exp(-log_s*gy/rate-logsumexp(-log_s*gy/rate))
        cx=np.exp(-np.logaddexp(0,log_s*cold_x/rate));cy=np.exp(-np.logaddexp(0,log_s*cold_y/rate))
        entropy=-sum(float(p[p>0]@np.log(p[p>0])) for p in (px,py))
        entropy+=float(np.logaddexp(0,-log_s*cold_x/rate).sum()+cx@(log_s*cold_x/rate))
        entropy+=float(np.logaddexp(0,-log_s*cold_y/rate).sum()+cy@(log_s*cold_y/rate))
        self.queries+=1
        detail={'rho':rate/2,'log_s':log_s,'bulk_heights':bulk.tolist(),'active_orders':profiles,
                'spectator_orders':spectator_profiles,'log_profile_product':log_b,'type_x_core':px.tolist(),
                'type_y_core':py.tolist(),'entropy':entropy,'cold_x_energies':cold_x.tolist(),'cold_y_energies':cold_y.tolist()}
        if not gradient:
            return rate/2,detail
        if entropy<1e-10:
            raise ArithmeticError('The physical entropy is too small for stable derivatives.')
        for place in range(2):
            derivative[place]+=log_s*np.stack((px,py,py,px))
            db[place,(0,3)]+=log_s*float(px@particles+cx.sum())
            db[place,(1,2)]+=log_s*float(py@particles+cy.sum())
        ds+=float(px@gx+py@gy+cx@cold_x+cy@cold_y)
        return rate/2,(derivative/(2*entropy),db/(2*entropy),ds/(2*entropy)),detail
