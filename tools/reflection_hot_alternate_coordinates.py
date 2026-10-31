#!/usr/bin/env python3
"""The regular T-pivot reflection cross basis and its bounded hot coordinates."""
from __future__ import annotations

from functools import lru_cache

from check_rectangular_defect_support import ZERO,ONE,Z,car_checks,pairing
from reflection_hot_block_profiles import validate_core
from reflection_hot_state_support import concentrate
from reflection_hot_state_coordinates import coordinate_basis_for_data,gaussian_coordinates_for_data


def full_images(interval,ramification=1,core=1):
    validate_core(interval,core)
    if not isinstance(ramification,int)or ramification<1:
        raise ValueError('Use a positive integer substitution degree.')
    m,p=interval//2,ramification
    c,s,u=ONE+Z**(2*p),Z**(4*p),Z**(2*p)
    f=Z**(2*interval*p);b=ONE+f;answer={}
    def vector():return [ZERO]*interval
    for j in range(core):
        d=u*s**j
        column=[c*s**(r-(j&r))for r in range(m)]
        transpose=[c*s**(j-(j&r))for r in range(m)]
        creation,annihilation=vector(),vector();creation[j]=ONE
        for r in range(m):creation[m+r]=column[r]/b
        answer[f'left_first{j}_create']=(creation,annihilation)
        creation,annihilation=vector(),vector();annihilation[j]=ONE
        answer[f'left_first{j}_annihilate']=(creation,annihilation)
        creation,annihilation=vector(),vector();creation[m+j]=ONE
        answer[f'left_second{j}_create']=(creation,annihilation)
        creation,annihilation=vector(),vector();annihilation[m+j]=ONE
        for r in range(m):annihilation[r]=transpose[r]/b
        answer[f'left_second{j}_annihilate']=(creation,annihilation)
        creation,annihilation=vector(),vector()
        for r in range(m):creation[r]=d*column[r]
        answer[f'right_first{j}_create']=(creation,annihilation)
        creation,annihilation=vector(),vector();annihilation[m+j]=ONE/(f*d)
        for r in range(m):annihilation[r]=transpose[r]/(d*b)
        answer[f'right_first{j}_annihilate']=(creation,annihilation)
        creation,annihilation=vector(),vector();creation[j]=d
        for r in range(m):creation[m+r]=f*d*column[r]/b
        answer[f'right_second{j}_create']=(creation,annihilation)
        creation,annihilation=vector(),vector()
        for r in range(m):annihilation[m+r]=transpose[r]/(f*d)
        answer[f'right_second{j}_annihilate']=(creation,annihilation)
    return answer


@lru_cache(maxsize=32)
def core_data(interval,place='zero',ramification=1,core=1):
    if place not in ('zero','infinity'):raise ValueError('Choose zero or infinity.')
    before=full_images(interval,ramification,core)
    _,images,active,orders,trace=concentrate(interval,ramification,core,True,before)
    left=('x','y')if place=='zero'else('u','x');right=('u','v')if place=='zero'else('v','y')
    names={};renamed={}
    for side,letters in (('left',left),('right',right)):
        names[side]=[]
        for which,letter in zip(('first','second'),letters):
            for j in range(core):
                prefix=f'{side}_{letter}{j}';names[side].append(prefix)
                for kind in ('create','annihilate'):
                    creation,annihilation=images[f'{side}_{which}{j}_{kind}']
                    renamed[prefix+'_'+kind]={'creation':[creation[i].pair()for i in active],
                                             'annihilation':[annihilation[i].pair()for i in active]}
    return {'M':interval,'place':place,'ramification':ramification,'core':core,
            'cross_pivot':'T','internal_gaussians_removed':True,
            'active_indices':active,'active_orders':[orders[i]for i in active],
            'spectator_orders':[order for i,order in enumerate(orders)if i not in active],
            'active_images':renamed,'core_names':names,'regular_shears':trace,
            'logical_core_indices':[interval//2-1-j for j in range(core)],
            'scalar_order':(1 if place=='zero'else-1)*ramification*interval*interval//2}


@lru_cache(maxsize=32)
def coordinate_basis(interval,place,side,ramification=1,core=1):
    return coordinate_basis_for_data(core_data(interval,place,ramification,core),side)


@lru_cache(maxsize=8)
def gaussian_coordinates(interval,place,ramification=1,core=1,profile_precision=512):
    data=core_data(interval,place,ramification,core)
    if len(data['active_indices'])>11:raise ValueError('The parity exporter supports at most eleven active modes.')
    left,right=[coordinate_basis(interval,place,side,ramification,core)for side in ('left','right')]
    if place=='infinity':
        zero=gaussian_coordinates(interval,'zero',ramification,core,profile_precision)
        for current,original in ((left,zero['left']),(right,zero['right'])):
            assert current['vacuum_rows']==original['vacuum_rows']
            assert current['vacuum_columns']==original['vacuum_columns']
            for a,b in zip(current['names'],original['names']):
                for kind in ('create','annihilate'):
                    assert current['regular_images'][a+'_'+kind]==original['regular_images'][b+'_'+kind]
        return {**zero,'data':data,'left':left,'right':right,
                'metadata':{**zero['metadata'],'reused_identical_zero_matrix':True}}
    return gaussian_coordinates_for_data(data,left,right,profile_precision)


def check_support(interval,ramification=1,core=1):
    before=full_images(interval,ramification,core)
    _,after,active,orders,trace=concentrate(interval,ramification,core,True,before)
    canonical=0
    for side in ('left','right'):
        names=[f'{side}_{which}{j}'for which in ('first','second')for j in range(core)]
        canonical+=car_checks(before,names)+car_checks(after,names)
    for first,a in before.items():
        for second,b in before.items():assert pairing(a,b)==pairing(after[first],after[second])
    for name,(creation,annihilation)in after.items():
        if name.startswith('right'):
            creation=[value*Z**(-order)for value,order in zip(creation,orders)]
            annihilation=[value*Z**order for value,order in zip(annihilation,orders)]
        assert all(value.order()>=0 for value in creation+annihilation if value)
    return {'interval':interval,'ramification':ramification,'core':core,'active_modes':len(active),
            'spectator_modes':interval-len(active),'regular_shears':len(trace),
            'canonical_pairings':canonical,'all_cross_pairings':len(before)**2}
