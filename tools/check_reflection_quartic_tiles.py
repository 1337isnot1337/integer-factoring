#!/usr/bin/env python3
"""Exact full-profile and typed-denominator controls for reflection quartic words."""
import argparse
from hashlib import sha256
from itertools import combinations,permutations,product
import json
from pathlib import Path
import random

from check_rectangular_defect_support import ONE,Z,pdivmod,pmul
from check_reflection_unrestricted_gaussian import from_polynomial,rational_smith
from reflection_dual_rail_gaussian import laurent_kernel
from reflection_gaussian_physical_frames import norm,add,mul
from reflection_quartic_tiles import (
    charge_vector,chi_power,complete_profiles,construct,field_multiply,
    field_power,full_rank_certificate,gate_data,
)
from physical_split_crossings import flatten


def polynomial_product(first,second):
    low,high,shift=mul(norm(*first,0),norm(*second,0))
    return low<<shift,high<<shift


def endpoint(value):
    return (value[0].bit_count()&1)+2*(value[1].bit_count()&1)


def edge_controls():
    checked=0
    for M in (4,8,16,32):
        charge=charge_vector(M);loss=M//2-1
        kernel=laurent_kernel(M)
        for i in range(2*M):
            for j in range(i+1,2*M):
                entry=kernel[i][j]
                if not entry:
                    continue
                shift=min(entry)
                value=sum(1<<(power-shift) for power in entry)
                order=0
                while True:
                    quotient,remainder=pdivmod(value,3)
                    if remainder:
                        break
                    order+=1;value=quotient
                assert charge[i]+charge[j]>=-loss*order
                checked+=1
    return checked


def two_copy_control(seed):
    q=seed['q_x']
    divisor=chi_power(2*seed['endpoint_pole_order'])
    count=0
    def at(x,y,u,v):
        return seed['values'][x+q*y+q*q*u+q**3*v]
    for x in product(range(q),repeat=2):
        for v in set(permutations(x)):
            for y in product(range(q),repeat=2):
                for u in set(permutations(y)):
                    numerator=polynomial_product(at(x[0],y[0],u[0],v[0]),
                                                 at(x[1],y[1],u[1],v[1]))
                    quotient=[]
                    for plane in numerator:
                        value,remainder=pdivmod(plane,divisor)
                        assert remainder==0
                        assert value.bit_length()-1<=2*seed['effective_cleared_degree']
                        quotient.append(value)
                    assert endpoint(quotient)==int(x==v and y==u)
                    count+=1
    return count


def exact_fock_profile(seed,side,heights):
    matrix,rows,columns=flatten(seed,side)
    rowlegs=(0,1) if side==0 else (0,2)
    collegs=(2,3) if side==0 else (1,3)
    rw=[sum(heights[leg][letter] for leg,letter in zip(rowlegs,pair)) for pair in rows]
    cw=[sum(heights[leg][letter] for leg,letter in zip(collegs,pair)) for pair in columns]
    clearing=seed['common_w_clearing'];pole=seed['endpoint_pole_order']
    scalar=Z**(clearing if side else -clearing)/(ONE+Z**(-1 if side else 1))**pole
    literal=[[from_polynomial(low,infinity=bool(side))*scalar*Z**(rw[i]+cw[j])
              for j,(low,high) in enumerate(row)] for i,row in enumerate(matrix)]
    assert all(high==0 for row in matrix for low,high in row)
    return rational_smith(literal)


def non_gaussian_control(seed,gate):
    modes=gate['modes'];zero=norm(*seed['values'][0],0)
    value=mul(zero,norm(*seed['values'][sum(1<<i for i in modes)],0))
    i,j,k,l=modes
    for first,second in (((i,j),(k,l)),((i,k),(j,l)),((i,l),(j,k))):
        a=norm(*seed['values'][sum(1<<v for v in first)],0)
        b=norm(*seed['values'][sum(1<<v for v in second)],0)
        value=add(value,mul(a,b))
    numerator=gate_data(seed['interval'],gate)[4]
    assert value==mul(numerator,mul(zero,zero)) and (value[0]|value[1])
    return list(value)


def verify():
    rng=random.Random(552861)
    omega=field_power(2,5)
    assert field_multiply(omega,omega)==omega^1
    for value in range(1,16):
        assert field_multiply(value,field_power(value,14))==1
    for a,b,c in product(range(16),repeat=3):
        assert field_multiply(field_multiply(a,b),c)==field_multiply(a,field_multiply(b,c))
    edges=edge_controls()
    cases=[
      (4,[{'modes':[0,2,4,6],'powers':[-2,2],'gamma':1}]),
      (4,[{'modes':[1,3,4,6],'kinds':['create','annihilate','create','annihilate'],'powers':[-4,1],'gamma':2}]),
      (4,[{'modes':[1,3,4,6],'powers':[3],'gamma':1}]),
      (4,[{'modes':[1,3,5,7],'kinds':['create','create','annihilate','annihilate'],'powers':[-2],'gamma':2}]),
      (4,[{'modes':[1,3,4,6],'powers':[3],'pole':1,'gamma':1}]),
      (4,[{'modes':[1,3,4,6],'powers':[4],'pole':1,'gamma':2}]),
      (4,[{'modes':[1,3,5,7],'kinds':['create','create','annihilate','annihilate'],'powers':[-1],'pole':3,'gamma':1}]),
      (4,[{'modes':[0,2,5,7],'kinds':['create','create','annihilate','annihilate'],'powers':[-2],'pole':1},
          {'modes':[0,1,3,4],'kinds':['annihilate','create','create','annihilate'],'powers':[3],'pole':1}]),
      (8,[{'modes':[0,1,2,3],'powers':[-3,2],'gamma':2}]),
      (8,[{'modes':[2,6,9,13],'powers':[-2]},{'modes':[3,7,8,12],'powers':[1]}]),
      (8,[{'modes':[2,6,10,14],'kinds':['create','create','annihilate','annihilate'],'powers':[-1],'pole':1},
          {'modes':[3,7,11,15],'kinds':['create','create','annihilate','annihilate'],'powers':[2],'pole':2}]),
    ]
    records=[];rational_profiles=typed_entries=non_gaussian=0
    for M,gates in cases:
        seed=construct(M,gates);q=seed['q_x']
        if seed['endpoint_pole_order']:
            seed['full_rank_witness']=full_rank_certificate(seed)
        heights=[[[rng.randrange(-8,9) for _ in range(q)] for leg in range(4)] for side in range(2)]
        profiles=complete_profiles(seed,heights)
        exact=[]
        if M==4:
            typed_entries+=two_copy_control(seed)
            if all(high==0 for low,high in seed['values']):
                for side in (0,1):
                    assert exact_fock_profile(seed,side,heights[side])==profiles[side]
                    exact.append(side);rational_profiles+=1
        residual=None
        if len(gates)==1 and gates[0].get('kinds',['create']*4)==['create']*4 and seed['strict_endpoint']:
            residual=non_gaussian_control(seed,gates[0]);non_gaussian+=1
        records.append({'M':M,'gates':gates,'strict_endpoint':seed['strict_endpoint'],
                        'type_projection_required':seed['type_projection_required'],
                        'degree':seed['degree'],'common_w_clearing':seed['common_w_clearing'],
                        'endpoint_pole_order':seed['endpoint_pole_order'],
                        'effective_cleared_degree':seed['effective_cleared_degree'],
                        'effective_degree_bound':seed['effective_degree_bound'],
                        'heights':heights,'restored_complete_profiles':profiles,
                        'independent_rational_profile_sides':exact,
                        'full_rank_witness':seed.get('full_rank_witness'),
                        'nonzero_quartic_residual':residual})
    rejected=0
    for gate in [
      {'modes':[1,3,4,6],'powers':[0],'pole':2},
      {'modes':[0,2,4,6],'powers':[0]},
    ]:
        try:construct(4,[gate])
        except ValueError:rejected+=1
        else:raise AssertionError('An unsupported endpoint condition was accepted.')
    return {'status':'Exact growing quartic and typed-pole controls passed',
            'gate_words':len(records),'complete_profiles':2*len(records),
            'independent_rational_fock_profiles':rational_profiles,
            'kernel_edge_charge_checks':edges,'two_copy_typed_entries':typed_entries,
            'nonzero_quartic_residuals':non_gaussian,'rejected_endpoint_conditions':rejected,
            'records':records}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);args=p.parse_args()
    result=verify();root=Path(__file__).resolve().parents[1]
    sources=['research/reflection-charged-quartic-tiles.md','tools/reflection_quartic_tiles.py',
             'tools/check_reflection_quartic_tiles.py','tools/reflection_dual_rail_gaussian.py',
             'tools/reflection_gaussian_physical_frames.py']
    result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in sources}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('records','source_sha256')},sort_keys=True),flush=True)


if __name__=='__main__':main()
