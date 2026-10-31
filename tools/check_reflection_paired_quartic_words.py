#!/usr/bin/env python3
"""Exact controls for growing high-pole paired quartic words."""
import argparse
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import random

from check_rectangular_defect_support import pdivmod
from check_reflection_quartic_tiles import (two_copy_control,exact_fock_profile,
    polynomial_product,endpoint)
from reflection_dual_rail_gaussian import laurent_kernel,tensor
from reflection_gaussian_physical_frames import norm,add
from reflection_quartic_tiles import complete_profiles,chi_power
from reflection_paired_quartic_words import construct,paired_gates,validate_word,chi_order
from weyl_lattice_filter_synthesis import polynomial_multiply


def gauged_kernel(M,orientation):
    m=M//2;positive,negative=orientation
    charge=[2 if i//m==positive else -2 if i//m==negative else 0 for i in range(2*M)]
    kernel=laurent_kernel(M);result=[]
    for i,row in enumerate(kernel):
        output=[]
        for j,entry in enumerate(row):
            if not entry:output.append((0,0,0));continue
            shift=min(entry);value=sum(1<<(power-shift) for power in entry)
            exponent=charge[i]+charge[j]
            if exponent>=0:value=polynomial_multiply(value,chi_power(exponent))
            else:
                value,remainder=pdivmod(value,chi_power(-exponent));assert not remainder
            output.append(norm(value,0,shift))
        result.append(output)
    return result


def kernel_controls():
    rows=intersections=0;records=[]
    for M in (4,8,16,32,64):
        m=M//2
        for orientation in ((0,3),(3,0),(1,2),(2,1)):
            kernel=gauged_kernel(M,orientation);positive,negative=orientation
            bit=1
            while bit<m:
                pairs=[(i,i^bit) for i in range(m) if not i&bit];gain=4*bit
                for i,j in pairs:
                    first=negative*m+i;second=negative*m+j
                    absent={positive*m+i,positive*m+j}
                    for column in range(2*M):
                        if column in absent:continue
                        value=add(kernel[first][column],kernel[second][column])
                        assert value[1]==0
                        if value[0]:assert chi_order(value[0])>=gain
                        rows+=1
                for (i,j),(k,l) in combinations(pairs,2):
                    value=(0,0,0)
                    for a,b in product((i,j),(k,l)):
                        value=add(value,kernel[negative*m+a][negative*m+b])
                    assert value==(0,0,0);intersections+=1
                records.append({'interval':M,'orientation':orientation,'xor_bit':bit,
                    'pairs':len(pairs),'pfaffian_gain_per_pair':gain,'pole_bound':gain+7})
                bit*=2
    return rows,intersections,records


def subword_coefficient_controls(M,xor_bit,orientation):
    m=M//2;positive,negative=orientation;source=tensor(M)['values']
    pairs=[(i,i^xor_bit) for i in range(m) if not i&xor_bit]
    checked=0
    for subset in range(1,1<<len(pairs)):
        selected=[pair for index,pair in enumerate(pairs) if subset>>index&1]
        occupied=sum(1<<(negative*m+i) for pair in selected for i in pair)
        empty=sum(1<<(positive*m+i) for pair in selected for i in pair)
        for state,(value,high) in enumerate(source):
            if not value or state&occupied!=occupied or state&empty:continue
            assert not high
            charge=2*sum(((state>>(positive*m+i))&1)-((state>>(negative*m+i))&1) for i in range(m))
            assert chi_order(value)+charge>=len(selected)*4*xor_bit
            checked+=1
    return checked


def sampled_two_copy(seed,rng,count=512):
    q=seed['q_x'];checked=0;divisor=chi_power(2*seed['endpoint_pole_order'])
    at=lambda x,y,u,v:seed['values'][x+q*y+q*q*u+q**3*v]
    for _ in range(count):
        x=[rng.randrange(q),rng.randrange(q)];y=[rng.randrange(q),rng.randrange(q)]
        v=x[:];u=y[:];rng.shuffle(v);rng.shuffle(u)
        value=polynomial_product(at(x[0],y[0],u[0],v[0]),at(x[1],y[1],u[1],v[1]))
        quotient=[]
        for plane in value:
            result,remainder=pdivmod(plane,divisor);assert not remainder
            assert result.bit_length()-1<=2*seed['effective_cleared_degree'];quotient.append(result)
        assert endpoint(quotient)==int(x==v and y==u);checked+=1
    return checked


def gaussian_residual(seed):
    values=seed['values'];vacuum=values[0]
    for modes in combinations(range(2*seed['interval']),4):
        i,j,k,l=modes
        residual=polynomial_product(vacuum,values[sum(1<<v for v in modes)])
        for a,b,c,d in ((i,j,k,l),(i,k,j,l),(i,l,j,k)):
            other=polynomial_product(values[(1<<a)|(1<<b)],values[(1<<c)|(1<<d)])
            residual=tuple(x^y for x,y in zip(residual,other))
        if any(residual):return {'modes':modes,'residual':residual}
    return None


def verify():
    rng=random.Random(938173);row_checks,intersections,kernels=kernel_controls()
    subwords=0
    for M,bit,orientation in ((4,1,(0,3)),(4,1,(2,1)),(8,2,(0,3)),(8,2,(1,2))):
        subwords+=subword_coefficient_controls(M,bit,orientation)
    cases=[(4,paired_gates(4,1,power=0,gamma=1)),
           (4,paired_gates(4,1,power=-5,gamma=2)),
           (4,paired_gates(4,1,orientation=(1,2),power=2,gamma=3)),
           (4,paired_gates(4,1,orientation=(3,0),power=-2,gamma=1)),
           (8,paired_gates(8,2,power=0,gamma=1)),
           (8,paired_gates(8,2,orientation=(2,1),power=-3,gamma=2)),
           (8,paired_gates(8,1,power=1,gamma=1))]
    records=[];profiles=independent=typed=residuals=0
    for M,gates in cases:
        seed=construct(M,gates);q=seed['q_x']
        heights=[[[rng.randrange(-5,6) for _ in range(q)] for _ in range(4)] for side in (0,1)]
        actual=complete_profiles(seed,heights);profiles+=2
        exact=[]
        if M==4:
            typed+=two_copy_control(seed)
            if all(not high for low,high in seed['values']):
                for side in (0,1):
                    assert exact_fock_profile(seed,side,heights[side])==actual[side]
                    independent+=1;exact.append(side)
            residual=gaussian_residual(seed)
            assert residual is not None;residuals+=1
        else:
            typed+=sampled_two_copy(seed,rng)
            residual=None
        assert seed['endpoint_pole_order']<=M-1
        assert seed['cancelled_common_chi_order']>=seed['raw_endpoint_pole_order']-(M-1)
        records.append({'interval':M,'gates':gates,'raw_pole':seed['raw_endpoint_pole_order'],
            'cancelled_chi_content':seed['cancelled_common_chi_order'],'remaining_pole':seed['endpoint_pole_order'],
            'degree':seed['degree'],'effective_degree':seed['effective_cleared_degree'],
            'effective_degree_bound':seed['effective_degree_bound'],'profiles':actual,
            'independent_rational_sides':exact,'non_gaussian_residual':residual})
    rejected=0
    for gates in (paired_gates(4,1,pole=12),
                  paired_gates(8,2,indices=[0,0]),
                  paired_gates(8,1,indices=[0])+paired_gates(8,2,indices=[1])):
        try:validate_word(4 if len(gates)==1 else 8,gates)
        except ValueError:rejected+=1
        else:raise AssertionError('An unproved word was accepted.')
    # The first excluded pole can fail at a balanced three-copy endpoint.
    source=tensor(4)['values'];raw=[(polynomial_multiply(a,chi_power(12)),0) for a,b in source]
    mask=(1<<0)|(1<<1)|(1<<6)|(1<<7);occupied=(1<<6)|(1<<7)
    for state,(value,_) in enumerate(source):
        if state&mask==occupied:
            a,b=raw[state^mask];raw[state^mask]=(a^value,b)
    at=lambda x,y,u,v:raw[x+4*y+16*u+64*v]
    value=polynomial_product(polynomial_product(at(3,0,0,0),at(0,0,1,1)),at(1,1,0,3))
    quotient=[]
    for plane in value:
        result,remainder=pdivmod(plane,chi_power(36));assert not remainder;quotient.append(result)
    assert endpoint(quotient)==1
    assert sorted((3,0,1))==sorted((0,1,3)) and sorted((0,0,1))==sorted((0,1,0))
    return {'status':'Exact amplified-pole paired quartic controls passed',
        'kernel_row_congruences':row_checks,'exact_zero_special_intersections':intersections,
        'original_subword_coefficient_margins':subwords,'complete_fock_profiles':profiles,
        'independent_rational_profiles':independent,'typed_two_copy_controls':typed,
        'nonzero_gaussian_residuals':residuals,'rejected_invalid_words':rejected,
        'sharp_boundary_three_copy_failures':1,'kernel_cases':kernels,'word_cases':records}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('evidence/reflection-paired-quartic-words.json'))
    args=parser.parse_args();result=verify()
    names=[Path(__file__).name,'reflection_paired_quartic_words.py','reflection_quartic_tiles.py',
        'check_reflection_quartic_tiles.py','reflection_dual_rail_gaussian.py',
        'reflection_gaussian_physical_frames.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp']
    result['source_sha256']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    note=Path(__file__).resolve().parents[1]/'research'/'reflection-paired-quartic-poles.md'
    if note.exists():result['note_sha256']=sha256(note.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in
        ('kernel_cases','word_cases','source_sha256','note_sha256')}))


if __name__=='__main__':main()
