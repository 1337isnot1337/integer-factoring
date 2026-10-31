#!/usr/bin/env python3
"""Exact original-coefficient and endpoint controls for mixed paired words."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_reflection_quartic_tiles import two_copy_control,exact_fock_profile
from check_reflection_paired_quartic_words import sampled_two_copy,gauged_kernel
from reflection_paired_quartic_words import paired_gates,chi_order,construct as one_direction
from reflection_mixed_paired_words import construct,validate_mixed_word
from reflection_dual_rail_gaussian import tensor
from reflection_gaussian_physical_frames import norm,add
from reflection_quartic_tiles import chi_power,complete_profiles
from check_rectangular_defect_support import pdivmod


def mixed_gates(M,bit,reverse_indices,base_orientation=(0,3),gamma=1):
    gates=[];m=M//2
    for i in range(m):
        if i&bit:continue
        reverse=i in reverse_indices
        orientation=base_orientation[::-1] if reverse else base_orientation
        gates+=paired_gates(M,bit,orientation=orientation,indices=[i],
            pole=4*bit+(-1 if reverse else 7),power=0,gamma=gamma)
    return gates


def kernel_control(M,bit,reverse_indices,orientation):
    m=M//2;positive,negative=orientation;kernel=gauged_kernel(M,orientation)
    pairs=[(i,i^bit) for i in range(m) if not i&bit]
    absent=set();scaled=set();special=[]
    for i,j in pairs:
        reverse=i in reverse_indices
        source=positive if reverse else negative
        missing=negative if reverse else positive
        absent.update((missing*m+i,missing*m+j));special.append((source*m+i,source*m+j))
        if reverse:scaled.update((source*m+i,source*m+j))
    active=[i for i in range(2*M) if i not in absent];matrix={};divisions=rows=intersections=0
    for i in active:
        for j in active:
            value=kernel[i][j];power=4*(int(i in scaled)+int(j in scaled))
            if value[0] and power:
                quotient,remainder=pdivmod(value[0],chi_power(power));assert not remainder
                value=norm(quotient,0,value[2]);divisions+=1
            matrix[i,j]=value
    for first,second in special:
        for col in active:
            value=add(matrix[first,col],matrix[second,col])
            if value[0]:assert chi_order(value[0])>=4*bit
            rows+=1
    for a,(i,j) in enumerate(special):
        for k,l in special[a+1:]:
            value=(0,0,0)
            for u,v in product((i,j),(k,l)):value=add(value,matrix[u,v])
            assert value==(0,0,0);intersections+=1
    return divisions,rows,intersections


def source_controls(M,bit,reverse_indices,orientation):
    m=M//2;positive,negative=orientation;source=tensor(M)['values']
    pairs=[(i,i^bit) for i in range(m) if not i&bit];checked=0
    for subset in range(1,1<<len(pairs)):
        occupied=empty=reverse_count=pair_count=0
        for index,(i,j) in enumerate(pairs):
            if not subset>>index&1:continue
            reverse=i in reverse_indices;pair_count+=1;reverse_count+=reverse
            leg=positive if reverse else negative;opposite=negative if reverse else positive
            occupied|=(1<<(leg*m+i))|(1<<(leg*m+j))
            empty|=(1<<(opposite*m+i))|(1<<(opposite*m+j))
        for state,(value,high) in enumerate(source):
            if not value or state&occupied!=occupied or state&empty:continue
            assert not high
            energy=2*sum(((state>>(positive*m+i))&1)-((state>>(negative*m+i))&1) for i in range(m))
            assert chi_order(value)+energy>=8*reverse_count+4*bit*pair_count;checked+=1
    return checked


def verify():
    rng=random.Random(972398);divisions=rows=intersections=source_checks=0
    kernel_records=[]
    for M in (4,8,16,32,64):
        m=M//2;bit=m//2;indices=[i for i in range(m) if not i&bit]
        reverse={i for j,i in enumerate(indices) if j%2}
        for orientation in ((0,3),(3,0),(1,2),(2,1)):
            a,b,c=kernel_control(M,bit,reverse,orientation)
            divisions+=a;rows+=b;intersections+=c
            kernel_records.append({'interval':M,'xor_bit':bit,'base_orientation':orientation,
                'reverse_pair_starts':sorted(reverse)})
            if M<=8:source_checks+=source_controls(M,bit,reverse,orientation)
    cases=[(4,paired_gates(4,1,orientation=(3,0),pole=3,power=0),(0,3)),
           (4,paired_gates(4,1,orientation=(2,1),pole=3,power=-2),(1,2)),
           (8,mixed_gates(8,2,{1},(0,3),1),(0,3)),
           (8,mixed_gates(8,2,{0},(1,2),2),(1,2)),
           (8,mixed_gates(8,1,{2},(3,0),3),(3,0))]
    records=[];profiles=typed=independent=0
    for M,gates,orientation in cases:
        seed=construct(M,gates,orientation);q=seed['q_x']
        heights=[[[rng.randrange(-4,5) for _ in range(q)] for _ in range(4)] for side in (0,1)]
        actual=complete_profiles(seed,heights);profiles+=2
        exact=[]
        if M==4:
            comparison=one_direction(M,gates)
            assert seed['values']==comparison['values']
            assert seed['endpoint_pole_order']==comparison['endpoint_pole_order']
            typed+=two_copy_control(seed)
            for side in (0,1):
                assert exact_fock_profile(seed,side,heights[side])==actual[side]
                independent+=1;exact.append(side)
        else:typed+=sampled_two_copy(seed,rng)
        records.append({'interval':M,'gates':gates,'base_orientation':orientation,
            'directions':seed['directions'],'per_gate_margins':seed['per_gate_margins'],
            'raw_pole':seed['raw_endpoint_pole_order'],'remaining_pole':seed['endpoint_pole_order'],
            'cancelled_chi':seed['cancelled_common_chi_order'],'effective_degree':seed['effective_cleared_degree'],
            'profiles':actual,'independent_rational_sides':exact})
    rejected=0
    invalid=[paired_gates(4,1)+paired_gates(4,1,orientation=(3,0),pole=3),
             paired_gates(4,1,orientation=(3,0),pole=4)]
    for gates in invalid:
        try:validate_mixed_word(4,gates)
        except ValueError:rejected+=1
        else:raise AssertionError('An overlapping or over-budget word was accepted.')
    return {'status':'Exact mixed-direction paired-word controls passed',
        'initial_X_row_divisions':divisions,'scaled_row_congruences':rows,
        'exact_zero_mixed_intersections':intersections,'original_subword_margins':source_checks,
        'complete_local_profiles':profiles,'independent_rational_profiles':independent,
        'typed_two_copy_coefficients':typed,'rejected_invalid_words':rejected,
        'kernel_cases':kernel_records,'word_cases':records}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('evidence/reflection-mixed-paired-words.json'))
    args=p.parse_args();result=verify()
    names=[Path(__file__).name,'reflection_mixed_paired_words.py','reflection_paired_quartic_words.py',
        'check_reflection_paired_quartic_words.py','reflection_quartic_tiles.py',
        'reflection_dual_rail_gaussian.py','check_reflection_quartic_tiles.py',
        'native_matrix_smith_f4.py','native_matrix_smith_f4.cpp']
    result['source_sha256']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    note=Path(__file__).resolve().parents[1]/'research'/'reflection-mixed-paired-quartics.md'
    if note.exists():result['note_sha256']=sha256(note.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('kernel_cases','word_cases','source_sha256','note_sha256')}))


if __name__=='__main__':main()
