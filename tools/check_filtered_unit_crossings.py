#!/usr/bin/env python3
"""Exact controls for the complete filtered-unit crossing constructor."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
from functools import reduce
import json
from pathlib import Path
import random

from filtered_unit_crossings import (coefficient_polynomials,coefficient_ranks,
    dense_coefficient,group_kernel,multiply_field,parameter_positions,rank,tensor_kernels)

ROOT=Path(__file__).resolve().parents[1]


def multiply_exterior(first,second,width):
    out={}
    for (i,I),a in first.items():
        for (j,J),b in second.items():
            if not I&J:
                key=(i+j,I|J);out[key]=out.get(key,0)^multiply_field(a,b)
    return {key:value for key,value in out.items() if value}


def selector_curve(width,mask):
    curve={(0,0):1};state=0
    for p in range(1,1<<width):
        if mask>>(p-1)&1:
            state^=p;factor={(0,0):1};J=p
            while J:
                factor[1,J]=1;J=(J-1)&p
            curve=multiply_exterior(curve,factor,width)
    assert all(j<=J.bit_count() for j,J in curve)
    parameters=[curve.get((j,J),0) for J,j in parameter_positions(width)]
    reconstructed=coefficient_polynomials(width,state,parameters)
    assert curve=={(j,J):value for j,row in enumerate(reconstructed) for J,value in enumerate(row) if value}
    return state,parameters


def exact_endpoints(kernels):
    width=len(kernels[0])-1;dimension=len(kernels)
    for state,coeffs in enumerate(kernels):
        assert coeffs[-1]==[int(w==0) for w in range(dimension)]
        assert [reduce(int.__xor__,[coeffs[k][w] for k in range(width+1)],0)
                for w in range(dimension)]==[int(w==state) for w in range(dimension)]
        augmentation=[]
        for row in coeffs:
            value=0
            for scalar in row:value^=scalar
            augmentation.append(value)
        assert augmentation==[0]*width+[1]


def literal_square(coefficients):
    width=len(coefficients)-1;dimension=len(coefficients[0])
    answer=[[0]*dimension for _ in range(2*width+1)]
    for j,row in enumerate(coefficients):
        for k,other in enumerate(coefficients):
            for x,a in enumerate(row):
                if a:
                    for y,b in enumerate(other):
                        if b:answer[j+k][x^y]^=multiply_field(a,b)
    expected=[[0]*dimension for _ in range(2*width+1)];expected[-1][0]=1
    assert answer==expected


def main():
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=Path('/tmp/filtered-unit-crossings.json'));args=pa.parse_args()
    rng=random.Random(197721);constructor_states=square_controls=dense_ranks=0;profiles=[]
    for width in range(1,7):
        dimension=1<<width;positions=parameter_positions(width)
        assert len(positions)==((width-2)*(1<<(width-1))+1)
        for trial in range(4):
            params=[[rng.randrange(4) for _ in positions] for _ in range(dimension)]
            kernels=tensor_kernels(width,params);exact_endpoints(kernels);constructor_states+=dimension
            for state in range(min(dimension,4)):
                literal_square(kernels[state]);square_controls+=1
            if width<=3:
                rows,columns=coefficient_ranks(kernels)
                for power in range(width+1):
                    for cut,expected in enumerate([rows[power],columns[power]]):
                        assert rank(dense_coefficient(kernels,power,cut))==expected;dense_ranks+=1
                profiles.append({'width':width,'row_ranks':rows,'column_ranks':columns})
    selector_controls=0;zero_selectors=set()
    for width in range(1,4):
        for mask in range(1<<((1<<width)-1)):
            state,params=selector_curve(width,mask);selector_controls+=1
            if width==3 and state==0:zero_selectors.add(tuple(tuple(x) for x in group_kernel(width,state,params)))
    witness=group_kernel(3,0,[0,0,0,0,1])
    assert tuple(tuple(x) for x in witness) not in zero_selectors and len(zero_selectors)==16
    family_profiles=Counter()
    for free in product(range(4),repeat=4):
        kernels=tensor_kernels(2,[[x] for x in free]);exact_endpoints(kernels)
        rows,columns=coefficient_ranks(kernels)
        for power in range(3):
            for cut,expected in enumerate([rows[power],columns[power]]):
                assert rank(dense_coefficient(kernels,power,cut))==expected;dense_ranks+=1
        family_profiles[tuple(rows+columns)]+=1
    candidate=[[0]*5 for _ in range(8)]
    candidate[3][0]=candidate[5][1]=candidate[6][2]=1
    candidate[7]=[1]*5
    kernels=tensor_kernels(3,candidate);exact_endpoints(kernels)
    rows,columns=coefficient_ranks(kernels)
    assert rows==[1,0,28,64] and columns==[64,0,28,1]
    for power in range(4):
        for cut,expected in enumerate([rows[power],columns[power]]):
            assert rank(dense_coefficient(kernels,power,cut))==expected;dense_ranks+=1
    report={'scope':'Exact construction and finite controls only; no improved factoring score or certified weighted pressure.',
        'random_constructor_states':constructor_states,'literal_group_polynomial_squares':square_controls,
        'independent_dense_cut_ranks':dense_ranks,'all_selector_products_through_width_three':selector_controls,
        'binary_zero_state_selector_kernels':len(zero_selectors),
        'strict_enlargement_witness':{'width':3,'state':0,'free_parameters':[0,0,0,0,1],'group_coefficients':witness},
        'complete_width_two_GF4_families':256,
        'width_two_rank_profile_counts':[{'row_ranks':list(key[:3]),'column_ranks':list(key[3:]),'count':value} for key,value in sorted(family_profiles.items())],
        'finite_candidate':{'affine_coefficients':candidate,'row_ranks':rows,'column_ranks':columns,'histogram_cost':44},
        'random_profiles':profiles,
        'source_sha256':{path:sha256((ROOT/path).read_bytes()).hexdigest() for path in ['tools/filtered_unit_crossings.py','tools/check_filtered_unit_crossings.py']}}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if isinstance(v,(int,str))}),flush=True)
if __name__=='__main__':main()
