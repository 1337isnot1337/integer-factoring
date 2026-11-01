#!/usr/bin/env python3
"""Complete F3 minor controls for rounded weak pole potentials.

The test is a scalar rational crossing with extra pole at w=2. Its weak
paired potential is necessarily (1/2,1/2), not in Z[1/3]. Exact sparse
determinants independently compute every weighted Smith prefix.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path


def add(a,b,sign=1):
    c=dict(a)
    for i,v in b.items():
        c[i]=(c.get(i,0)+sign*v) % 3
        if not c[i]:
            del c[i]
    return c


def multiply(a,b):
    c={}
    for i,x in a.items():
        for j,y in b.items():
            c[i+j]=(c.get(i+j,0)+x*y) % 3
    return {i:v for i,v in c.items() if v}


def power(a,k):
    if k<0:
        raise ValueError('A negative exponent is not a polynomial.')
    c={0:1}
    for _ in range(k):
        c=multiply(c,a)
    return c


def determinant(matrix):
    n=len(matrix)
    result={}
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term={0:1}
        for i,j in enumerate(p):
            term=multiply(term,matrix[i][j])
        result=add(result,term,sign)
    return result


def profile(tensor,infinity=False,heights=None):
    if heights is None:
        heights=[[0,0] for _ in range(4)]
    pairs=list(product(range(2),repeat=2))
    if infinity:
        matrix=[[tensor[x,y,u,v] for y,v in pairs] for x,u in pairs]
        row=[heights[0][x]+heights[2][u] for x,u in pairs]
        col=[heights[1][y]+heights[3][v] for y,v in pairs]
    else:
        matrix=[[tensor[x,y,u,v] for u,v in pairs] for x,y in pairs]
        row=[heights[0][x]+heights[1][y] for x,y in pairs]
        col=[heights[2][u]+heights[3][v] for u,v in pairs]
    prefixes=[0]
    checked=0
    for k in range(1,5):
        candidates=[]
        for rows,cols in product(combinations(range(4),k),repeat=2):
            value=determinant([[matrix[i][j] for j in cols] for i in rows])
            checked+=1
            if value:
                candidates.append((-max(value) if infinity else min(value))+
                                  sum(row[i] for i in rows)+sum(col[j] for j in cols))
        assert candidates
        prefixes.append(min(candidates))
    return [b-a for a,b in zip(prefixes,prefixes[1:])],checked


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    A,B={0:1,1:1},{0:2,1:1} # A=w-2, B=w-1 over F3.
    valuations={(1,1,0,0):-1,(0,0,1,1):1,(1,0,1,0):0,(0,1,0,1):0}
    source={}
    for entry in product(range(2),repeat=4):
        x,y,u,v=entry
        if x==v and y==u:
            source[entry]=A
        elif entry in valuations:
            source[entry]=multiply(B,power(A,valuations[entry]+1))
        else:
            source[entry]={}
    # source=(w-2)*K. This clears one original scalar denominator.
    original,minor_count=[],0
    for infinity in (False,True):
        orders,checked=profile(source,infinity)
        original.append([v+int(infinity) for v in orders])
        minor_count+=checked
    for entry,nu in valuations.items():
        x,y,u,v=entry
        assert nu+Fraction(x-v+y-u,2)==0
    # The four weak inequalities force f+g=1 and f-g=0.
    assert Fraction(1,2).denominator % 3 != 0
    records=[]
    endpoint_entries=0
    for D in (1,3,9,27):
        H=D//2
        tensor={}
        residual=[]
        for entry in product(range(2),repeat=4):
            x,y,u,v=entry
            if x==v and y==u:
                value=power(A,4)
            elif entry in valuations:
                charge=x-v+y-u
                exponent=D*valuations[entry]+H*charge
                residual.append(exponent)
                assert exponent>=-4
                value=multiply({0:2,D:1},power(A,exponent+4))
            else:
                value={}
            tensor[entry]=value
            assert sum(value.values()) % 3 == int(x==v and y==u)
            endpoint_entries+=1
        assert min(residual)==-1 # Rounding without the charged correction fails.
        heights=[[0,H],[0,H],[0,-H],[0,-H]]
        zero,checked=profile(tensor,False)
        minor_count+=checked
        strong,checked=profile(tensor,True,heights)
        minor_count+=checked
        assert zero==[D*v for v in original[0]]
        assert strong==[D*v-4 for v in original[1]]
        assert all(heights[0][i]+heights[3][i]==0 for i in range(2))
        assert all(heights[1][i]+heights[2][i]==0 for i in range(2))
        records.append({'Frobenius_degree':D,'rounded_potential':H,
                        'uncleared_extra_pole_valuations':residual,
                        'degree':max(max(v,default=0) for v in tensor.values()),
                        'compensating_infinity_heights':heights,
                        'zero_orders':zero,'infinity_orders':strong,
                        'margin_penalty_per_unit_scale':str(Fraction(4,D)),
                        'full_polynomial_tensor':{','.join(map(str,k)):v for k,v in tensor.items()}})
    result={'status':'pass','field_characteristic':3,'physical_alphabet':2,
            'original_zero_infinity_orders':original,'complete_minor_determinants':minor_count,
            'independent_complete_orders':40,'strict_endpoint_entries':endpoint_entries,
            'required_weak_potential':['1/2','1/2'],'extra_pole_total_degree':1,
            'records':records,'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact finite control of rounded weak-potential strictification, '
                    'including its charged degree error. Not a favorable crossing or new score.'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','complete_minor_determinants',
                     'independent_complete_orders','strict_endpoint_entries')}))


if __name__=='__main__':
    main()
