#!/usr/bin/env python3
"""Exact controls for the all-field Boolean-incidence support barrier."""
import argparse
from functools import lru_cache
import hashlib,json
from pathlib import Path
from random import Random


def binary_rank(rows):
    pivots={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p in pivots:row^=pivots[p]
            else:pivots[p]=row;break
    return len(pivots)


@lru_cache(maxsize=None)
def convolution_rank(word,width):
    d=1<<width
    return binary_rank([sum(((word>>(a^b))&1)<<b for b in range(d)) for a in range(d)])


def cuts(kernel,width):
    d=1<<width
    columns=[sum(((kernel[s]>>w)&1)<<s for s in range(d)) for w in range(d)]
    return (sum(convolution_rank(row,width) for row in kernel),
            sum(convolution_rank(column,width) for column in columns))


def support_bounds(kernel,width,incidence=False):
    d=1<<width;r,c=cuts(kernel,width)
    if not any(kernel):return r,c
    u=sum(bool(row) for row in kernel)
    union=0
    for row in kernel:union|=row
    v=union.bit_count();m=sum(row.bit_count() for row in kernel)
    assert r*m>=d*u*u and c*m>=d*v*v and r*c>=d*d
    if incidence:
        assert all(not ((row>>w)&1) or not (w&~s) for s,row in enumerate(kernel) for w in range(d))
        diagonal=[s for s,row in enumerate(kernel) if row>>s&1]
        rr=sum(1<<(width-s.bit_count()) for s in diagonal)
        cc=sum(1<<s.bit_count() for s in diagonal)
        assert r>=rr and c>=cc and r*c>=d*len(diagonal)**2
    return r,c


def prime_rank(matrix,p):
    rows=[r[:] for r in matrix];height=len(rows);width=len(rows[0]);answer=0
    for col in range(width):
        chosen=next((i for i in range(answer,height) if rows[i][col]%p),None)
        if chosen is None:continue
        rows[answer],rows[chosen]=rows[chosen],rows[answer]
        inverse=pow(rows[answer][col],-1,p)
        rows[answer]=[(v*inverse)%p for v in rows[answer]]
        for i in range(answer+1,height):
            scalar=rows[i][col]
            if scalar:rows[i]=[(a-scalar*b)%p for a,b in zip(rows[i],rows[answer])]
        answer+=1
        if answer==height:break
    return answer


def direct_tensor_cuts(kernel,width):
    d=1<<width;rows=[0]*(d*d);cols=[0]*(d*d)
    for x in range(d):
        for y in range(d):
            for u in range(d):
                v=x^y^u
                entry=(kernel[x^y]>>(x^u))&1
                rows[x*d+y]|=entry<<(u*d+v)
                cols[x*d+u]|=entry<<(y*d+v)
    return binary_rank(rows),binary_rank(cols)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    args=parser.parse_args();rng=Random(20261030125)
    arbitrary=incidence=random_binary=prime_controls=dense=0
    for value in range(1<<16):
        support_bounds([(value>>(4*s))&15 for s in range(4)],2)
        arbitrary+=1
    for width in (1,2):
        d=1<<width;positions=[(s,w) for s in range(d) for w in range(d) if not w&~s]
        assert len(positions)==3**width
        for mask in range(1<<len(positions)):
            kernel=[0]*d
            for i,(s,w) in enumerate(positions):kernel[s]|=((mask>>i)&1)<<w
            support_bounds(kernel,width,True);incidence+=1
    for width in range(1,6):
        d=1<<width
        for sample in range(40):
            kernel=[sum(rng.randrange(2)<<w for w in range(d) if not w&~s) for s in range(d)]
            support_bounds(kernel,width,True);random_binary+=1
            if width<=3 and sample<6:
                assert direct_tensor_cuts(kernel,width)==cuts(kernel,width);dense+=1
    for width in range(1,5):
        d=1<<width
        for _ in range(8):
            kernel=[[rng.randrange(3) if not w&~s else 0 for w in range(d)] for s in range(d)]
            r=sum(prime_rank([[kernel[s][a^b] for b in range(d)] for a in range(d)],3) for s in range(d))
            c=sum(prime_rank([[kernel[a^b][w] for b in range(d)] for a in range(d)],3) for w in range(d))
            diagonal=[s for s in range(d) if kernel[s][s]]
            assert r>=sum(1<<(width-s.bit_count()) for s in diagonal)
            assert c>=sum(1<<s.bit_count() for s in diagonal)
            assert r*c>=d*len(diagonal)**2
            prime_controls+=1
    for width in range(1,6):
        d=1<<width;kernel=[(1<<d)-1]*d
        assert support_bounds(kernel,width)==(d,d)
        assert d*d<d**3
    result={'status':'PASS','scope':'All-field support proof; no improved factoring exponent',
            'all_binary_kernels_at_four_labels':arbitrary,'exhaustive_binary_incidence_kernels':incidence,
            'random_binary_incidence_controls':random_binary,'prime_three_incidence_controls':prime_controls,
            'independent_dense_tensor_controls':dense,'off_incidence_negative_controls':5,
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
