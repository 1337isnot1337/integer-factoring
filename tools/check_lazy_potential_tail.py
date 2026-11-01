#!/usr/bin/env python3
"""Exact rational controls for the uniform lazy-potential tail proof."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from math import comb,factorial
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/lazy-potential-tail-controls.json'


def report():
    e07=sum((Q(7,10)**j/factorial(j) for j in range(7)),Q(0))
    e3=sum((Q(3)**j/factorial(j) for j in range(9)),Q(0))
    assert e07>2 and e3>19
    assert Q(4,5)*Q(23,6)*2==Q(92,15)
    assert Q(4,5)*Q(23,6)**2+Q(142,15)==Q(191,9)<22
    assert Q(3,2)*Q(7,10)**2>Q(7,10)
    assert 1-Q(9,64)==Q(55,64)>Q(1,2)
    rademacher=[]
    for n in range(1,257):
        distribution={2*j-n:comb(n,j) for j in range(n+1)}
        tail=0;worst=None;which=None
        for k in range(n-1,-1,-1):
            tail+=distribution.get(k+1,0)
            A=Q(3*(k+1)**2,2*n)+70
            exponent=A.numerator//A.denominator
            # exp(-A) <= 2^(-floor(A)), since e>2.
            ratio=Q(tail*(1 << exponent),1 << n)
            assert ratio>=1
            if worst is None or ratio<worst:worst,which=ratio,k
        rademacher.append(dict(n=n,all_k_count=n,smallest_dyadic_ratio=str(worst),at_k=which))
    walks=[]
    for p in (Q(1,4),Q(1,3),Q(3,10),Q(4,15),Q(1,2)):
        num,den=p.numerator,p.denominator
        for w in (64,100,256):
            distribution={0:1}
            for _ in range(w):
                changed={}
                for z,c in distribution.items():
                    for shift,weight in ((1,num),(-1,num),(0,den-2*num)):
                        if weight:changed[z+shift]=changed.get(z+shift,0)+c*weight
                distribution=changed
            total=den**w;kmax=w//4-1
            tail=sum(c for z,c in distribution.items() if z>kmax+1)
            worst=None;which=None;witness=None
            for k in range(kmax,-1,-1):
                tail+=distribution.get(k+1,0)
                A=Q(9*(k+1)**2,2*w)+70
                exponent=A.numerator//A.denominator+1
                ratio=Q(tail*(1 << exponent),total)
                assert ratio>=1
                if worst is None or ratio<worst:
                    worst,which,witness=ratio,k,dict(tail_numerator=tail,denominator=total,
                                                    dyadic_exponent=exponent)
            walks.append(dict(p=str(p),w=w,k_count=kmax+1,smallest_dyadic_ratio=str(worst),
                              at_k=which,witness=witness))
    source='tools/check_lazy_potential_tail.py'
    return dict(scope='Exact finite controls for an ordinary all-parameter lazy-walk proof; no arbitrary-algorithm lower bound or stronger factoring score.',
                exponential_series_certificates=[dict(x='7/10',terms=7,sum=str(e07),strict_lower=2),
                                                  dict(x='3',terms=9,sum=str(e3),strict_lower=19)],
                completed_square=dict(coefficient='4/5',center='23/6',constant='191/9'),
                rademacher_controls=rademacher,lazy_walk_controls=walks,
                source_sha256={source:sha256((ROOT/source).read_bytes()).hexdigest()})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],rademacher_k_cases=sum(x['all_k_count'] for x in result['rademacher_controls']),
                          lazy_walk_k_cases=sum(x['k_count'] for x in result['lazy_walk_controls']),
                          exponential_rational_certificates=2)))


if __name__=='__main__':main()
