#!/usr/bin/env python3
"""Exact quotient-mixing and short-product injection controls; no orbit search."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from math import gcd,isqrt,prod
from pathlib import Path
from binary_quadratic_forms import compose_forms,reduce_form
from prime_form_bias import primes_below


def quotient_control():
    n=16;choices=list(range(1,n));L=3
    character=lambda x,c:(-1)**((x&c).bit_count()%2)
    eta=max(abs(sum(F(character(x,c),n-1) for x in choices)) for c in range(1,n))
    assert eta==F(1,15)
    counts=Counter(a^b^c for a,b,c in product(choices,repeat=L))
    beta=4*eta;rows=[]
    for mask in (1,3,7,15):
        index=1<<mask.bit_count()
        hits=[sum(v for g,v in counts.items() if g&mask==coset) for coset in range(index)]
        assert index*beta**L<=F(1,2)
        for count in hits:
            probability=F(count,len(choices)**L)
            assert probability>=F(1,index)-eta**L>=F(1,2*index)
            # Each color has 15 equal-weight raw choices in its disjoint copy.
            assert count>=F(len(choices)**L,2*index)
        rows.append(dict(index=index,coset_counts=hits,minimum_probability=str(F(min(hits),15**L)),guaranteed_probability=str(F(1,2*index))))
    return dict(group='(Z/2)^4',colors=L,raw_choices_per_color=15,per_color_bias=str(eta),coarse_color_bias=str(beta),raw_tuples=15**L,records=rows)


def injection_control():
    D=-1000036000099;bound=isqrt(-D)//4;eligible=[]
    for p in primes_below(101):
        if gcd(p,D)!=1:continue
        roots=[b for b in range(-p+1,p+1) if (b*b-D)%(4*p)==0]
        if roots:eligible.append((p,roots))
        if len(eligible)==6:break
    assert len(eligible)==6
    colors=[[(p,b) for p,roots in eligible[2*j:2*j+2] for b in roots] for j in range(3)]
    forms=set();records=[]
    for choices in product(*colors):
        norm=prod(p for p,_ in choices)
        assert norm<=bound and gcd(norm,D)==1
        B=sum(b*(norm//p)*pow(norm//p,-1,p) for p,b in choices)%norm
        if (B-D)%2:B+=norm
        if B>norm:B-=2*norm
        assert -norm<B<=norm and (B*B-D)%(4*norm)==0
        literal=(norm,B,(B*B-D)//(4*norm))
        assert literal[2]>literal[0] and reduce_form(literal)==literal
        current=(1,1,(1-D)//4)
        for p,b in choices:current=compose_forms(current,reduce_form((p,b,(b*b-D)//(4*p))))
        assert current==literal and current not in forms
        forms.add(current);records.append(dict(oriented_primes=choices,reduced_form=current))
    assert len(forms)==64
    return dict(discriminant=D,norm_bound=bound,underlying_primes=[p for p,_ in eligible],colors=3,oriented_tuples=64,distinct_reduced_classes=len(forms),records=records)


def main():
    assert F(1,2)*F(1,32)==F(1,64)
    q=quotient_control();a=injection_control()
    result=dict(scope='Exact controls for quotient mixing and injective distinct-prime products below the reduction threshold. They do not infer the analytic estimate, probabilistic partition lemma or uniform asymptotic density from finite data.',quotient=q,injection=a,threshold_exponent='1/64')
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(quotient_tuples=q['raw_tuples'],quotient_subgroups=len(q['records']),distinct_arithmetic_classes=a['distinct_reduced_classes'],threshold_exponent='1/64')))


if __name__=='__main__':main()
