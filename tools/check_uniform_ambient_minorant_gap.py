#!/usr/bin/env python3
"""Exact band and composition controls for the uniform ambient obstruction."""
import argparse
from fractions import Fraction as F
import json
from math import gcd,isqrt
from pathlib import Path
from binary_quadratic_forms import compose_forms,reduce_form
from prime_form_bias import primes_below
from check_prime_translate_obstruction import smooth_remainder


def band_control():
    # Coefficients of t,t^2 on the three intervals, from the original w terms.
    weights=[F(1),-F(21,16),F(27,64)]
    slopes=[F(1),F(4,3),F(16,9)]
    pieces=[]
    for active in (3,2,1):
        pieces.append([sum(weights[j]*slopes[j] for j in range(active)),
                       -sum(weights[j]*slopes[j]**2 for j in range(active))])
    assert pieces==[[0,0],[-F(3,4),F(4,3)],[1,-1]]
    integral=sum(weights[j]/slopes[j]/6 for j in range(3))
    assert integral==F(259,6144)>0
    assert sum(abs(w) for w in weights)==F(175,64)<3
    return dict(piece_coefficients=[[str(v) for v in p] for p in pieces],integral=str(integral),absolute_combination_mass=str(sum(abs(w) for w in weights)))


def prime_forms(D,lower,upper):
    for p in primes_below(upper):
        if p<=lower or gcd(p,D)!=1:continue
        for B in range(-p+1,p+1):
            if (B*B-D)%(4*p)==0:
                yield p,reduce_form((p,B,(B*B-D)//(4*p)))


def geometry_control(D):
    assert D<0 and D%4 in (0,1)
    Z=16;Y=13;cap=15;size=-D;X=isqrt(size)//(8*Z);H=isqrt(size//3)
    assert 9*X>32*Z and X>2
    A,target=next(prime_forms(D,F(9,8)*Z,2*Z))
    assert A>Z>=max(Y,cap) and A<2*Z
    records=[]
    for p,form in prime_forms(D,F(9,16)*X,X):
        assert p>A and gcd(A*p,size)==1 and 16*(A*p)**2<size
        translated=compose_forms(target,form);a,B,C=translated
        assert a==A*p and a<=H and gcd(a,size)==1
        assert smooth_remainder(a,Y)>cap
        assert abs(B)!=a and B!=0 and a!=C  # No ambiguous-class shortcut.
        records.append(dict(prime=p,source_form=form,translated_form=translated,residual=smooth_remainder(a,Y)))
        if len(records)==4:break
    assert len(records)==4
    return dict(discriminant=D,Z=Z,smooth_bound=Y,cofactor_cap=cap,auxiliary_prime=A,target_form=target,band_cutoff=X,full_height=H,oriented_translates=records)


def main():
    # Different fundamental/nonmaximal orders; no specially chosen middle=1 form.
    controls=[geometry_control(D) for D in (-1000000007,-4000000028,
        -1000036000099,-4000144000396,-9000324000891)]
    result=dict(scope='Exact band identity and fixed arbitrary-discriminant geometry controls. These finite cutoffs do not certify the uniform analytic estimates or their asymptotic thresholds.',band=band_control(),controls=controls,oriented_translates_checked=sum(len(r['oriented_translates']) for r in controls),binary_decay_denominator=64)
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(discriminants=len(controls),oriented_translates_checked=result['oriented_translates_checked'],band_integral=result['band']['integral'],binary_decay_denominator=64)))


if __name__=='__main__':main()
