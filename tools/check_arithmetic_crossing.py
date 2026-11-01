#!/usr/bin/env python3
"""Independent finite controls for arithmetic crossing initial forms."""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from math import comb,factorial
from pathlib import Path
from fused_rmatrix_probe import (fused_polynomial,arithmetic_initial_form,
    initial_form,cut_blocks,smith,translated,ptrim,profile_initial)


def rank_mod(matrix,p):
    a=[[x%p for x in row] for row in matrix];m=len(a);n=len(a[0]) if a else 0;r=0
    for j in range(n):
        pi=next((i for i in range(r,m) if a[i][j]),None)
        if pi is None:continue
        a[r],a[pi]=a[pi],a[r];inv=pow(a[r][j],-1,p)
        a[r]=[x*inv%p for x in a[r]]
        for i in range(r+1,m):
            co=a[i][j]
            if co:a[i]=[(x-co*y)%p for x,y in zip(a[i],a[r])]
        r+=1
        if r==m:break
    return r


def lower(vector,q):
    out=[0]*(q*q);k=q-1
    for i in range(q):
        for j in range(q):
            co=vector[i*q+j]
            if i<k:out[(i+1)*q+j]+=co*(k-i)
            if j<k:out[i*q+j+1]+=co*(k-j)
    return out


def source_control(k,raw):
    q=k+1;n=q*q;basis=[];checks=0
    for spin in range(k+1):
        w=k-spin;vector=[0]*n
        for i in range(w+1):vector[i*q+w-i]=(-1)**i*comb(w,i)
        polynomial=[1]
        for l in range(1,k+1):
            out=[0]*(len(polynomial)+1)
            for i,c in enumerate(polynomial):out[i]+=c*(l if l<=spin else -l);out[i+1]+=c
            polynomial=out
        for step in range(2*spin+1):
            assert any(vector);basis.append(vector)
            for degree in range(k+1):
                result=[Fraction(0)]*n
                for (x,y,u,v),entry in raw.items():result[x*q+y]+=entry[degree]*vector[u*q+v]
                assert result==[polynomial[degree]*v for v in vector]
                checks+=n
            assert [vector[j*q+i] for i in range(q) for j in range(q)]==[(-1)**(k-spin)*v for v in vector]
            vector=lower(vector,q)
        assert not any(vector)
    assert len(basis)==n
    assert rank_mod(basis,1000003)==n
    return checks


def toeplitz_rank(block,precision,p):
    n=len(block);matrix=[[0]*(n*precision) for _ in range(n*precision)]
    for output_degree in range(precision):
        for input_degree in range(output_degree+1):
            degree=output_degree-input_degree
            for i in range(n):
                for j in range(n):
                    entry=block[i][j]
                    if degree<len(entry):matrix[output_degree*n+i][input_degree*n+j]=entry[degree]
    return rank_mod(matrix,p)


def audit_initial(q,p,heights=None):
    k=q-1;raw,denominator=fused_polynomial(k)
    source_checks=source_control(k,raw)
    values,degree,slope=initial_form(raw,k,p,heights)
    # Check the endpoint after z -> 1-z using ordinary binomial sums.
    endpoints=0
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    source=values.get((x,y,u,v),[]);result=[0]*(degree+1)
                    for d,c in enumerate(source):
                        for j in range(d+1):result[j]=(result[j]+c*(-1)**j*comb(d,j))%p
                    assert sum(result)%p==int(x==v and y==u)
                    endpoints+=1
    toeplitz=0;all_profiles=[]
    for cut in (0,1):
        for point in [1,'infinity']:
            orders=[]
            for block in cut_blocks(values,q,cut):
                if point=='infinity':local=[[(entry+[0]*(degree+1-len(entry)))[::-1] for entry in row] for row in block]
                else:local=[[translated(entry,point,p) for entry in row] for row in block]
                spectrum=smith(local,p,len(local)*degree+1)
                for precision in range(1,max(spectrum)+2):
                    assert toeplitz_rank(local,precision,p)==sum(max(0,precision-v) for v in spectrum)
                    toeplitz+=1
                orders+=spectrum
            all_profiles.append({'cut':cut,'point':point,'orders':sorted(orders)})
    result=profile_initial(values,degree,q,p)
    assert result['parameter_pairs_checked']==result['parameter_pairs_excluded_exactly']
    return {'q':q,'prime':p,'degree':degree,'critical_slope':str(slope),'heights':heights,
        'source_denominator':denominator,'source_spectral_entry_checks':source_checks,
        'endpoint_checks':endpoints,'independent_toeplitz_checks':toeplitz,
        'profiles':all_profiles,'best_local_profile':result}


def controls():
    # A ramified slope 1/2 is legal even though no such scalar lies in Q_3.
    q=2;p=3;raw={}
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    raw[x,y,u,v]=(Fraction(x==v and y==u),Fraction(0),Fraction(int(x==u and y==v),3))
    value,degree,slope=arithmetic_initial_form(raw,q,p)
    assert degree==2 and slope==Fraction(1,2)
    for key,entry in raw.items():assert value.get(key,[])==ptrim([int(entry[0]),0,int(entry[2]*3)])
    # Malformed constant terms and a composite coefficient modulus fail.
    bad=dict(raw);bad[0,0,0,0]=(Fraction(2),)
    failures=0
    for args in ((bad,q,p),(raw,q,4)):
        try:arithmetic_initial_form(*args)
        except AssertionError:failures+=1
        else:raise AssertionError('negative control accepted')
    return {'fractional_slope':'1/2','negative_controls':failures}


def main():
    a=argparse.ArgumentParser();a.add_argument('--output',default='arithmetic-crossing-audit.json');args=a.parse_args()
    cases=[audit_initial(2,2),audit_initial(3,3),audit_initial(5,2),
        audit_initial(6,2,([0,0,0,0,1,1],[0,0,0,0,1,1]))]
    report={'status':'passed','controls':controls(),'cases':cases,
        'scope':'exact construction, independent source eigenspaces, literal endpoints, and ordinary Toeplitz ranks; no improved factoring exponent'}
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
