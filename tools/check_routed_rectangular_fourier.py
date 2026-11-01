#!/usr/bin/env python3
"""Exact rectangular Fourier, carry-relay, and local reconstruction controls.

Small controls store explicit arrays. The companion proof gives streaming
table regeneration and its asymptotic workspace bound.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from check_small_word_fourier_constructor import make_crt_parameters,full_dft,integer_convolution,interpolate
from check_rough_length_fourier_constructor import branches,ceildiv,input_word

def messages(p,groups,word):
    rows=p.t//groups;root=pow(p.omega,rows,p.q)
    return [[pow(p.omega,r*k,p.q)*sum(word[j]*pow(root,((j-r)//rows)*k,p.q)
             for j in range(r,len(word),rows))%p.q for k in range(groups)] for r in range(rows)]

def transform(p,groups,x,y):
    rows=p.t//groups;root=pow(p.omega,groups,p.q)
    xm,ym=messages(p,groups,x),messages(p,groups,y)
    xhat,yhat=[0]*p.t,[0]*p.t;z=[[0]*rows for _ in range(groups)]
    for k in range(groups):
        products=[]
        for v in range(rows):
            xhat[k+groups*v]=sum(xm[r][k]*pow(root,r*v,p.q) for r in range(rows))%p.q
            yhat[k+groups*v]=sum(ym[r][k]*pow(root,r*v,p.q) for r in range(rows))%p.q
            products.append(xhat[k+groups*v]*yhat[k+groups*v]%p.q)
        for r in range(rows):
            z[k][r]=pow(rows,-1,p.q)*sum(products[v]*pow(root,(-r*v)%rows,p.q) for v in range(rows))%p.q
    root=pow(p.omega,rows,p.q);coefficient=[0]*p.t
    for r in range(rows):
        for i in range(groups):
            coefficient[r+rows*i]=pow(groups,-1,p.q)*sum(z[k][r]*pow(p.omega,(-r*k)%p.t,p.q)*
                         pow(root,(-i*k)%groups,p.q) for k in range(groups))%p.q
    return coefficient,xhat,yhat,xm,ym,z

def route(position,rows,columns,groups):
    assert 0<=position<rows*columns-1
    block=rows//groups;first=(position%rows)//block;second=((position+1)%rows)//block
    if first==second:return (first,)
    relay=(position//rows)%groups;path=[first]
    if relay not in (first,second):path.append(relay)
    path.append(second)
    return tuple(path)

def routing_controls():
    cases=segments=relays=0
    for exponent in range(1,7):
        groups=1<<exponent
        for row_factor in (1,2,4):
            for column_factor in (1,2,4):
                rows,columns=groups*row_factor,groups*column_factor
                counts=Counter();local_relays=Counter()
                for position in range(rows*columns-1):
                    path=route(position,rows,columns,groups)
                    for a,b in zip(path,path[1:]):
                        assert a!=b
                        counts[tuple(sorted((a,b)))]+=1;segments+=1
                    if len(path)==3:local_relays[path[1]]+=1;relays+=1
                for a in range(groups):
                    for b in range(a+1,groups):
                        assert counts[a,b]==4*columns//groups-(a==0 and b==groups-1)
                    assert local_relays[a]==columns//groups*(groups-2)
                assert sum(local_relays.values())==columns*(groups-2)
                cases+=1
    return dict(parameter_cases=cases,routed_segments=segments,relay_vertices=relays)

def phase_coefficient(bit):return (bit*bit+bit+2)%3

def phase_from_rows(p,rows,local):
    result=0
    for r,word in local.items():
        for i,value in enumerate(word):
            for bit in range(p.d):result+=(value>>bit&1)*phase_coefficient((r+rows*i)*p.d+bit)
    return result%3

def make_boundary(p,groups,x,y):
    rows=p.t//groups;block=rows//groups
    coefficient,xhat,yhat,xm,ym,z=transform(p,groups,x,y)
    assert xhat==full_dft(p,x) and yhat==full_dft(p,y)
    assert coefficient==integer_convolution(x,y,p.t)
    target=sum(value*p.base**i for i,value in enumerate(coefficient));carries=[0]
    for value in coefficient:
        next_carry,digit=divmod(value+carries[-1],p.base)
        assert next_carry<=p.carry_cap
        carries.append(next_carry)
    assert carries[-1]==0
    boundary={}
    for r in range(rows):
        for k in range(groups):
            if r//block!=k:
                boundary['x',r,k]=xm[r][k];boundary['y',r,k]=ym[r][k];boundary['z',k,r]=z[k][r]
    contributions=[phase_from_rows(p,rows,{r:x[r::rows] for r in range(g*block,(g+1)*block)}) for g in range(groups)]
    for position in range(p.t-1):
        path=route(position,rows,groups,groups)
        if len(path)==1:continue
        phase=sum(contributions[:path[0]+1])%3 if position<rows-1 else 0
        packed=carries[position+1]+(phase<<p.carry_bits)
        for leg in range(len(path)-1):boundary['c',position,leg]=packed
    return boundary,target,sum(contributions)%3

def query_group(p,groups,group,lengths,target,boundary):
    """Literal small-instance table query, including all boundary checks."""
    rows=p.t//groups;block=rows//groups;local_range=range(group*block,(group+1)*block)
    root=pow(p.omega,rows,p.q)
    def arithmetic(key):
        value=boundary[key];assert 0<=value<p.q
        return value
    local=[]
    for family,bits in zip(('x','y'),lengths):
        recovered={};length=ceildiv(bits,p.d)
        for r in local_range:
            count=max(0,ceildiv(length-r,rows));outside=[k for k in range(groups) if k!=group]
            assert count<=len(outside)
            chosen=outside[:count]
            coefficients=interpolate([pow(root,k,p.q) for k in chosen],
                 [arithmetic((family,r,k))*pow(pow(p.omega,r*k,p.q),-1,p.q)%p.q for k in chosen],p.q)
            assert all(0<=value<p.base for value in coefficients)
            if count and r+rows*(count-1)==length-1:
                top_bits=bits-(length-1)*p.d
                assert 1<<(top_bits-1)<=coefficients[-1]<1<<top_bits
            for k in outside:
                expected=pow(p.omega,r*k,p.q)*sum(value*pow(root,i*k,p.q) for i,value in enumerate(coefficients))%p.q
                assert arithmetic((family,r,k))==expected
            recovered[r]=coefficients
        local.append(recovered)
    incoming=[]
    for family,recovered in zip(('x','y'),local):
        incoming.append([pow(p.omega,r*group,p.q)*sum(value*pow(root,i*group,p.q)
                         for i,value in enumerate(recovered[r]))%p.q if r in recovered
                         else arithmetic((family,r,group)) for r in range(rows)])
    middle_root=pow(p.omega,groups,p.q)
    product=[(sum(incoming[0][r]*pow(middle_root,r*v,p.q) for r in range(rows))*
              sum(incoming[1][r]*pow(middle_root,r*v,p.q) for r in range(rows)))%p.q for v in range(rows)]
    middle=[pow(rows,-1,p.q)*sum(product[v]*pow(middle_root,(-r*v)%rows,p.q) for v in range(rows))%p.q for r in range(rows)]
    for r in range(rows):
        if r not in local_range:assert arithmetic(('z',group,r))==middle[r]
    coefficients={}
    for r in local_range:
        for i in range(groups):
            coefficients[r+rows*i]=pow(groups,-1,p.q)*sum(
                (middle[r] if k==group else arithmetic(('z',k,r)))*pow(p.omega,(-r*k)%p.t,p.q)*
                pow(root,(-i*k)%groups,p.q) for k in range(groups))%p.q
    def unpack(position,leg):
        value=boundary['c',position,leg];assert 0<=value<1<<(p.carry_bits+2)
        carry=value&((1<<p.carry_bits)-1);phase=value>>p.carry_bits
        assert carry<=p.carry_cap and phase<3
        if position>=rows-1:assert phase==0
        return carry,phase
    contribution=phase_from_rows(p,rows,local[0]);terminal=None
    for i in range(groups):
        start=i*rows+group*block;stop=start+block
        if start==0:carry,phase=0,0
        else:
            previous=route(start-1,rows,groups,groups);carry,phase=unpack(start-1,len(previous)-2)
        for position in range(start,stop):
            carry,digit=divmod(coefficients[position]+carry,p.base)
            assert carry<=p.carry_cap
            assert digit==(target>>(position*p.d))&(p.base-1)
        if i==0:
            phase=(phase+contribution)%3
            if group==groups-1:terminal=phase;phase=0
        if stop==p.t:assert carry==0
        else:assert unpack(stop-1,0)==(carry,phase)
    for position in range(p.t-1):
        path=route(position,rows,groups,groups)
        if len(path)==3 and path[1]==group:assert unpack(position,0)==unpack(position,1)
    return terminal

def transform_controls():
    cases=queries=negative=0;records=[]
    for a,n,groups in ((4,15,4),(8,120,4),(8,120,8),(16,500,4),(16,500,8),(16,500,16)):
        p,primes=make_crt_parameters(a,ceildiv(n+1,a*a));rows=p.t//groups;s=p.t*p.d//groups
        assert rows%groups==0 and p.carry_bits+2<=p.f
        lengths=sorted({(s+1,n-s),(n-s,s+1),(n//2,n-n//2)});local_cases=0
        for k,l in lengths:
            assert (k,l) in branches(n,s)
            for xvalue,yvalue in (((1<<k)-1,(1<<l)-1),((1<<(k-1))+1,(1<<(l-1))+1)):
                x,y=input_word(xvalue,k,p),input_word(yvalue,l,p)
                boundary,target,phase=make_boundary(p,groups,x,y);assert target==xvalue*yvalue
                for group in range(groups):
                    terminal=query_group(p,groups,group,(k,l),target,boundary)
                    assert terminal==(phase if group==groups-1 else None)
                    queries+=1
                rejected=0
                for group in range(groups):
                    try:query_group(p,groups,group,(k,l),target+1,boundary)
                    except AssertionError:rejected+=1
                assert rejected;negative+=1
                relay=next((j for j in range(p.t-1) if len(route(j,rows,groups,groups))==3),None)
                if relay is not None:
                    altered=dict(boundary);altered['c',relay,0]^=1
                    try:query_group(p,groups,route(relay,rows,groups,groups)[1],(k,l),target,altered)
                    except AssertionError:negative+=1
                    else:raise AssertionError('Corrupted relay was accepted.')
                local_cases+=1;cases+=1
        width=6*p.t//(groups*groups)*p.f+4*(p.carry_bits+2)
        records.append(dict(A=a,R=rows,S=groups,P=groups,d=p.d,f=p.f,carry_bits_with_phase=p.carry_bits+2,
                            width=width,crt_primes=primes,arithmetic_modulus=str(p.q),cases=local_cases))
    return dict(transform_pairs=cases,complete_group_queries=queries,rejected_false_targets_or_relays=negative,parameters=records)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=dict(status='passed',factoring_exponent_improved=False,routing=routing_controls(),transforms=transform_controls(),
                source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                 for name in ('check_routed_rectangular_fourier.py','check_small_word_fourier_constructor.py','check_rough_length_fourier_constructor.py')},
                scope='Exact constructor identities, balanced carry routes, complete local boundary reconstruction, and scalar phase propagation; asymptotic bounds are proved in the companion note.')
    if args.output:args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
