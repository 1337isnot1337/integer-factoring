#!/usr/bin/env python3
"""Exact finite-group controls for symmetric autocorrelation on a target cycle.

The implementation uses groups of exponent 2,3,4 or6, whose symmetric
character sums have rational real parts. The accompanying proof is general.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import product
from math import lcm

COS={2:(F(1),F(-1)),3:(F(1),F(-1,2),F(-1,2)),
     4:(F(1),F(0),F(-1),F(0)),6:(F(1),F(1,2),F(-1,2),F(-1),F(-1,2),F(1,2))}


def elements(shape):return list(product(*(range(d) for d in shape)))
def subtract(a,b,shape):return tuple((x-y)%d for x,y,d in zip(a,b,shape))
def phase(k,x,shape):return sum((F(a*b,d) for a,b,d in zip(k,x,shape)),F(0))%1
def cosine(theta,order):
    index=(theta%1)*order
    assert index.denominator==1
    return COS[order][int(index)]


def analyze(shape,A,g):
    A=set(A);zero=(0,)*len(shape);group=elements(shape);order=lcm(*shape)
    assert order in COS and A and zero in A
    assert {subtract(zero,a,shape) for a in A}==A
    cycle=[];x=zero
    while not cycle or x!=zero:
        cycle.append(x);x=tuple((a+b)%d for a,b,d in zip(x,g,shape))
        assert len(cycle)<=len(group)
    assert len(cycle)>1
    counts=Counter(subtract(a,b,shape) for a in A for b in A)
    f={x:F(counts[x],len(A)) for x in group}
    spectrum=defaultdict(F)
    for chi in group:
        sigma=sum((cosine(phase(chi,a,shape),order) for a in A),F(0))/len(A)
        spectrum[phase(chi,g,shape)]+=sigma*sigma
    c0=spectrum[F(0)];assert c0>0
    pairs=sum(subtract(a,b,shape) in cycle for a in A for b in A)
    assert c0==F(len(group)*pairs,len(cycle)*len(A)**2)
    for j,x in enumerate(cycle):
        actual=F(len(A),len(group))*sum((c*cosine(j*theta,order) for theta,c in spectrum.items()),F(0))
        assert actual==f[x]
    windows=[]
    for W in range(1,2*len(cycle)+1):
        nonzero={min(theta,1-theta):F(0) for theta in spectrum if theta}
        for theta,c in spectrum.items():
            if theta:nonzero[min(theta,1-theta)]+=c
        error=sum((c*min(F(1),1/(2*W*eta)) for eta,c in nonzero.items()),F(0))
        lower=F(len(A)*W,len(group))*(c0-error)
        actual=[sum((f[cycle[(start+j)%len(cycle)]] for j in range(W)),F(0))
                for start in range(len(cycle))]
        assert min(actual)>=lower
        if error<c0:assert min(actual)>0
        # Exact layer-cake integration of the stepwise cumulative spectral mass.
        u0=F(1,2*W);points=sorted({u0,F(1,2)}|{eta for eta in nonzero if u0<eta<F(1,2)})
        integral=F(0)
        for left,right in zip(points,points[1:]):
            midpoint=(left+right)/2
            mass=sum((c for eta,c in nonzero.items() if eta<=midpoint),F(0))
            integral+=mass*(1/left-1/right)
        assert error==sum(nonzero.values(),F(0))/W+integral/(2*W)
        windows.append(dict(W=W,error=error,main=c0,certified=error<c0,
                            minimum_actual_sum=min(actual),lower_bound=lower))
    return dict(shape=list(shape),A=sorted(A),target=g,cycle=cycle,spectrum=dict(spectrum),
                orbit_mean_multiplier=c0,coset_collision_pairs=pairs,windows=windows,
                f_on_cycle=[f[x] for x in cycle])


def encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {str(k):encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v) for v in value]
    return value
