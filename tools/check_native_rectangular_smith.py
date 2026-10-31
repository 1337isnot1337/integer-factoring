#!/usr/bin/env python3
"""Exact profile replays and API controls for the optional rectangular engine."""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path
import random
from native_rectangular_smith import RectangularOracle
from weighted_smith_oracle import WeightedSmithOracle


def matrix(values,A,B,cut):
    out=[[0]*(A*B) for _ in range(A*B)]
    for x,y,u,v in product(range(A),range(B),range(B),range(A)):
        out[B*x+(y if cut==0 else u)][A*(u if cut==0 else y)+v]=values[x+A*(y+B*(u+B*v))]
    return out


def exact_orders(values,A,B,params,cut,denominator):
    alph=[A,B,B,A];offsets=[];stride=0
    for q in alph:offsets.append(stride);stride+=q-1
    def height(side,leg,a):return params[side*stride+offsets[leg]+a-1] if a else 0
    orders=[]
    for side in range(2):
        selected=cut if side==0 else 1-cut
        row_legs=(0,1) if selected==0 else(0,2)
        col_legs=(2,3) if selected==0 else(1,3)
        rr=[height(side,row_legs[0],x)+height(side,row_legs[1],y) for x in range(A) for y in range(B)]
        cc=[height(side,col_legs[0],u)+height(side,col_legs[1],v) for u in range(B) for v in range(A)]
        oracle=WeightedSmithOracle(matrix(values,A,B,selected),infinity=(side==0))
        current,_=oracle.profile(rr,cc,params[-1])
        if side==0:current=[v+denominator*params[-1] for v in current]
        orders.extend(current)
    return orders


def logsum(values):
    top=max(values);return top+math.log(sum(math.exp(v-top) for v in values))


def check():
    rng=random.Random(941);records=[];maximum_error=0.;last=None
    for A,B,d in ((2,2,0),(2,3,0),(3,2,1),(3,3,2)):
        values=[0]*(A*A*B*B)
        for x,y,u,v in product(range(A),range(B),range(B),range(A)):
            p=rng.randrange(16);values[x+A*(y+B*(u+B*v))]=int(x==v and y==u)^(p<<1)^(p<<2)
        oracle=RectangularOracle(values,A,B,d,True)
        params=[Fraction(rng.randrange(-150,151),509) for _ in range(oracle.parameters-1)]+[Fraction(1,7)]
        px=[Fraction(i+1,A*(A+1)//2) for i in range(A)];py=[Fraction(i+1,B*(B+1)//2) for i in range(B)]
        for cut in (0,1):
            expected=exact_orders(values,A,B,params,cut,d)
            cost,gradient,actual=oracle.evaluate(params,px,py,cut,True)
            error=max(abs(float(a)-b) for a,b in zip(expected,actual));maximum_error=max(maximum_error,error);assert error<2e-10
            mean=0.
            for side in range(2):
                for leg in range(4):
                    p=px if leg in(0,3) else py
                    mean+=sum(float(p[a]*params[side*oracle.stride+oracle.offsets[leg]+a-1]) for a in range(1,oracle.alphabets[leg]))
            exact_cost=logsum([-float(v) for v in expected[:oracle.n]])+logsum([-float(v) for v in expected[oracle.n:]])+mean
            assert abs(exact_cost-cost)<2e-10
            for j in (0,len(params)//2,len(params)-1):
                xp=list(map(float,params));xm=xp.copy();xp[j]+=1e-7;xm[j]-=1e-7
                difference=(oracle.evaluate(xp,px,py,cut)[0]-oracle.evaluate(xm,px,py,cut)[0])/2e-7
                assert abs(difference-gradient[j])<2e-6,(A,B,cut,j,difference,gradient[j])
            records.append({'A':A,'B':B,'denominator_degree':d,'cut':cut,'exact_orders':[str(v) for v in expected]})
        if d==0:
            rho,gx,gy,_,_,_=oracle.types(params)
            cost,_=oracle.evaluate(params,gx,gy)
            entropy=-sum(p*math.log(p) for p in gx+gy)
            assert abs(cost/(2*entropy)-rho)<2e-10
        last=(oracle,values,params,px,py)
    oracle,values,params,px,py=last
    newer=RectangularOracle(values,oracle.A,oracle.B,oracle.denominator_degree,True)
    try:oracle.evaluate(params,px,py)
    except RuntimeError:pass
    else:raise AssertionError('Inactive seed was reused without a generation error.')
    oracle.activate();oracle.evaluate(params,px,py)
    try:RectangularOracle([1<<512]+[0]*15,2,2)
    except ValueError:pass
    else:raise AssertionError('An oversized polynomial was accepted.')
    try:RectangularOracle([0]*16,2,2,endpoint_check=True)
    except ArithmeticError:pass
    else:raise AssertionError('A wrong strict endpoint was accepted.')
    try:oracle.evaluate(params,px,py)
    except RuntimeError:pass
    else:raise AssertionError('A failed seed load did not invalidate the prior generation.')
    oracle.activate();oracle.evaluate(params,px,py)
    zero=RectangularOracle([0]*16,2,2,endpoint_check=False)
    try:zero.evaluate([0.]*8+[.1],[.5,.5],[.5,.5])
    except ArithmeticError:pass
    else:raise AssertionError('An unresolved singular profile was accepted.')
    return {'status':'pass','exact_fraction_free_replays':len(records),
            'maximum_order_error':maximum_error,'paired_gibbs_identity':True,'gradient_coordinates_checked':24,
            'active_seed_guard':True,'failed_load_invalidates_prior_seed':True,
            'oversized_input_rejected':True,'wrong_endpoint_rejected':True,
            'singular_profile_rejected':True,'replays':records,
            'scope':'Native weights are exploratory; exact controls do not certify arbitrary future queries.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args();result=check();text=json.dumps(result,indent=2)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text+'\n')
    print(text)

if __name__=='__main__':main()
