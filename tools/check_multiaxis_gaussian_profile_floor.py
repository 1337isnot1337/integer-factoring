#!/usr/bin/env python3
"""Certify the shared-minor Gaussian pressure floors through eight axes.

All enclosures use integers at a fixed binary scale. No floating-point
optimizer or transcendental library function is used by the certificate.
The dimension-nine witness satisfies the relaxed budgets only; it is not
a polynomial kernel or a new factoring algorithm.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
import hashlib
import json

BITS=128
S=1<<BITS

def ceildiv(a,b):
    return -(-a//b)

@lru_cache(None)
def exp_negative(x):
    assert x>=0
    num,den=x.numerator,x.denominator
    h=0
    while 2*num>den:
        den*=2; h+=1
    # Even alternating Taylor sum is an upper bound, next odd sum lower.
    term_lo=term_hi=S
    sum_lo=sum_hi=S
    for j in range(1,65):
        term_lo=term_lo*num//(den*j)
        term_hi=ceildiv(term_hi*num,den*j)
        if j&1:
            sum_lo-=term_hi; sum_hi-=term_lo
        else:
            sum_lo+=term_lo; sum_hi+=term_hi
    next_hi=ceildiv(term_hi*num,den*65)
    lo=max(0,sum_lo-next_hi); hi=min(S,sum_hi)
    for _ in range(h):
        lo=lo*lo//S; hi=ceildiv(hi*hi,S)
    return lo,hi

def log_one_plus(lo,hi):
    assert 0<=lo<=hi<=S
    zl=lo*S//(2*S+lo)
    zu=ceildiv(hi*S,2*S+hi)
    z2l=zl*zl//S; z2u=ceildiv(zu*zu,S)
    pl,pu=zl,zu
    sl=su=0
    for j in range(64):
        sl+=pl//(2*j+1)
        su+=ceildiv(pu,2*j+1)
        pl=pl*z2l//S; pu=ceildiv(pu*z2u,S)
    # The true z<=1/3. The remaining logarithm series is <3^(-128).
    return 2*sl,2*su+ceildiv(S,3**128)

@lru_cache(None)
def softplus_negative(x):
    return log_one_plus(*exp_negative(x))

def log_positive_fraction(value):
    assert value>0
    k=0
    while value<1:
        value*=2
        k-=1
    while value>=2:
        value/=2
        k+=1
    u=value-1
    lo,hi=log_one_plus(u.numerator*S//u.denominator,
                       ceildiv(u.numerator*S,u.denominator))
    two_lo,two_hi=log_one_plus(S,S)
    if k>=0:
        return lo+k*two_lo,hi+k*two_hi
    return lo+k*two_hi,hi+k*two_lo

def profile(d):
    ans=[(1,Fraction(1,d))]
    ans.extend((1<<(d-q),Fraction(2,q-1)-Fraction(1,q)) for q in range(d,1,-1))
    assert sum(c for c,_ in ans)==1<<(d-1)
    assert sum(c*u for c,u in ans)==1<<(d-1)
    return ans

def certificate(d,rho):
    m=1<<(d-1); prof=profile(d)
    @lru_cache(None)
    def lhs_lo(x):
        return rho.denominator*sum(c*softplus_negative(x*u)[0] for c,u in prof)
    @lru_cache(None)
    def rhs_hi(x):
        return m*rho.numerator*softplus_negative(x/(d*rho))[1]
    tail_u=32*(1-rho)/(d*rho)
    tail_sum=sum((tail_u**i)/factorial(i) for i in range(33))
    assert tail_sum>2*m*rho
    stack=[(0,0)]; leaves=[]; gaps=[]
    while stack:
        depth,index=stack.pop()
        a=Fraction(32*index,1<<depth); b=Fraction(32*(index+1),1<<depth)
        gap=lhs_lo(b)-rhs_hi(a)
        if gap>0:
            leaves.append((depth,index));gaps.append(gap)
        else:
            assert depth<24,(d,depth,index)
            stack.extend(((depth+1,2*index+1),(depth+1,2*index)))
    depth=max(d for d,_ in leaves)
    prev=0
    for dep,idx in leaves:
        assert idx<<(depth-dep)==prev
        prev=(idx+1)<<(depth-dep)
    assert prev==1<<depth
    assert d*rho.numerator**2>rho.denominator**2
    packed=json.dumps(leaves,separators=(',',':')).encode()
    return dict(dimension=d,rho=[rho.numerator,rho.denominator],
                profile=[[c,u.numerator,u.denominator] for c,u in prof],
                threshold_integer_gap=d*rho.numerator**2-rho.denominator**2,
                leaves=leaves,leaf_count=len(leaves),depth=depth,
                minimum_interval_gap=min(gaps),
                cover_sha256=hashlib.sha256(packed).hexdigest(),
                tail_exponent=[tail_u.numerator,tail_u.denominator],
                tail_taylor_terms=33)

def independent_enclosure_controls():
    # Direct rational Taylor sums use no fixed-point recurrences or scaling.
    exp_cases=[Fraction(0),Fraction(1,8),Fraction(1,3),Fraction(1,2),
               Fraction(1),Fraction(2),Fraction(7,3),Fraction(8),
               Fraction(17,2),Fraction(32)]
    for x in exp_cases:
        term=Fraction(1); upper=term
        for j in range(1,513):
            term*=-x/j
            upper+=term
        lower=upper-term*x/513
        lo,hi=exp_negative(x)
        assert Fraction(lo,S)<=lower<=upper<=Fraction(hi,S)
    log_cases=[Fraction(0),Fraction(1,17),Fraction(1,7),Fraction(1,2),
               Fraction(2,3),Fraction(1)]
    for u in log_cases:
        z=u/(2+u)
        lower=2*sum(z**(2*j+1)/(2*j+1) for j in range(128))
        upper=lower+2*z**257/(257*(1-z*z))
        lo,hi=log_one_plus(u.numerator*S//u.denominator,
                           ceildiv(u.numerator*S,u.denominator))
        assert Fraction(lo,S)<=lower<=upper<=Fraction(hi,S)
    prefixes=0
    for d in range(2,13):
        expanded=[u for count,u in profile(d) for _ in range(count)]
        assert expanded==sorted(expanded)
        for q in range(1,d+1):
            j=1<<(d-q)
            assert sum(expanded[:j])==Fraction(j,q)
            prefixes+=1
    return dict(direct_rational_exp_controls=len(exp_cases),
                direct_rational_log_controls=len(log_cases),
                exact_profile_prefix_controls=prefixes)

def dimension_nine_relaxation():
    d=9; rho=Fraction(1,3)
    prof=profile(d); m=1<<(d-1)
    degree=2520
    integer_profile=[[c,int(degree*u)] for c,u in prof]
    assert all(degree*u==int(degree*u) for _,u in prof)
    t=Fraction(447,448); p=Fraction(13,100)
    log_z_hi=0
    for count,order in integer_profile:
        power=t**order
        bounds=(power.numerator*S//power.denominator,
                ceildiv(power.numerator*S,power.denominator))
        log_z_hi+=count*log_one_plus(*bounds)[1]
    log_t_lo,_=log_positive_fraction(t)
    _,log_p_hi=log_positive_fraction(p)
    _,log_one_minus_p_hi=log_positive_fraction(1-p)
    entropy_lo=-p*log_p_hi-(1-p)*log_one_minus_p_hi
    numerator_hi=log_z_hi-Fraction(m*degree,d)*p*log_t_lo
    gap=m*rho*entropy_lo-numerator_hi
    assert gap>0
    return dict(scope='Relaxed shared-minor budgets only; no realizing kernel',
                dimension=d,degree=degree,integer_profile=integer_profile,
                rational_scale=[t.numerator,t.denominator],rho=[1,3],
                physical_frequency=[p.numerator,p.denominator],
                strict_negative_pressure_gap=[gap.numerator,gap.denominator])

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    controls=independent_enclosure_controls()
    results=[]
    for d,r in [(2,Fraction(3,4)),(3,Fraction(3,5)),(4,Fraction(21,40)),(5,Fraction(9,20)),(6,Fraction(5,12)),(7,Fraction(19,50)),(8,Fraction(177,500))]:
        v=certificate(d,r);results.append(v)
    result=dict(status='passed',scope='Common-scale unweighted Gaussian profiles with a common degree cap',
                arithmetic_bits=BITS,compact_interval=[0,32],
                controls=controls,certificates=results,
                dimension_nine_relaxation=dimension_nine_relaxation(),
                factoring_score='1/2; unchanged')
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    summary=dict(result)
    summary['certificates']=[{k:v for k,v in result.items() if k not in ('leaves','profile')}
                             for result in results]
    print(json.dumps(summary,indent=2))
