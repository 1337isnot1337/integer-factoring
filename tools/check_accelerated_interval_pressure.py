#!/usr/bin/env python3
"""Exact logarithmic pressure certificates without visiting every mode.

An integrated finite geometric identity gives an Euler expansion for
sum_j log(1+x*q**j), with a positive explicit remainder. Every computed
endpoint is an integer multiple of 2**(-precision). No floating-point
arithmetic participates in a certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import time


class FixedIntervals:
    def __init__(self, precision):
        if precision<32:raise ValueError('Use at least 32 fractional bits.')
        self.precision=precision;self.scale=1<<precision
        self.zero=(0,0);self.one=(self.scale,self.scale)

    @staticmethod
    def ceiling(a,b):return -((-a)//b)

    def rational(self,numerator,denominator=1):
        if denominator<=0:raise ValueError('The denominator must be positive.')
        numerator*=self.scale
        return numerator//denominator,self.ceiling(numerator,denominator)

    @staticmethod
    def add(a,b):return a[0]+b[0],a[1]+b[1]

    @staticmethod
    def subtract(a,b):return a[0]-b[1],a[1]-b[0]

    def multiply(self,a,b):
        products=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
        return min(products)//self.scale,self.ceiling(max(products),self.scale)

    def divide_positive(self,a,b):
        if a[0]<0 or b[0]<=0:raise ArithmeticError('Positive interval division requires resolved signs.')
        return a[0]*self.scale//b[1],self.ceiling(a[1]*self.scale,b[0])

    @staticmethod
    def integer_multiple(a,coefficient):
        return (a[0]*coefficient,a[1]*coefficient) if coefficient>=0 else (a[1]*coefficient,a[0]*coefficient)

    def integer_division(self,a,denominator):
        if denominator<=0:raise ValueError('The denominator must be positive.')
        return a[0]//denominator,self.ceiling(a[1],denominator)

    def power(self,value,exponent):
        if exponent<0:raise ValueError('Use a nonnegative power.')
        result=self.one
        while exponent:
            if exponent&1:result=self.multiply(result,value)
            exponent>>=1
            if exponent:value=self.multiply(value,value)
        return result

    def fractions(self,value):return tuple(Fraction(endpoint,self.scale) for endpoint in value)


def geometric_moments(arithmetic,x,q,count,terms):
    """A_k=sum_(j<count) (x*q**j)**k/k, k=1,...,terms."""
    if count<0 or terms<1:raise ValueError('Invalid product size or series length.')
    if count==0:return [arithmetic.zero]*terms
    if not(0<=x[0]<=x[1]<=arithmetic.scale and 0<=q[0]<=q[1]<arithmetic.scale):
        raise ValueError('Require 0<=x<=1 and a resolved 0<=q<1.')
    q_count=arithmetic.power(q,count)
    x_power,q_power,q_count_power=x,q,q_count
    result=[]
    for k in range(1,terms+1):
        numerator=arithmetic.multiply(x_power,arithmetic.subtract(arithmetic.one,q_count_power))
        denominator=arithmetic.integer_multiple(arithmetic.subtract(arithmetic.one,q_power),k)
        result.append(arithmetic.divide_positive(numerator,denominator))
        x_power=arithmetic.multiply(x_power,x)
        q_power=arithmetic.multiply(q_power,q)
        q_count_power=arithmetic.multiply(q_count_power,q_count)
    return result


def euler_log_product(arithmetic,x,q,count,terms,return_differences=False):
    """Enclose the log product; the omitted tail is at most A_1/2**terms."""
    moments=geometric_moments(arithmetic,x,q,count,terms)
    if not moments:return arithmetic.zero
    first_upper=moments[0][1];differences=[];total=arithmetic.zero
    for order in range(terms):
        # D_m=sum_j integral_0^(x*q^j)(1-t)^m dt lies in[0,A_1].
        value=max(0,moments[0][0]),min(first_upper,moments[0][1])
        if value[0]>value[1]:raise ArithmeticError('The exact Euler positivity bounds became inconsistent.')
        differences.append(value)
        total=arithmetic.add(total,arithmetic.integer_division(value,1<<(order+1)))
        moments=[arithmetic.subtract(a,b) for a,b in zip(moments,moments[1:])]
    tail=arithmetic.ceiling(first_upper,1<<terms)
    interval=total[0],total[1]+tail
    return (interval,differences) if return_differences else interval


def log_pressure_interval(n,a,b,p,q,terms=None,precision=None):
    if min(n,a,b,q)<=0 or p<=q:raise ValueError('Supply a positive interval witness with r>1.')
    if n%2 and p>2*q:raise ValueError('For odd intervals this log-root formula requires 1<r<=2.')
    size=max(n,a,b,p,q).bit_length()
    if terms is None:terms=2*max(n,a,b).bit_length()+80
    if precision is None:precision=4*terms+2*size+128
    arithmetic=FixedIntervals(precision)
    inverse_root=arithmetic.rational(q,p)
    numerator=euler_log_product(arithmetic,arithmetic.power(inverse_root,a),
                                 arithmetic.power(inverse_root,2*a),n,terms)
    total=arithmetic.integer_multiple(numerator,4*b)
    if n%2:
        denominator=euler_log_product(arithmetic,arithmetic.power(inverse_root,8*b),
                                       arithmetic.power(inverse_root,8*b),(n-1)//2,terms)
        log_root=euler_log_product(arithmetic,arithmetic.rational(p-q,q),
                                    arithmetic.zero,1,terms)
        log_two=euler_log_product(arithmetic,arithmetic.one,arithmetic.zero,1,terms)
        total=arithmetic.add(total,arithmetic.integer_multiple(log_root,a*b))
        total=arithmetic.subtract(total,arithmetic.integer_multiple(log_two,a))
    else:
        denominator=euler_log_product(arithmetic,arithmetic.power(inverse_root,4*b),
                                       arithmetic.power(inverse_root,8*b),n//2,terms)
    total=arithmetic.subtract(total,arithmetic.integer_multiple(denominator,2*a))
    return arithmetic,total,terms


def exact_controls():
    checks=0
    for precision in (48,127):
        arithmetic=FixedIntervals(precision)
        for x,q,n in ((Fraction(1),Fraction(0),1),(Fraction(1,2),Fraction(2,3),3),
                      (Fraction(7,9),Fraction(4,5),8),(Fraction(1,17),Fraction(31,32),5)):
            terms=16
            ex=arithmetic.rational(x.numerator,x.denominator)
            eq=arithmetic.rational(q.numerator,q.denominator)
            moments=geometric_moments(arithmetic,ex,eq,n,terms)
            _,differences=euler_log_product(arithmetic,ex,eq,n,terms,True)
            for k,value in enumerate(moments,1):
                exact=sum((x*q**j)**k/Fraction(k) for j in range(n))
                lo,hi=arithmetic.fractions(value);assert lo<=exact<=hi;checks+=1
            for m,value in enumerate(differences):
                exact=sum((1-(1-x*q**j)**(m+1))/Fraction(m+1) for j in range(n))
                lo,hi=arithmetic.fractions(value);assert lo<=exact<=hi;checks+=1
    from check_large_interval_pressure import pressure_ratio,BinaryIntervals
    sign_cases=[]
    for witness in ((4,7,3,33,32),(7,7,3,33,32),(63,17,6,177,176),(511,99,35,5756,5755)):
        arithmetic,interval,terms=log_pressure_interval(*witness)
        low,high=BinaryIntervals.fraction_bounds(pressure_ratio(*witness))
        if high<1:assert interval[1]<0;sign=-1
        elif low>1:assert interval[0]>0;sign=1
        else:raise AssertionError('The independent direct product did not resolve its sign.')
        sign_cases.append({'witness':list(witness),'exact_sign':sign,'terms':terms})
    return {'independent_fraction_moments_and_differences':checks,
            'independent_direct_product_signs':sign_cases}


def pell_witness(level):
    if level<5:raise ValueError('The simple root rule is intended for Pell levels at least 5.')
    a,b=3,1
    for _ in range(1,level):a,b=3*a+8*b,a+3*b
    assert a*a-8*b*b==1
    bits=a.bit_length()
    while 3*((1<<bits)-1)<a*bits:bits+=1
    n=(1<<bits)-1
    root_denominator=a*n//bits
    return n,a,b,root_denominator+1,root_denominator


def certify(witness,level=None):
    started=time.monotonic();n,a,b,p,q=witness
    arithmetic,interval,terms=log_pressure_interval(*witness)
    assert interval[0]<=interval[1]<0,(witness,interval)
    rho=Fraction(a,4*b)
    result={'modes':n,'mode_count_binary_exponent':(n+1).bit_length()-1 if n&(n+1)==0 else None,
            'R':[a,b],'rho_upper':[rho.numerator,rho.denominator],'root':[p,q],
            's_root_power':a,'euler_terms':terms,'fractional_bits':arithmetic.precision,
            'log_pressure_ratio_interval':{'lower_numerator':str(interval[0]),
                                           'upper_numerator':str(interval[1]),
                                           'denominator_binary_exponent':arithmetic.precision},
            'strict_pressure_inequality':True,'post_half_target_reached':rho*rho<Fraction(1,2),
            'rational_mode_probabilities':{'formula':'r^(-b*(8*i-4*(n-1)))/(1+r^(-b*(8*i-4*(n-1))))',
                                           'mode_range':[0,n-1],'r':[p,q],'b':b},
            'elapsed_seconds':time.monotonic()-started}
    if level is not None:result['pell_level']=level
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--levels',type=int,nargs='+',default=[5,6,7,8,12,20,40,80])
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();controls=exact_controls();certificates=[]
    for level in args.levels:
        certificate=certify(pell_witness(level),level);certificates.append(certificate)
        print(json.dumps({'level':level,'mode_bits':certificate['mode_count_binary_exponent'],
                          'rho_upper':certificate['rho_upper'],'seconds':certificate['elapsed_seconds']}),flush=True)
    result={'status':'exact finite interval witnesses; no improved factoring exponent',
            'method':'integrated Euler expansion with outward fixed-point arithmetic',
            'controls':controls,'certificates':certificates,
            'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in (Path(__file__).name,'check_large_interval_pressure.py')}}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'certificates':len(certificates),'all_strict':True}))


if __name__=='__main__':main()
