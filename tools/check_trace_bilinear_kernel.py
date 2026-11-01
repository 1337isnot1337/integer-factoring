#!/usr/bin/env python3
"""Binary trace-kernel ranks, Vandermonde coefficients, and a large symbolic witness."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path


def polynomial_remainder(a,b):
    while a and a.bit_length()>=b.bit_length():
        a^=b<<(a.bit_length()-b.bit_length())
    return a


def irreducibility_by_division(p):
    degree=p.bit_length()-1
    assert p&1
    checked=0
    for e in range(1,degree//2+1):
        for middle in range(1<<(e-1)):
            divisor=(1<<e)|(middle<<1)|1
            assert polynomial_remainder(p,divisor), (p,divisor)
            checked+=1
    return checked


class Field:
    def __init__(self,polynomial):
        self.polynomial=polynomial
        self.degree=polynomial.bit_length()-1
        self.order=1<<self.degree
        self.alpha=polynomial_remainder(2,polynomial)

    def mul(self,a,b):
        answer=0
        while b:
            if b&1:
                answer^=a
            b>>=1
            a<<=1
            if a&self.order:
                a^=self.polynomial
        return answer

    def power(self,a,n):
        answer=1
        while n:
            if n&1:
                answer=self.mul(answer,a)
            a=self.mul(a,a)
            n>>=1
        return answer

    def trace(self,a):
        answer=0
        for _ in range(self.degree):
            answer^=a
            a=self.mul(a,a)
        assert answer in (0,1)
        return answer


def rank(rows):
    basis={}
    for v in rows:
        while v:
            pivot=v.bit_length()-1
            if pivot in basis:
                v^=basis[pivot]
            else:
                basis[pivot]=v
                break
    return len(basis)


def convolution_rank(values):
    d=len(values)
    rows=[sum(values[x^y]<<y for y in range(d)) for x in range(d)]
    return rank(rows)


def trace_polynomial(field,k):
    modulus=field.order-1
    polynomial={0:1}
    for j in range(k):
        factor=[(1<<i,field.power(field.alpha,j*(1<<i)))
                for i in range(field.degree)]
        result={}
        for a,c in polynomial.items():
            for b,e in factor:
                exponent=(a+b)%modulus
                result[exponent]=result.get(exponent,0)^field.mul(c,e)
        polynomial={a:c for a,c in result.items() if c}
    return polynomial


def profile(field,k):
    t,d=field.degree,field.order
    alpha_powers=[field.power(field.alpha,j) for j in range(k)]
    values=[]
    for z in range(d):
        value=1
        for a in alpha_powers:
            value&=field.trace(field.mul(a,z))
        values.append(value)
    kernel=[[values[field.mul(s,w)] for w in range(d)] for s in range(d)]
    ordinary_rank=rank([sum(v<<j for j,v in enumerate(row)) for row in kernel])
    row_blocks=[convolution_rank(row) for row in kernel]
    column_blocks=[convolution_rank([kernel[s][w] for s in range(d)]) for w in range(d)]
    assert row_blocks==[0]+[1<<k]*(d-1)
    assert column_blocks==row_blocks
    polynomial=trace_polynomial(field,k)
    assert ordinary_rank==len(polynomial)
    assert ordinary_rank>=math.comb(t,k)
    for z in range(1,d):
        actual=0
        for exponent,coefficient in polynomial.items():
            actual^=field.mul(coefficient,field.power(z,exponent))
        assert actual==values[z]
    conjugates=[field.power(field.alpha,1<<i) for i in range(t)]
    checked=0
    for subset in itertools.combinations(range(t),k):
        expected=1
        for i,j in itertools.combinations(subset,2):
            expected=field.mul(expected,conjugates[i]^conjugates[j])
        exponent=sum(1<<i for i in subset)
        assert expected and polynomial[exponent]==expected
        checked+=1
    if t<=3:
        row=[0]*(d*d)
        col=[0]*(d*d)
        for x,y,u,v in itertools.product(range(d),repeat=4):
            if x^y^u^v==0 and kernel[x^y][x^u]:
                row[x*d+y]|=1<<(u*d+v)
                col[x*d+u]|=1<<(y*d+v)
        assert rank(row)==sum(row_blocks)
        assert rank(col)==sum(column_blocks)
    # Independent coordinate monomial controls, with explicit support loss.
    a=min(2,t)
    b=min(1,t)
    controlled=[[value if (s&((1<<a)-1))==((1<<a)-1)
                 and (w&((1<<b)-1))==((1<<b)-1) else 0
                 for w,value in enumerate(row)] for s,row in enumerate(kernel)]
    controlled_r=sum(convolution_rank(row) for row in controlled)
    controlled_c=sum(convolution_rank([controlled[s][w] for s in range(d)]) for w in range(d))
    assert controlled_r*(1<<a)<=d*(1<<(k+b))
    assert controlled_c*(1<<b)<=d*(1<<(k+a))
    assert controlled_r*controlled_c<=(d*(1<<k))**2
    return {"t":t,"k":k,"field_polynomial":field.polynomial,
            "ordinary_kernel_rank":ordinary_rank,"binomial_rank_lower":math.comb(t,k),
            "row_cut_rank":sum(row_blocks),"column_cut_rank":sum(column_blocks),
            "top_coefficients_checked":checked,
            "controlled_prefactor_degrees":[a,b],
            "controlled_cut_ranks":[controlled_r,controlled_c],
            "kernel_sha256":hashlib.sha256(bytes(v for row in kernel for v in row)).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    fields={2:7,3:11,4:19,5:37,6:67,7:131}
    profiles=[]
    for t,p in fields.items():
        irreducibility_by_division(p)
        field=Field(p)
        for k in range(1,min(t,4)):
            result=profile(field,k)
            profiles.append(result)
            print(json.dumps(result),flush=True)
    t,k=33,11
    polynomial=(1<<33)|(1<<10)|1
    divisors=irreducibility_by_division(polynomial)
    assert divisors==65535
    field=Field(polynomial)
    conjugates=[field.power(field.alpha,1<<i) for i in range(t)]
    assert len(set(conjugates))==t
    assert field.power(field.alpha,1<<t)==field.alpha
    assert field.trace(1)==1
    d=1<<t
    cut_rank=(d-1)*(1<<k)
    kernel_rank_lower=math.comb(t,k)
    assert kernel_rank_lower**2>1<<(t+2*k)
    assert cut_rank**2<d*kernel_rank_lower**2
    witness={"t":t,"k":k,"label_dimension":d,
             "irreducible_polynomial":polynomial,"polynomial_hex":hex(polynomial),
             "candidate_divisors_checked":divisors,"distinct_conjugates":len(set(conjugates)),
             "exact_row_cut_rank":cut_rank,"exact_column_cut_rank":cut_rank,
             "ordinary_kernel_rank_lower":kernel_rank_lower,
             "kernel_rank_lower_squared":kernel_rank_lower**2,
             "binomial_comparison_power":1<<(t+2*k),
             "exact_cut_rank_product":cut_rank**2,
             "rejected_bound_lower":d*kernel_rank_lower**2,
             "full_large_tensor_materialized":False}
    report={"status":"passed","profiles":profiles,"large_symbolic_witness":witness,
            "scope":"Binary candidate pieces and a rank-inequality counterexample; no crossing synthesis."}
    if args.output:
        temporary=args.output.with_suffix(args.output.suffix+".tmp")
        with temporary.open("w") as out:
            json.dump(report,out,indent=2)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        temporary.replace(args.output)
    print(json.dumps(witness,indent=2))


if __name__=="__main__":
    main()
