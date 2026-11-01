#!/usr/bin/env python3
"""Exact all-scale certificate and finite controls for a fixed cube projection."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from cube_homogeneous_gaussian import (kernel,augmentation_crossblock,
    exact_monomial_smith,determinantal_orders,gaussian_tensor,
    gaussian_cut_factors,matrix_product,axis_matrix,serializable)

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/fixed-cube-projection-floor.json'


def multiply(a,b):
    result=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):result[i+j]+=x*y
    return result


def power(a,n):
    result=[1]
    for _ in range(n):result=multiply(result,a)
    return result


def binomial_polynomial(exponent):
    return [1]+[0]*(exponent-1)+[1]


def partition_inequality():
    left=multiply(power(binomial_polynomial(3),2),power(binomial_polynomial(9),2))
    right=power(binomial_polynomial(4),3)
    difference=[left[i]-(right[i] if i<len(right) else 0) for i in range(len(left))]
    assert difference[:3]==[0,0,0]
    q=difference[3:]
    assert q==[2,-3,0,1,0,-3,2,0,0,3,0,0,2,0,0,1,0,0,2,0,0,1]
    degree=64
    bernstein=[sum((Fraction(c*comb(k,j),comb(degree,j))
                    for j,c in enumerate(q) if j<=k),Fraction())
               for k in range(degree+1)]
    assert min(bernstein)>0
    # Independent reconstruction in the ordinary monomial basis.
    reconstruction=[Fraction() for _ in range(degree+1)]
    for k,b in enumerate(bernstein):
        for j in range(degree-k+1):
            reconstruction[k+j]+=b*comb(degree,k)*comb(degree-k,j)*((-1)**j)
    assert reconstruction==list(map(Fraction,q))+[Fraction()]*(degree+1-len(q))
    assert min(bernstein)==Fraction(1589593605121,28736878651740)
    return dict(identity='(1+u^3)^2(1+u^9)^2-(1+u^4)^3=u^3*Q(u)',
                Q=q,bernstein_degree=degree,
                bernstein_coefficients=[[x.numerator,x.denominator] for x in bernstein],
                minimum_bernstein=[min(bernstein).numerator,min(bernstein).denominator],
                ordinary_basis_reconstruction=True,
                positive_interval='0<=u<=1; reciprocal extension handles all r>1')


def projection_controls():
    records=[];minors=0;factorizations=0;endpoints=0
    for a,b in ((1,0),(1,1),(2,1),(3,1),(3,2),(2,0)):
        def substitute(polynomial):
            result=set()
            for i,j,k in polynomial:
                exponent=(a*i+b*k,a*j+(a-b)*k,0)
                if exponent in result:result.remove(exponent)
                else:result.add(exponent)
            return frozenset(result)
        projected=[[substitute(p) for p in row] for row in kernel()]
        tensor=gaussian_tensor(projected)
        for mask,polynomial in enumerate(tensor):
            expected=int(all(((mask>>v)&1)==((mask>>(v^7))&1) for v in range(8)))
            assert len(polynomial)%2==expected
            assert all(i+j+k==3*a*mask.bit_count() for i,j,k in polynomial)
            endpoints+=1
        profiles=[]
        for axis in (0,1):
            block=augmentation_crossblock(projected,axis)
            orders,count=determinantal_orders(block,axis);minors+=count
            exact=exact_monomial_smith(block,axis)
            s=min(3*a,2*a+2*(b if axis==0 else a-b))
            predicted=[s,6*a-s,6*a+s,12*a-s]
            assert orders==exact['orders']==predicted
            left,middle,right,cross=gaussian_cut_factors(projected,axis)
            assert matrix_product(matrix_product(left,middle),right)==axis_matrix(tensor,axis)
            factorizations+=1
            profiles.append(dict(axis=axis,orders=orders,determinantal_minors=count,
                                 exact_finite_Laurent=exact))
        records.append(dict(a=a,b=b,homogeneous_degree=6*a,profiles=profiles))
    return dict(records=records,determinantal_minors=minors,
                physical_factorizations=factorizations,endpoint_assignments=endpoints)


def report():
    sources=['tools/check_fixed_cube_projection_floor.py','tools/cube_homogeneous_gaussian.py']
    return serializable(dict(
        scope='All-scale common-type Smith-profile envelope floor rho>=3/4 for fixed homogeneous monomial projections of the existing cube; not a lower bound for other representations or factoring.',
        partition_certificate=partition_inequality(),projections=projection_controls(),
        target_excluded=dict(floor=[3,4],squared_floor=[9,16],target_squared=[1,2],
                             exact_gap=[1,16]),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],
        projections=len(result['projections']['records']),
        minors=result['projections']['determinantal_minors'],
        factorizations=result['projections']['physical_factorizations'],
        endpoints=result['projections']['endpoint_assignments'],
        bernstein_minimum=result['partition_certificate']['minimum_bernstein'],
        floor=result['target_excluded']['floor'])))


if __name__=='__main__':main()
