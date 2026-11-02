#!/usr/bin/env python3
"""Generic distinct-color pressure obstruction, with no large-q matrix screen."""
import argparse
from hashlib import sha256
from fractions import Fraction
from itertools import permutations,product
import json
from math import comb,factorial
from pathlib import Path
import cube_homogeneous_gaussian as cube


def colored_entry(left,right,axis,gaussian):
    first,second=cube.axis_parts(axis);masks=[0]*6
    for labels,vertices in ((left,first),(right,second)):
        for c,v in zip(labels,vertices):masks[c]|=1<<v
    value=cube.ONE
    for mask in masks[1:]:
        value=cube.multiply(value,gaussian[mask])
        if not value:break
    return value


def structural_control():
    kernel=cube.kernel();gaussian=cube.gaussian_tensor(kernel)
    assert gaussian[0]==cube.ONE
    assert all(not g for mask,g in enumerate(gaussian) if mask.bit_count()%2)
    assert all(sum(term)==3*mask.bit_count() for mask,g in enumerate(gaussian) for term in g)
    assert kernel[0][7]==cube.monomial(2,2,2)
    words=list(permutations((1,2,3,4)));records=[]
    for axis in range(3):
        first,second=cube.axis_parts(axis)
        assert all(cube.valuation(kernel[a][b],axis)==2 for a in first for b in second)
        endpoint_checks=0
        for left in words:
            for right in words:
                actual=colored_entry(left,right,axis,gaussian)
                expected=cube.ONE
                for i,color in enumerate(left):
                    expected=cube.multiply(expected,kernel[first[i]][second[right.index(color)]])
                assert actual==expected and cube.valuation(actual,axis)==8
                endpoint=all(right[second.index(v^7)]==left[i] for i,v in enumerate(first))
                assert len(actual)%2==endpoint;endpoint_checks+=1
        left=words[0]
        antipodal=tuple(left[first.index(v^7)] for v in second)
        assert colored_entry(left,antipodal,axis,gaussian)==cube.monomial(8,8,8)
        closed_checks=0
        for right in product(range(6),repeat=4):
            if sorted(right)!=[1,2,3,4]:
                assert not colored_entry(left,right,axis,gaussian)
                closed_checks+=1
        minor_checks=0
        for c,d in product(range(6),repeat=2):
            left=[0]*4;right=[0]*4;left[first.index(0)]=c;right[second.index(7)]=d
            actual=colored_entry(left,right,axis,gaussian)
            expected=cube.ONE if c==d==0 else kernel[0][7] if c==d else cube.ZERO
            assert actual==expected;minor_checks+=1
        records.append(dict(axis=axis,sector_order=24,sector_endpoint_checks=endpoint_checks,closed_sector_zero_checks=closed_checks,complement_minor_checks=minor_checks,smallest_sector_order=8))
    return records


def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def power(a,k):
    c=[1]
    for _ in range(k):c=mul(c,a)
    return c


def pressure_certificate():
    # d=255+r, A(r)=24*binom(d,4). P=24^5*(L_d(s^3)^5-(1+d*s^10)^12).
    A=[1]
    for c in (252,253,254,255):A=mul(A,[c,1])
    d=[255,1];P=[[0]*121 for _ in range(21)]
    for a in range(6):
        for b in range(6-a):
            c=5-a-b;degree=6*a+24*b
            co=factorial(5)//(factorial(a)*factorial(b)*factorial(c))*24**(a+c)
            for j,v in enumerate(mul(power(d,a),power(A,b))):P[j][degree]+=co*v
    for b in range(13):
        for j,v in enumerate(power(d,b)):P[j][10*b]-=24**5*comb(12,b)*v
    assert all(not any(p[:6]) for p in P)
    transformed=[]
    for p in P:
        q=p[6:]
        transformed.append([sum(q[j]*comb(114-j,k-j) for j in range(k+1)) for k in range(115)])
    for parameter,point in ((0,Fraction(1,2)),(1,Fraction(2,3)),(17,Fraction(3,5))):
        colors=255+parameter
        explicit=24**5*((1+colors*point**6+comb(colors,4)*point**24)**5-(1+colors*point**10)**12)
        expanded=sum(v*parameter**j*point**k for j,row in enumerate(P) for k,v in enumerate(row))
        assert explicit==expanded
        u=point/(1-point)
        check=sum(v*parameter**j*u**k for j,row in enumerate(transformed) for k,v in enumerate(row))
        assert check==(1+u)**114*explicit/point**6
    assert all(v>=0 for row in transformed for v in row)
    assert min(transformed[0])>0
    assert 3*3**2>5**2
    return dict(base_nonzero_colors=255,degree_in_parameter=20,degree_in_s=120,divided_power=6,mobius_degree=114,coefficients_checked=2415,independent_polynomial_identity_checks=3,constant_parameter_min=min(transformed[0]),nonnegative_coefficients=2415,coefficient_sha256=sha256(json.dumps(transformed,separators=(',',':')).encode()).hexdigest(),rho_floor=[3,5],threshold_check=dict(left=27,right=25),scope='All integer q>=256, all common physical types, common-scale unweighted three-axis Smith-profile pressure.')


def main():
    controls=structural_control();certificate=pressure_certificate()
    result=dict(scope='Generic direct-sector and complementary-minor identities plus a uniform integer pressure certificate. No specialized rank extrapolation, large-q matrix screen, stronger factoring exponent, or independent mathematical review.',original_color_checkpoint='1772970a7dbaab83eeb0607014a1111e46718d86',axes=controls,certificate=certificate,kernel_source_sha256=sha256(Path(cube.__file__).read_bytes()).hexdigest())
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(axes=len(controls),sector_endpoint_checks=sum(r['sector_endpoint_checks'] for r in controls),closed_sector_zero_checks=sum(r['closed_sector_zero_checks'] for r in controls),minor_checks=sum(r['complement_minor_checks'] for r in controls),coefficients_checked=2415,rho_floor=[3,5],minimum_q=256)))


if __name__=='__main__':main()
