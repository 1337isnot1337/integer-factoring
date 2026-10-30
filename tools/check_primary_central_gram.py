#!/usr/bin/env python3
"""Exact checks for the primitive restriction of the central homogeneous piece.

This checks finite algebraic identities and rational ranks. The all-weight
argument and its restricted scope are in research/primary-central-gram.md.
No stronger asymptotic factoring exponent is asserted.
"""

from fractions import Fraction as F
from math import comb, factorial, prod
from functools import lru_cache
import json
from check_partition_channels import operator, partitions, partition_denominator
from check_homogeneous_duals import add_scaled, shift, x_operator


def rank_q(a):
    a=[[F(x) for x in r] for r in a]
    if not a: return 0
    rank=0
    for c in range(len(a[0])):
        pivot=next((r for r in range(rank,len(a)) if a[r][c]),None)
        if pivot is None: continue
        a[rank],a[pivot]=a[pivot],a[rank]
        lead=a[rank][c]
        a[rank][c:]=[v/lead for v in a[rank][c:]]
        for r in range(rank+1,len(a)):
            if not a[r][c]: continue
            mult=a[r][c]
            a[r][c:]=[v-mult*w for v,w in zip(a[r][c:],a[rank][c:])]
        rank+=1
        if rank==len(a): break
    return rank


def scaled(a,c): return {x:v*c for x,v in a.items() if v*c}


def project(j,a):
    if j<2: return a.copy()
    result=a.copy()
    for k in range(1,j):
        chi=j*(j-1)-k*(k-1)
        old=result
        result=x_operator(old)
        add_scaled(result,old,-chi)
        result=scaled(result,F(-1,chi))
    return result


@lru_cache(None)
def primitive_basis(h):
    basis=partitions(h)
    vectors=[]
    columns=[]
    for mu in basis:
        v=project(h,{mu:F(1)})
        if not v: continue
        col=[v.get(nu,0) for nu in basis]
        if rank_q(columns+[col])>len(columns):
            columns.append(col)
            vectors.append(v)
    assert all(not shift(v,-1) for v in vectors)
    expected=1 if h==1 else len(partitions(h))-len(partitions(h-1))
    assert len(vectors)==expected,(h,len(vectors),expected)
    return vectors


def act(a,y,energy):
    result={}
    for out in partitions(energy):
        v=sum(av*yv*operator(mu,out,nu) for mu,av in a.items() for nu,yv in y.items())
        if v: result[out]=v
    return result


def dot(a,b):
    return sum(v*b.get(mu,0)*partition_denominator(mu) for mu,v in a.items())


def signs(a):
    return {mu:v*(-1)**(sum(mu)-len(mu)) for mu,v in a.items()}


def weight_coeff(h,j):
    if j==0: return F(1)
    return sum(F((-1)**(j+n)*comb(2*h,j+n)*prod(range(j,j+n))**2,
                 factorial(n)*prod(range(2*j,2*j+n))) for n in range(2*h-j+1))


def fusion(h):
    basis=primitive_basis(h)
    dim=len(basis)
    pairs=[(a,y) for a in basis for y in basis]
    images=[[act(a,y,i) for a,y in pairs] for i in range(2*h+1)]
    channels=[[project(j,x) for x in images[j]] for j in range(2*h+1)]
    coeffs=[weight_coeff(h,j) for j in range(2*h+1)]
    assert all((-1)**j*c>0 for j,c in enumerate(coeffs))
    # The D covariance and reconstruction of each fusion vector.
    for i in range(1,2*h+1):
        for ay in range(dim*dim):
            assert shift(images[i][ay],-1)==scaled(images[i-1][ay],i-1),(h,i,ay,'D')
            expected={}
            for j in range(1,i+1):
                n=i-j
                v=channels[j][ay]
                for _ in range(n): v=shift(v,1)
                factor=F(prod(range(j,j+n)),factorial(n)*prod(range(2*j,2*j+n)))
                add_scaled(expected,v,factor)
            assert images[i][ay]==expected,(h,i,ay,'ladder',images[i][ay],expected)
    # Primitive skew symmetry; compare exact full Q0 with both Gram sums.
    for j in range(2*h+1):
        for a in range(dim):
            for y in range(dim):
                assert channels[j][y*dim+a]==scaled(channels[j][a*dim+y],(-1)**j),(h,j,a,y,'skew')
    full=[]
    swapped=[]
    for ay,(a,y) in enumerate(pairs):
        row=[]
        for bu,(b,u) in enumerate(pairs):
            gram=sum((-1)**i*comb(2*h,i)*dot(images[i][ay],images[i][bu]) for i in range(2*h+1))
            primary=sum(coeffs[j]*dot(channels[j][ay],channels[j][bu]) for j in range(2*h+1))
            assert gram==primary,(h,ay,bu,'gram')
            # Direct source operator products, independently of the adjoint identity.
            direct=0
            for i in range(2*h+1):
                result=act(signs(b),images[i][ay],h)
                direct+=(-1)**i*comb(2*h,i)*dot(u,result)
            assert direct==gram,(h,ay,bu,'adjoint',direct,gram)
            row.append(direct)
        full.append(row)
    for a in range(dim):
        for y in range(dim): swapped.append(full[y*dim+a])
    for ay in range(dim*dim):
        for bu in range(dim*dim):
            psd=sum((-1)**j*coeffs[j]*dot(channels[j][ay],channels[j][bu]) for j in range(2*h+1))
            assert swapped[ay][bu]==psd
    # At h=4 there is a nonzero odd primitive channel. Omitting the row
    # swap therefore leaves a negative direction, not a positive Gram form.
    if h == 4:
        left, right = 1, dim
        antisymmetric = full[left][left] + full[right][right] - full[left][right] - full[right][left]
        assert antisymmetric < 0
    rank=rank_q(full)
    stacked=[[channels[j][ay].get(mu,0) for j in range(2*h+1) for mu in partitions(j)] for ay in range(dim*dim)]
    stacked_rank=rank_q(stacked)
    assert rank==stacked_rank
    individual=[rank_q([[channel.get(mu,0) for mu in partitions(j)] for channel in channels[j]]) for j in range(2*h+1)]
    return {'h':h,'primitive_dimension':dim,'central_primary_rank':rank,'stacked_fusion_rank':stacked_rank,'individual_primary_fusion_ranks':individual,'coefficient_magnitudes':[str((-1)**j*c) for j,c in enumerate(coeffs)]}

def verify():
    records = [fusion(h) for h in range(1, 7)]
    return {
        'scope': 'exact finite primitive central-piece identities and ranks; no asymptotic image estimate',
        'weights': records,
        'direct_central_entries': sum(item['primitive_dimension']**4 for item in records),
        'covariance_and_ladder_vectors': sum(2*item['h']*item['primitive_dimension']**2 for item in records),
        'skew_symmetry_vectors': sum((2*item['h']+1)*item['primitive_dimension']**2 for item in records),
        'negative_control': 'omitting the row swap leaves a negative direction at h=4',
        'new_factoring_exponent_claimed': False,
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
