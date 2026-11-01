#!/usr/bin/env python3
"""Bound all algebraic local determinant orders using rational-place profiles."""
from numbers import Integral

def all_place_bound(profiles,degree,field_size):
    if type(degree)is not int or degree<0:raise ValueError('Nonnegative polynomial degree required.')
    if type(field_size)is not int or field_size<2 or field_size&(field_size-1):raise ValueError('This interface uses binary finite fields.')
    if len(profiles)!=field_size+1:raise ValueError('Include every rational point and infinity.')
    n=len(profiles[0])
    if not n or any(len(p)!=n or any(type(x)is not int or x<0 for x in p) for p in profiles):raise ValueError('Complete nonnegative local profiles required.')
    sums=list(map(sum,profiles));remaining=n*degree-sum(sums)
    if remaining<0:raise ValueError('Local orders exceed the homogeneous determinant degree.')
    return {'matrix_order':n,'degree':degree,'homogeneous_determinant_degree':n*degree,'rational_order_sums':sums,'unaccounted_degree':remaining,'nonrational_point_order_bound':remaining//2,'all_algebraic_point_order_bound':max(max(sums),remaining//2)}


def _times_x(row):
    a,b,c=row;return c,a^c,b

def scalar_row(row,c):
    out=(0,0,0);term=row
    while c:
        if c&1:out=tuple(a^b for a,b in zip(out,term))
        c>>=1;term=_times_x(term)
    return out

def rank_f8_rows(rows):
    """Independent exact Gaussian rank on three bit planes per row."""
    basis={};inverse=(0,1,5,6,7,2,3,4)
    for row in rows:
        while row[0]|row[1]|row[2]:
            bits=row[0]|row[1]|row[2];j=(bits&-bits).bit_length()-1
            c=sum(((v>>j)&1)<<a for a,v in enumerate(row))
            if j in basis:row=tuple(a^b for a,b in zip(row,basis[j][c]))
            else:
                row=scalar_row(row,inverse[c]);basis[j]=[scalar_row(row,c)for c in range(8)];break
    return len(basis)


def toeplitz_order_sum(coefficients,known_orders):
    """Replay the determinant order using a finite block-Toeplitz kernel."""
    if not known_orders or not len(coefficients) or any(type(x)is not int or x<0 for x in known_orders):raise ValueError('Complete nonnegative orders and coefficients required.')
    n=len(known_orders);K=max(known_orders)+1
    if any(len(A)!=n or any(len(row)!=n for row in A)for A in coefficients):raise ValueError('Square coefficient matrices required.')
    layers=[]
    for A in coefficients[:K]:
        rows=[]
        for row in A:
            if any(not isinstance(c,Integral) or not 0<=int(c)<8 for c in row):raise ValueError('Invalid F8 coefficient.')
            rows.append(tuple(sum(((int(c)>>a)&1)<<j for j,c in enumerate(row))for a in range(3)))
        layers.append(rows)
    def rows():
        for s in range(K):
            for i in range(n):
                row=[0,0,0]
                for d in range(min(s+1,len(layers))):
                    for a in range(3):row[a]^=layers[d][i][a]<<((s-d)*n)
                yield tuple(row)
    rank=rank_f8_rows(rows());defect=n*K-rank
    if defect!=sum(known_orders):raise ArithmeticError('Independent Toeplitz determinant-order replay failed.')
    return {'jet_length':K,'toeplitz_order':n*K,'rank':rank,'kernel_dimension':defect}
