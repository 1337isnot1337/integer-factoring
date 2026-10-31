#!/usr/bin/env python3
"""Exact three-CNOT crossing modules and all-power coefficient obstruction."""
from __future__ import annotations
import argparse,json,math
from functools import reduce
from pathlib import Path
from check_nilpotent_flag_crossings import (
    compose_columns,column_cut,inverse,check_orders,reverse_at_degree,det_poly,
)


def crossing(q):
    if q<2 or q&(q-1):raise ValueError('q must be a power of two')
    n=q*q;identity=[1<<i for i in range(n)]
    a=[1<<((u^v)*q+v) for u in range(q) for v in range(q)]
    b=[1<<(u*q+(u^v)) for u in range(q) for v in range(q)]
    ab=compose_columns(a,b);ba=compose_columns(b,a);flip=compose_columns(a,ba)
    assert compose_columns(a,a)==identity==compose_columns(b,b)
    assert flip==[1<<((i%q)*q+i//q) for i in range(n)]
    coefficients=[identity,[x^y for x,y in zip(identity,b)],
                  [x^y for x,y in zip(ab,ba)],
                  [w^x^y^z for w,x,y,z in zip(b,ab,ba,flip)]]
    polynomial=[[sum(((coefficients[d][j]>>i)&1)<<d for d in range(4))
                 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            assert polynomial[i][j].bit_count()%2==int(j==(i%q)*q+i//q)
    return polynomial


def change_basis(polynomial,basis):
    n=len(basis)
    matrix=[[(v>>i)&1 for v in basis] for i in range(n)]
    inverse_matrix=inverse(matrix)
    return [[reduce(int.__xor__,
                    (polynomial[a][b] for a in range(n) for b in range(n)
                     if inverse_matrix[i][a] and matrix[b][j]),0)
             for j in range(n)] for i in range(n)]


def finite_lower_bound(t):
    products=[math.comb(t,k)*2**(t-k)*math.comb(2*t,3*k-t)
              for k in range(t+1) if 0<=3*k-t<=2*t]
    return 2**t+sum(math.isqrt(product) for product in products)


def verify():
    K=crossing(2)
    assert K==[[1,0,0,0],[0,9,4,12],[0,4,3,6],[0,12,6,11]]
    row_basis=[1,14,12,10]
    row_normal=change_basis(K,row_basis)
    assert row_normal==[[1,0,0,0],[0,1,0,0],[0,0,5,2],[0,0,8,5]]
    column_basis=[6,9,1,2]
    column_normal=change_basis(column_cut(K,2),column_basis)
    assert column_normal==[[4,12,0,4],[0,8,3,6],[0,0,2,6],[0,0,0,4]]
    assert det_poly(K)==1 and det_poly(column_cut(K,2))==1<<8
    row_orders=check_orders(reverse_at_degree(K,3),[0,3,3,6],7,[1,2,3,4,6,7])
    column_orders=check_orders(column_cut(K,2),[0,1,2,5],6,[1,2,3,5,6])
    finite=[]
    for t in range(1,101):
        bound=finite_lower_bound(t)
        assert bound**70>2**(99*t)
        finite.append({'tensor_exponent':t,'integer_cost_lower_bound':bound})
    # Constants used in the all-large-t proof, with exact integer comparisons.
    assert 99**2>2*70**2
    assert 9**105>2**332
    assert 6**3*16<16**3  # 6*2^(4/3)<16.
    assert 4**100>3**100*16*101*201
    assert 2*100**2-9*100-14>0
    # Literal larger-label endpoint controls of the same all-size formula.
    for q in [4,8]:crossing(q)
    return {'local_polynomial':K,'row_basis_columns':row_basis,
            'row_normal_form':row_normal,'column_basis_columns':column_basis,
            'column_triangular_form':column_normal,
            'row_infinity_orders':row_orders,'column_zero_orders':column_orders,
            'finite_tensor_lower_bounds':finite,
            'large_t_threshold':100,
            'claim_scope':'all ordinary univariate tensor-power coefficient decompositions; broader recombinations not covered'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=verify();rendered=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered,end='')
