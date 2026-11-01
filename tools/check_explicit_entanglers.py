#!/usr/bin/env python3
"""Replay exact finite-field universal-entangler certificates and a density control.

A certificate has no product-vector solution even over an algebraic closure.
The F4 control still has a normalized product-density scalarization.
No improving crossing seed or factoring bound is asserted.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

MUL=((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))
INV=(0,1,3,2)


def monomials(n,degree):
    if n==1:
        return [(degree,)]
    return [(i,)+tail for i in range(degree+1)
            for tail in monomials(n-1,degree-i)]


def xor_sum(values):
    out=0
    for value in values:
        out^=value
    return out


def matrix_product(A,B):
    return [[xor_sum(MUL[A[i][k]][B[k][j]] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def identity(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]


def inverse(A):
    n=len(A)
    M=[row[:]+I for row,I in zip(A,identity(n))]
    for column in range(n):
        pivot=next((j for j in range(column,n) if M[j][column]),None)
        if pivot is None:
            return None
        M[column],M[pivot]=M[pivot],M[column]
        M[column]=[MUL[INV[M[column][column]]][x] for x in M[column]]
        for j in range(n):
            if j!=column and M[j][column]:
                coefficient=M[j][column]
                M[j]=[x^MUL[coefficient][y] for x,y in zip(M[j],M[column])]
    return [row[n:] for row in M]


def tensor(A,B):
    return [[MUL[A[i][j]][B[k][l]]
             for j in range(len(A[0])) for l in range(len(B[0]))]
            for i in range(len(A)) for k in range(len(B))]


def multiply_row(row,scalar):
    a,b=row
    if scalar==0:return 0,0
    if scalar==1:return a,b
    if scalar==2:return b,a^b
    return a^b,a


def echelon_rank(rows):
    pivots={}
    for row in rows:
        while row[0] or row[1]:
            pivot=(row[0]|row[1]).bit_length()-1
            coefficient=((row[0]>>pivot)&1)+2*((row[1]>>pivot)&1)
            if pivot in pivots:
                a,b=multiply_row(pivots[pivot],coefficient)
                row=row[0]^a,row[1]^b
            else:
                pivots[pivot]=multiply_row(row,INV[coefficient])
                break
    return len(pivots)


def minor_equations(W,a=4,b=4):
    result=[]
    for i,j in itertools.combinations(range(a),2):
        for k,l in itertools.combinations(range(b),2):
            coefficients={}
            for P,Q in ((W[b*i+k],W[b*j+l]),(W[b*i+l],W[b*j+k])):
                for u,v in itertools.product(range(a*b),repeat=2):
                    scalar=MUL[P[u]][Q[v]]
                    if not scalar:continue
                    first=tuple(int(t==u//b)+int(t==v//b) for t in range(a))
                    second=tuple(int(t==u%b)+int(t==v%b) for t in range(b))
                    key=first,second
                    coefficients[key]=coefficients.get(key,0)^scalar
            result.append({key:scalar for key,scalar in coefficients.items() if scalar})
    return result


def replay(data,W):
    assert len(W)==16 and all(len(row)==16 for row in W)
    field=data['field']
    assert field in (2,4) and all(0<=x<field for row in W for x in row)
    W_inverse=inverse(W)
    assert W_inverse is not None
    assert matrix_product(W,W_inverse)==identity(16)
    certificate=data['certificate']
    assert certificate['degree']==3
    bases=monomials(4,3)
    low=monomials(4,1)
    index={m:i for i,m in enumerate(bases)}
    equations=minor_equations(W)
    assert len(equations)==36 and len(bases)==20
    selected=certificate['selected_rows']
    assert len(selected)==len(set(map(tuple,selected)))==400
    rows=[]
    for ei,ai,bi in selected:
        assert 0<=ei<36 and 0<=ai<4 and 0<=bi<4
        bits=[0,0]
        for (first,second),scalar in equations[ei].items():
            u=tuple(x+y for x,y in zip(first,low[ai]))
            v=tuple(x+y for x,y in zip(second,low[bi]))
            column=20*index[u]+index[v]
            if scalar&1:bits[0]^=1<<column
            if scalar&2:bits[1]^=1<<column
        rows.append(tuple(bits))
    assert echelon_rank(rows)==400
    return W_inverse


def check_density_control(data,W,W_inverse):
    first,second=data['density_factors_diagonal']
    A=[[first[i] if i==j else 0 for j in range(4)] for i in range(4)]
    B=[[second[i] if i==j else 0 for j in range(4)] for i in range(4)]
    D=tensor(A,B)
    assert xor_sum(D[i][i] for i in range(16))==1
    assert matrix_product(matrix_product(W_inverse,D),W)==D
    # Check the factorizing functional on all matrix-unit pairs.
    units=[]
    for i,j in itertools.product(range(4),repeat=2):
        U=[[0]*4 for _ in range(4)];U[i][j]=1;units.append(U)
    checks=0
    for P,Q in itertools.product(units,repeat=2):
        X=tensor(P,Q)
        expected=MUL[xor_sum(matrix_product(A,P)[i][i] for i in range(4))][xor_sum(matrix_product(B,Q)[i][i] for i in range(4))]
        for T in (X,matrix_product(matrix_product(W,X),W_inverse)):
            DT=matrix_product(D,T)
            assert xor_sum(DT[i][i] for i in range(16))==expected
            checks+=1
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=Path(__file__).resolve().parents[1]/'evidence'/'universal-entangler16.json')
    parser.add_argument('--control',type=Path,default=Path(__file__).resolve().parents[1]/'evidence'/'entangler-density-control.json')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();started=time.monotonic()
    binary=json.loads(args.binary.read_text())
    W2=[[(row>>i)&1 for i in range(16)] for row in binary['matrix_rows']]
    replay(binary,W2)
    control=json.loads(args.control.read_text())
    W4=control['matrix'];W4_inverse=replay(control,W4)
    functional_checks=check_density_control(control,W4,W4_inverse)
    output={'status':'passed','field2_macaulay_rank':400,'field4_macaulay_rank':400,
            'bilinear_scalarization_checks':functional_checks,
            'field2_certificate_sha256':hashlib.sha256(args.binary.read_bytes()).hexdigest(),
            'field4_control_sha256':hashlib.sha256(args.control.read_bytes()).hexdigest(),
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'seconds':time.monotonic()-started,
            'scope':'Explicit product-ket exclusions over algebraic closures, and an explicit product-density scalarization surviving such an exclusion. No crossing improvement.'}
    print(json.dumps(output,indent=2,sort_keys=True))
    if args.output:args.output.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':
    main()
