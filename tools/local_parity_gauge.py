#!/usr/bin/env python3
"""Exact constant local parity gauges for polynomial binary-leg F4 tensors.

Coefficients 0,1,2,3 denote 0,1,alpha,1+alpha in F4. Returned gauges use
F16 = F2[t]/(t^4+t+1), with F4 embedded by alpha -> t^2+t (the integer 6).
The linear test is necessary and sufficient over coefficient-field extensions.
It does not permit parameter-dependent rational changes of basis.
"""
from __future__ import annotations

F4_MUL = ((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))
F4_INV = (0,1,3,2)
F4_EMBED = (0,1,6,7)


def _f16_multiply(a,b):
    result=0
    while b:
        if b&1:result^=a
        b>>=1;a<<=1
        if a&16:a^=19
    return result


F16_MUL=tuple(tuple(_f16_multiply(a,b) for b in range(16)) for a in range(16))
F16_INV=(0,)+tuple(next(b for b in range(1,16) if F16_MUL[a][b]==1) for a in range(1,16))


def validate(tensor,pairs=None):
    size=len(tensor)
    if size<2 or size&(size-1):raise ValueError('Tensor size must be a power of two, at least two.')
    n=size.bit_length()-1
    if not any(tensor):raise ValueError('The tensor must be nonzero.')
    for polynomial in tensor:
        if any(c not in (1,2,3) for c in polynomial.values()):raise ValueError('Use nonzero F4 coefficients, omitting zero terms.')
    if pairs is not None:
        if sorted(i for pair in pairs for i in pair)!=list(range(n)):
            raise ValueError('Pairs must partition the physical legs.')
        for mask,polynomial in enumerate(tensor):
            endpoint=0
            for c in polynomial.values():endpoint^=c
            expected=int(all((mask>>i&1)==(mask>>j&1) for i,j in pairs))
            if endpoint!=expected:raise ValueError('The specified product-of-EQ endpoint failed.')
    return n


def equation(tensor,exponent,mask,pairs=None):
    n=len(tensor).bit_length()-1;c=tensor[mask].get(exponent,0)
    if pairs is None:
        row=[c]
        row += [tensor[mask^(1<<i)].get(exponent,0) if not(mask>>i&1) else 0 for i in range(n)]
        row += [tensor[mask^(1<<i)].get(exponent,0) if mask>>i&1 else 0 for i in range(n)]
    else:
        row=[]
        for i,j in pairs:
            ti=tensor[mask^(1<<i)].get(exponent,0);tj=tensor[mask^(1<<j)].get(exponent,0)
            row.append((ti if not(mask>>i&1) else 0)^(tj if mask>>j&1 else 0))
            row.append((ti if mask>>i&1 else 0)^(tj if not(mask>>j&1) else 0))
    row.append(c if mask.bit_count()%2 else 0)
    return row


def _add_scaled(first,second,c):
    result=dict(first)
    for key,value in second.items():
        result[key]=result.get(key,0)^F4_MUL[c][value]
        if not result[key]:del result[key]
    return result


def solve(tensor,pairs=None):
    """Solve the 2n+1-variable test, or its n-variable EQ-endpoint reduction."""
    n=validate(tensor,pairs);variables=2*n+1 if pairs is None else n
    basis={};rows=0
    for exponent in sorted({e for polynomial in tensor for e in polynomial}):
        for mask in range(len(tensor)):
            row=equation(tensor,exponent,mask,pairs);combination={(exponent,mask):1};rows+=1
            for pivot in range(variables):
                c=row[pivot]
                if not c:continue
                if pivot not in basis:
                    inv=F4_INV[c]
                    basis[pivot]=([F4_MUL[inv][x] for x in row],
                                  {key:F4_MUL[inv][x] for key,x in combination.items()})
                    break
                other,old_combination=basis[pivot]
                row=[a^F4_MUL[c][b] for a,b in zip(row,other)]
                combination=_add_scaled(combination,old_combination,c)
            else:
                if row[-1]:
                    inv=F4_INV[row[-1]]
                    certificate=[dict(exponent=list(key[0]),mask=key[1],
                                      coefficient=F4_MUL[inv][value])
                                 for key,value in sorted(combination.items())]
                    return dict(equivalent=False,variables=variables,rank=len(basis),tested_equations=rows,
                                obstruction=certificate)
    solution=[0]*variables
    for pivot,(row,_) in sorted(basis.items(),reverse=True):
        value=row[-1]
        for j in range(pivot+1,variables):value^=F4_MUL[row[j]][solution[j]]
        solution[pivot]=value
    return dict(equivalent=True,variables=variables,rank=len(basis),tested_equations=rows,
                solution_dimension=variables-len(basis),solution=solution)


def verify_obstruction(tensor,certificate,pairs=None):
    n=validate(tensor,pairs);variables=2*n+1 if pairs is None else n;total=[0]*(variables+1)
    for item in certificate:
        exponent=tuple(item['exponent']);mask=item['mask'];c=item['coefficient']
        if not 0<=mask<len(tensor) or c not in (1,2,3):raise ValueError('Invalid certificate term.')
        row=equation(tensor,exponent,mask,pairs)
        total=[a^F4_MUL[c][b] for a,b in zip(total,row)]
    if total!=[0]*variables+[1]:raise AssertionError('The claimed inconsistency certificate is not 0=1.')
    return dict(certificate_terms=len(certificate),verified_zero_coefficients=variables,right_hand_side=1)


def transpose(matrix):return [list(row) for row in zip(*matrix)]


def inverse(matrix):
    a,b=matrix[0];c,d=matrix[1];determinant=F16_MUL[a][d]^F16_MUL[b][c]
    if not determinant:raise ValueError('The local matrix is singular.')
    s=F16_INV[determinant]
    return [[F16_MUL[s][d],F16_MUL[s][b]],[F16_MUL[s][c],F16_MUL[s][a]]]


def _idempotent_gauge(operator,root):
    b=[[operator[i][j]^(root if i==j else 0) for j in range(2)] for i in range(2)]
    image=next([b[0][j],b[1][j]] for j in range(2) if b[0][j] or b[1][j])
    kernel=next([b[0][j]^(1 if j==0 else 0),b[1][j]^(1 if j==1 else 0)]
                for j in range(2) if (b[0][j]^(1 if j==0 else 0)) or (b[1][j]^(1 if j==1 else 0)))
    return inverse([[kernel[0],image[0]],[kernel[1],image[1]]])


def find_gauge(tensor,pairs=None):
    """Return either an exact 0=1 certificate or explicit constant F16 gauges."""
    result=solve(tensor,pairs)
    if not result['equivalent']:return result
    n=len(tensor).bit_length()-1;solution=result['solution'];operators=[None]*n
    if pairs is None:
        for i in range(n):
            a=solution[0] if i==0 else 0
            operators[i]=[[F4_EMBED[a],F4_EMBED[solution[1+i]]],
                          [F4_EMBED[solution[1+n+i]],F4_EMBED[a^1]]]
    else:
        for p,(i,j) in enumerate(pairs):
            b,c=solution[2*p:2*p+2]
            operators[i]=[[0,F4_EMBED[b]],[F4_EMBED[c],1]]
            operators[j]=transpose(operators[i])
    determinants=[F16_MUL[a[0][0]][a[1][1]]^F16_MUL[a[0][1]][a[1][0]] for a in operators]
    total=0
    for d in determinants:total^=d
    if total:raise AssertionError('An annihilating trace-one family must have determinant sum zero.')
    roots=[next(r for r in range(16) if F16_MUL[r][r]^r==d) for d in determinants]
    gauges=[None]*n
    if pairs is None:
        root_sum=0
        for root in roots:root_sum^=root
        if root_sum not in (0,1):raise AssertionError('The root sum must lie in F2.')
        roots[0]^=root_sum
        for i in range(n):gauges[i]=_idempotent_gauge(operators[i],roots[i])
    else:
        for i,j in pairs:
            gauges[i]=_idempotent_gauge(operators[i],roots[i])
            gauges[j]=transpose(inverse(gauges[i]));roots[j]=roots[i]
    result.update(gauges_f16=gauges,scalar_shifts_f16=roots,field='F2[t]/(t^4+t+1)',f4_embedding=list(F4_EMBED))
    return result


def embed_tensor(tensor):
    return tuple({e:F4_EMBED[c] for e,c in polynomial.items()} for polynomial in tensor)


def apply_maps(tensor,maps,multiplication=F16_MUL):
    current=tuple(dict(polynomial) for polynomial in tensor)
    for vertex,matrix in enumerate(maps):
        out=[]
        for mask in range(len(current)):
            polynomial={};row=mask>>vertex&1
            for column in range(2):
                coefficient=matrix[row][column]
                if not coefficient:continue
                source=(mask&~(1<<vertex))|(column<<vertex)
                for exponent,value in current[source].items():
                    polynomial[exponent]=polynomial.get(exponent,0)^multiplication[coefficient][value]
            out.append({e:c for e,c in polynomial.items() if c})
        current=tuple(out)
    return current


def verify_gauge(tensor,record,pairs=None):
    n=validate(tensor,pairs)
    if not record.get('equivalent') or len(record.get('gauges_f16',[]))!=n:
        raise AssertionError('Missing constant local gauges.')
    gauges=record['gauges_f16']
    for matrix in gauges:inverse(matrix)
    transformed=apply_maps(embed_tensor(tensor),gauges)
    if any(polynomial for mask,polynomial in enumerate(transformed) if mask.bit_count()%2):
        raise AssertionError('The supplied gauges leave a nonzero odd coefficient.')
    if pairs is not None:
        for i,j in pairs:
            if gauges[j]!=transpose(inverse(gauges[i])):
                raise AssertionError('The paired gauges do not preserve EQ.')
        for mask,polynomial in enumerate(transformed):
            endpoint=0
            for c in polynomial.values():endpoint^=c
            if endpoint!=int(all((mask>>i&1)==(mask>>j&1) for i,j in pairs)):
                raise AssertionError('The transformed endpoint is not EQ.')
    return dict(local_invertible_maps=n,odd_coefficients_zero=True,eq_endpoint_preserved=pairs is not None,
                transformed_nonzero_coefficients=sum(len(p) for p in transformed))
