#!/usr/bin/env python3
"""Determinant-minor inverse systems and three prescribed F8 retractions.

Finite sparse map construction supports k<=4. The height-character algorithm
does not enumerate minors and works for arbitrary positive k.
"""
from itertools import combinations,permutations
from math import comb

from coherent_f8_maps import mul,validate_maps

INV=(0,1,5,6,7,2,3,4)


def minor_vector(k,rows,columns,weights=None):
    result={}
    for order in permutations(columns):
        word=0;coefficient=1
        for i,j in zip(rows,order):
            word|=1 << (k*i+j)
            if weights is not None:coefficient=mul(coefficient,weights[i][j])
        if coefficient:
            result[word]=result.get(word,0)^coefficient
            if not result[word]:del result[word]
    return result


def scalar_minor(matrix,rows,columns):
    result=0
    for order in permutations(columns):
        value=1
        for i,j in zip(rows,order):value=mul(value,matrix[i][j])
        result^=value
    return result


def cauchy_matrix(k):
    if not 1<=k<=4:raise ValueError('The prescribed F8 Cauchy points support k=1..4.')
    return [[INV[i ^ j] for j in range(k,2*k)] for i in range(k)]


def code(k,kind):
    if not isinstance(k,int) or not 1<=k<=4:raise ValueError('Sparse maps are implemented for k=1..4.')
    if kind not in ('diagonal','antidiagonal','cauchy'):raise ValueError('Unknown prescribed retraction.')
    matrix=cauchy_matrix(k);labels=[]
    for r in range(k+1):
        for rows in combinations(range(k),r):
            for cols in combinations(range(k),r):
                h=2*k*sum(rows)+2*sum(cols)+r
                labels.append(dict(rows=list(rows),columns=list(cols),height=h,charge=r))
    labels.sort(key=lambda x:(x['height'],x['charge'],x['rows'],x['columns']))
    embedding=[];retraction=[]
    for label in labels:
        rows,cols=label['rows'],label['columns']
        embedding.append(minor_vector(k,rows,cols))
        if kind=='cauchy':
            norm=scalar_minor(matrix,rows,cols)
            assert norm
            left={w:mul(INV[norm],v) for w,v in minor_vector(k,rows,cols,matrix).items()}
            label['pairing_norm']=norm
        else:
            chosen=cols if kind=='diagonal' else list(reversed(cols))
            left={sum(1 << (k*i+j) for i,j in zip(rows,chosen)):1}
        retraction.append(left)
    hs,charges=validate_maps(embedding,retraction,k*k)
    assert hs==[x['height'] for x in labels] and charges==[x['charge'] for x in labels]
    assert len(labels)==comb(2*k,k)
    return dict(k=k,kind=kind,modes=k*k,q=len(labels),top_height=k**3,
        labels=labels,embedding=embedding,retraction=retraction,
        cauchy_matrix=matrix if kind=='cauchy' else None)


def character(k):
    """Charge/height enumerator by elementary-symmetric dynamic programming."""
    if not isinstance(k,int) or k<1:raise ValueError('Use a positive integer k.')
    def elementary(heights):
        result=[{0:1}]+[{} for _ in range(k)]
        for index,h in enumerate(heights):
            for r in range(index+1,0,-1):
                for e,c in result[r-1].items():
                    result[r][e+h]=result[r].get(e+h,0)+c
        return result
    row=elementary([2*k*i for i in range(k)])
    col=elementary([2*j for j in range(k)])
    result={}
    for r in range(k+1):
        values={}
        for a,c in row[r].items():
            for b,d in col[r].items():values[a+b+r]=values.get(a+b+r,0)+c*d
        assert sum(values.values())==comb(k,r)**2
        result[r]=values
    assert sum(sum(v.values()) for v in result.values())==comb(2*k,k)
    return result


def serializable(value):
    if isinstance(value,dict):return {str(k):serializable(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serializable(v) for v in value]
    return value
