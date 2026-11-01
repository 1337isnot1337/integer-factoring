#!/usr/bin/env python3
"""Exact orthogonal retracts of homogeneous F4 differential kernels."""
from coherent_hot_compounds_f4 import mul,axpy,dot,basis,validate_maps


def inverse(matrix):
    n=len(matrix)
    if any(len(row)!=n for row in matrix):raise ValueError('Expected a square matrix.')
    A=[row[:]+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return None
        A[pivot],A[j]=A[j],A[pivot];c=mul(A[j][j],A[j][j]);A[j]=[mul(c,x) for x in A[j]]
        for i in range(n):
            if i!=j:
                c=A[i][j];A[i]=[x^mul(c,y) for x,y in zip(A[i],A[j])]
    return [row[n:] for row in A]


def rank(matrix):
    return len(basis([{j:x for j,x in enumerate(row) if x} for row in matrix]))


def subset_gram_certificate(coefficients):
    if not coefficients or not coefficients[0]:raise ValueError('Expected nonempty coefficients.')
    m=len(coefficients[0])
    if any(len(row)!=m or any(x not in range(4) for x in row) for row in coefficients):raise ValueError('Expected a rectangular F4 matrix.')
    records=[]
    for mask in range(1<<m):
        A=[[x for j,x in enumerate(row) if mask>>j&1] for row in coefficients]
        G=[[0 for _ in A] for _ in A]
        for i,row in enumerate(A):
            for j,other in enumerate(A):
                for a,b in zip(row,other):G[i][j]^=mul(a,b)
        r=rank(A);g=rank(G);records.append({'active_mask':mask,'rank':r,'gram_rank':g})
    return {'all_subsets_nondegenerate':all(r['rank']==r['gram_rank'] for r in records),'subsets':records}


def orthogonal_retraction(code):
    """Return L=(P^T P)^(-1)P^T, or None for a degenerate code."""
    p=code['embedding'];l=[{} for _ in p];cells={}
    for i,cell in enumerate(zip(code['heights'],code['charges'])):cells.setdefault(cell,[]).append(i)
    for indices in cells.values():
        G=[[dot(p[i],p[j]) for j in indices] for i in indices];V=inverse(G)
        if V is None:return None
        for i,row in zip(indices,V):
            for c,j in zip(row,indices):l[i]=axpy(l[i],c,p[j])
    validate_maps(p,l,code['modes'])
    return l
