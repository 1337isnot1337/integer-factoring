#!/usr/bin/env python3
"""Weighted Hilbert counts for kernels of disjoint paired derivatives."""
from collections import Counter
from math import comb
from coherent_hot_compounds_f4 import basis


def times(a,b):
    result=Counter()
    for (h,k),v in a.items():
        for (g,l),w in b.items():result[h+g,k+l]+=v*w
    return result


def hilbert(pairs,coefficients,modes):
    used=[i for pair in pairs for i in pair]
    if len(set(used))!=len(used) or any(not 0<=i<modes for i in used):raise ValueError('Pairs must be disjoint and in range.')
    m=len(pairs)
    if not coefficients or any(len(row)!=m or any(c not in range(4) for c in row) for row in coefficients):raise ValueError('Invalid F4 coefficient matrix.')
    heights={2*a+2*b+2 for a,b in pairs}
    if len(heights)!=1:raise ValueError('A common pair height is required for the weighted formula.')
    H=heights.pop();answer=Counter();ranks=[]
    for mask in range(1<<m):
        active=[i for i in range(m) if mask>>i&1]
        rank=len(basis([{j:row[i] for j,row in enumerate(coefficients) if row[i]} for i in active]));ranks.append(rank)
        free=len(active)-rank;part=Counter({(H*j,2*j):comb(free,j) for j in range(free+1)})
        for i,(a,b) in enumerate(pairs):
            if not mask>>i&1:part=times(part,Counter({(2*a+1,1):1,(2*b+1,1):1}))
        answer.update(part)
    for a in range(modes):
        if a not in used:answer=times(answer,Counter({(0,0):1,(2*a+1,1):1}))
    return answer,ranks
