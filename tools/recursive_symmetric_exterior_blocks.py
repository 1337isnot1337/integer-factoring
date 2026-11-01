#!/usr/bin/env python3
"""Polynomial-size simultaneous blocks for symmetric exterior multiplication.

Every block change is constant over F2 and works over any commutative
characteristic-two coefficient ring. Multiplicities are ordinary integers.
"""
from collections import Counter
from functools import lru_cache
from math import comb

@lru_cache(None)
def block_types(generators):
    if type(generators) is not int or generators<0:
        raise ValueError("Use a nonnegative integer generator count")
    if generators==0:return {():1}
    result=Counter()
    if generators&1:
        kind=1;pairs=generators//2
    else:
        kind=2;pairs=generators//2-1
    for active in range(pairs+1):
        factor=comb(pairs,active)*(1<<(pairs-active))
        for tail,multiplicity in block_types(active).items():
            result[(kind,)+tail]+=factor*multiplicity
    assert sum(m*(1<<sum(t)) for t,m in result.items())==1<<generators
    assert len(result)==(generators+1)//2
    return dict(result)

def coefficient_indices(kinds):
    if any(kind not in (1,2) for kind in kinds):
        raise ValueError("A block type uses only the symbols 1 and 2")
    result=[0]
    for level,kind in enumerate(kinds):
        old=result[:];weight=1<<level
        if kind==1:result=old+[None if x is None else x+weight for x in old]
        else:result=old+[None if x is None else x+weight for x in old]+[None]*len(old)+[None if x is None else x+2*weight for x in old]
    return result

def recursive_blocks(coefficients):
    if not coefficients:
        raise ValueError("Supply coefficients a_0 through a_r")
    for kinds,multiplicity in block_types(len(coefficients)-1).items():
        indices=coefficient_indices(kinds)
        assert all(i is None or i<len(coefficients) for i in indices)
        yield multiplicity,[0 if i is None else coefficients[i] for i in indices],kinds
