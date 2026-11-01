#!/usr/bin/env python3
"""Bare necessary box budgets applied to word-controlled adjugate crossings."""
from fractions import Fraction as F
from math import isqrt

from check_bare_seven_to_twelve_letters import constraints
from higher_jet_word_caps import word_cost_profile


def inequalities(degree,first,second):
    r,s = tuple(map(F,first)),tuple(map(F,second))
    q = isqrt(len(r))
    if (q < 2 or len(r) != q*q or len(s) != q*q or tuple(sorted(r)) != r
            or tuple(sorted(s)) != s or min(r+s) < 0 or degree < 0):
        raise ValueError('Use two complete nonnegative sorted square profiles and a degree bound.')
    result = []
    for side,(r,s) in enumerate(((r,s),(s,r))):
        for m in range(2,q+1):
            n = m*m
            source = r[:n]
            complement = tuple(sum(source)-value for value in source[::-1])
            word = word_cost_profile(m,s[:n])
            for name,a,b,N in constraints(m):
                lhs = sum(x*y for x,y in zip(a,complement))+sum(x*y for x,y in zip(b,word))
                rhs = F(N*(n-1))*degree
                result.append(dict(source_orientation=side,physical_size=m,box_orientation=name[0],
                                   alpha=name[1],beta=name[2],lhs=lhs,rhs=rhs,slack=rhs-lhs))
    return result


def first_jet_linear_caps(q):
    """A weaker linear consequence using H_i(s)>=s_2 for i>=2."""
    result = []
    for m in range(2,q+1):
        n = m*m
        for name,a,b,N in constraints(m):
            first = [N-value for value in a[::-1]]+[0]*(q*q-n)
            second = [0]*(q*q)
            second[1] = N-b[0]
            assert min(first+second) >= 0
            result.append(dict(physical_size=m,box=name,first=first,second=second,
                               degree_multiplier=N*(n-1)))
    return result
