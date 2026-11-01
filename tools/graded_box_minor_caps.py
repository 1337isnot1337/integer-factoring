#!/usr/bin/env python3
"""Necessary original and word-adjugate box budgets with physical energies."""
from fractions import Fraction as F
from itertools import accumulate

from higher_jet_word_caps import word_cost_profile


def data(deficit,first,second,energy_x,energy_y):
    x,y = tuple(sorted(map(F,energy_x))),tuple(sorted(map(F,energy_y)))
    r,s = tuple(map(F,first)),tuple(map(F,second))
    q,d = len(x),F(deficit)
    if (q < 2 or len(y) != q or x[0] or y[0] or d < 0 or len(r) != q*q or len(s) != q*q
            or r[0] or s[0] or tuple(sorted(r)) != r or tuple(sorted(s)) != s):
        raise ValueError('Use normalized equal paired alphabets and complete sorted profiles.')
    return q,d,r,s,x,y


def prefix(values):
    return (F(0),)+tuple(accumulate(values))


def box_values(values,sums,alpha,beta):
    """Weighted box prefix in O(1), without a dense coefficient vector."""
    size,rank = alpha*beta,alpha*alpha
    full,rest = divmod(size,rank)
    return rank*sums[full]+(rest*values[full] if rest else 0)


def original(deficit,first,second,energy_x,energy_y):
    q,d,r,s,x,y = data(deficit,first,second,energy_x,energy_y)
    pr,ps,px,py = map(prefix,(r,s,x,y))
    records = []
    for alpha in range(1,q+1):
        for beta in range(alpha,q+1):
            N = alpha*beta
            physical = min(beta*px[alpha]+alpha*py[beta],alpha*px[beta]+beta*py[alpha])
            rhs = N*d+physical
            for side,lhs in (('C',box_values(r,pr,alpha,beta)+ps[N]),
                             ('R',pr[N]+box_values(s,ps,alpha,beta))):
                records.append(dict(box_orientation=side,alpha=alpha,beta=beta,lhs=lhs,rhs=rhs,slack=rhs-lhs))
    return records


def adjugate(deficit,first,second,energy_x,energy_y):
    q,d,r,s,x,y = data(deficit,first,second,energy_x,energy_y)
    px,py = prefix(x),prefix(y)
    records = []
    for side,(r,s) in enumerate(((r,s),(s,r))):
        for m in range(2,q+1):
            n = m*m
            C = (n-1)*d+m*(px[m]+py[m])
            source = r[:n]
            complement = tuple(sum(source)-v for v in source[::-1])
            word = word_cost_profile(m,s[:n])
            pc,pw = prefix(complement),prefix(word)
            for alpha in range(1,m+1):
                for beta in range(alpha,m+1):
                    N = alpha*beta
                    physical = max(beta*(px[m]-px[m-alpha])+alpha*(py[m]-py[m-beta]),
                                   alpha*(px[m]-px[m-beta])+beta*(py[m]-py[m-alpha]))
                    rhs = N*C-physical
                    for orientation,lhs in (('C',box_values(complement,pc,alpha,beta)+pw[N]),
                                            ('R',pc[N]+box_values(word,pw,alpha,beta))):
                        records.append(dict(source_orientation=side,physical_size=m,box_orientation=orientation,
                                            alpha=alpha,beta=beta,lhs=lhs,rhs=rhs,slack=rhs-lhs))
    return records
