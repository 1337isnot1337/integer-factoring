#!/usr/bin/env python3
"""All-minor exact controls for the supplied-rank Laurent F4 interface."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations,permutations
import json
from pathlib import Path
import random

from native_matrix_smith_f4 import MatrixOracle
from weighted_smith_oracle_f4 import polynomial_multiply,scalar_multiply


def planes(value):
    lo=hi=0
    for k in range((value.bit_length()+1)//2):
        c=(value>>(2*k))&3;lo|=(c&1)<<k;hi|=(c>>1)<<k
    return lo,hi


def multiply(a,b):
    out=[[0]*len(b[0]) for _ in a]
    for i,row in enumerate(a):
        for k,v in enumerate(row):
            if v:
                for j,w in enumerate(b[k]):out[i][j]^=polynomial_multiply(v,w)
    return out


def determinant(matrix):
    n=len(matrix);value=0
    for perm in permutations(range(n)):
        term=1
        for i,j in enumerate(perm):term=polynomial_multiply(term,matrix[i][j])
        value^=term
    return value


def check():
    rng=random.Random(5186);records=[];maximum=0.;minors=0
    for n in (2,3,4):
        for rank in range(1,n+1):
            left=[[int(i==j) for j in range(n)] for i in range(n)]
            right=[[int(i==j) for j in range(n)] for i in range(n)]
            for i in range(n):
                for j in range(i+1,n):
                    left[i][j]=rng.randrange(1,64);right[j][i]=rng.randrange(1,64)
            diagonal=[[((1+rng.randrange(3))<<(2*rng.randrange(4))) if i==j and i<rank else 0
                       for j in range(n)] for i in range(n)]
            matrix=multiply(multiply(left,diagonal),right)
            row_shift=[rng.randrange(-3,4) for _ in range(n)]
            column_shift=[rng.randrange(-3,4) for _ in range(n)]
            shifts=[[row_shift[i]+column_shift[j] for j in range(n)] for i in range(n)]
            oracle=MatrixOracle([[planes(v) for v in row] for row in matrix],rank,shifts)
            rw=[Fraction(rng.randrange(-15,16),29) for _ in range(n)]
            cw=[Fraction(rng.randrange(-15,16),31) for _ in range(n)]
            scale=Fraction(2,7)
            for infinity in (False,True):
                sign=-1 if infinity else 1
                minima=[Fraction(0)]
                for k in range(1,rank+1):
                    choices=[]
                    for rows in combinations(range(n),k):
                        for cols in combinations(range(n),k):
                            value=determinant([[matrix[i][j] for j in cols] for i in rows]);minors+=1
                            if not value:continue
                            valuation=-((value.bit_length()-1)//2) if infinity else ((value & -value).bit_length()-1)//2
                            valuation+=sign*(sum(row_shift[i] for i in rows)+sum(column_shift[j] for j in cols))
                            choices.append(scale*valuation+sum(rw[i] for i in rows)+sum(cw[j] for j in cols))
                    minima.append(min(choices))
                expected=[minima[k+1]-minima[k] for k in range(rank)]
                actual=oracle.profile(rw,cw,scale,infinity)
                error=max(abs(float(a)-b) for a,b in zip(expected,actual['orders']));maximum=max(maximum,error)
                assert error<1e-10,(n,rank,infinity,expected,actual)
                for k in range(rank):
                    rows=actual['rows'][:k+1];cols=actual['cols'][:k+1]
                    value=determinant([[matrix[i][j] for j in cols] for i in rows])
                    valuation=-((value.bit_length()-1)//2) if infinity else ((value & -value).bit_length()-1)//2
                    valuation+=sign*(sum(row_shift[i] for i in rows)+sum(column_shift[j] for j in cols))
                    assert valuation==actual['base_prefix'][k]
                records.append({'matrix_order':n,'supplied_generic_rank':rank,'infinity':infinity,
                                'exact_orders':[str(v) for v in expected]})
    bad=MatrixOracle([[(1,0),(0,0)],[(0,0),(0,0)]],2)
    try:bad.profile()
    except ArithmeticError:pass
    else:raise AssertionError('An unresolved promised rank was silently accepted.')
    return {'status':'pass','profile_pairs':len(records),'all_minors_examined':minors,
            'maximum_floating_error':maximum,'missing_promised_rank_rejected':True,'records':records,
            'scope':'Generic rank is supplied and independently known from triangular factors; truncated zeros do not certify rank.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=check();text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text)


if __name__=='__main__':main()
