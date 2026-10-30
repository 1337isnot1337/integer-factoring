#!/usr/bin/env python3
"""Polynomial, overlapping, swap-commuting conjugations of the pairing pencil.

Each local shear is I+p(z)K, K=|00,11><01,10|+|11,00><10,01|.
It squares to I in characteristic two and commutes with the global swap.
The normalized row cut uses input order (v,u).
"""
import argparse
import json
import math
import os
import time
from pathlib import Path


def rank(rows):
    basis={}
    for row in rows:
        while row:
            pivot=row.bit_length()-1
            if pivot not in basis:
                basis[pivot]=row
                break
            row^=basis[pivot]
    return len(basis)


def trim(rows):
    while len(rows)>1 and not any(rows[-1]):
        rows.pop()


def profile(t, order, polynomial):
    n=1<<t
    exponents=[j for j in range(polynomial.bit_length()) if polynomial>>j&1]
    if not exponents:
        raise ValueError("Use a nonzero shear polynomial")
    degree_p=max(exponents)
    if order=="pairs":
        windows=list(range(0,t-1,2))
    elif order=="layers":
        windows=list(range(0,t-1,2))+list(range(1,t-1,2))
    else:
        windows=list(range(t-1))
    max_degree=t+2*degree_p*len(windows)
    columns=[[[0]*n for _ in range(n)] for _ in range(max_degree+1)]
    row_ranks=[0]*(max_degree+1)
    gates=[]
    for start in windows:
        mask=3<<start
        v=2 if order=="alternating" and start%2 else 1
        first=v<<start
        second=(v^3)<<start
        outside=[x for x in range(n) if x&mask==0]
        selector=sum(1<<x for x in outside)
        gates.append((mask,first,second,outside,selector))
    for g in range(n):
        offset=t-g.bit_count()
        kernels=[[] for _ in range(t+1)]
        s=g
        while True:
            sw=s.bit_count()
            for k in range(sw+1):
                if k&~sw==0:
                    kernels[offset+k].append(s)
            if s==0:
                break
            s=(s-1)&g
        rows=[[sum(1<<(x^s) for s in support) for x in range(n)]
              for support in kernels]
        for mask,first,second,outside,selector in gates:
            if g&mask!=mask:
                continue
            length=len(rows)
            rows.extend([[0]*n for _ in range(degree_p)])
            # The source row indices are disjoint from all target indices.
            # Therefore no in-place update changes a source for this left map.
            for k in range(length):
                for x in outside:
                    a=rows[k][x|first]
                    b=rows[k][x|second]
                    for j in exponents:
                        rows[k+j][x]^=a
                        rows[k+j][x|mask]^=b
            trim(rows)
            length=len(rows)
            rows.extend([[0]*n for _ in range(degree_p)])
            # Right multiplication similarly changes only target columns,
            # disjoint from the source columns selected by the bit masks.
            for k in range(length):
                for x,value in enumerate(rows[k]):
                    delta=((value&selector)<<first)^(((value>>mask)&selector)<<second)
                    for j in exponents:
                        rows[k+j][x]^=delta
            trim(rows)
        identity=[0]*n
        for k,matrix in enumerate(rows):
            row_ranks[k]+=rank(matrix)
            for x,row in enumerate(matrix):
                identity[x]^=row
                while row:
                    low=row&-row
                    row^=low
                    v=low.bit_length()-1
                    h=x^v^g
                    y=x^g
                    columns[k][h][x]|=1<<y
        assert identity==[1<<x for x in range(n)], (t,order,polynomial,g)
    while len(row_ranks)>1 and row_ranks[-1]==0:
        row_ranks.pop()
        columns.pop()
    column_ranks=[sum(rank(block) for block in degree) for degree in columns]
    products=[r*c for r,c in zip(row_ranks,column_ranks)]
    maximum=max(products)
    rank_sum=sum(math.sqrt(p) for p in products)
    return {"t":t,"order":order,"polynomial_mask":polynomial,
            "degree":len(row_ranks)-1,"row_ranks":row_ranks,
            "column_ranks":column_ranks,"maximum_rank_product":maximum,
            "rank_sum_estimate":rank_sum,
            "maximum_product_rho":math.log2(maximum)/(4*t),
            "rank_sum_rho":math.log2(rank_sum)/(2*t)}


def atomic_write(path,value):
    temporary=path.with_suffix(path.suffix+".tmp")
    with temporary.open("w") as out:
        json.dump(value,out,indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    temporary.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths",nargs="+",type=int,default=[2,4,6])
    parser.add_argument("--polynomials",nargs="+",type=int,default=[2,3,6])
    parser.add_argument("--orders",nargs="+",default=["chain","pairs","layers"])
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    results=[]
    for polynomial in args.polynomials:
        for order in args.orders:
            for t in args.lengths:
                start=time.monotonic()
                result=profile(t,order,polynomial)
                result["elapsed_seconds"]=time.monotonic()-start
                results.append(result)
                print(json.dumps(result),flush=True)
                atomic_write(args.output,results)
    assert len(json.loads(args.output.read_text()))==len(results)


if __name__=="__main__":
    main()
