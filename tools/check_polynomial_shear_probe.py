#!/usr/bin/env python3
"""Independent polynomial-matrix checks and complete 27-profile shear receipt."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
from polynomial_shear_probe import profile,rank


def multiply(a,b):
    c=0
    while b:
        bit=b&-b
        b^=bit
        c^=a<<(bit.bit_length()-1)
    return c


def full_profile(t,order,p):
    d=1<<t
    n=d*d
    matrix=[[0]*n for _ in range(n)]
    for x,y,v,u in itertools.product(range(d),repeat=4):
        value=1
        for bit in range(t):
            xb,yb,vb,ub=((q>>bit)&1 for q in (x,y,v,u))
            local=int(xb==vb and yb==ub)^(3*int(xb==ub and yb==vb))
            value=multiply(value,local)
        matrix[x*d+y][v*d+u]=value
    if order=="pairs":
        windows=list(range(0,t-1,2))
    elif order=="layers":
        windows=list(range(0,t-1,2))+list(range(1,t-1,2))
    else:
        windows=list(range(t-1))
    for start in windows:
        mask=3<<start
        first=1<<start
        second=2<<start
        arrows=[]
        for x,y in itertools.product(range(d),repeat=2):
            if x&mask or y&mask:
                continue
            arrows.append((x*d+(y|mask),(x|first)*d+(y|second)))
            arrows.append(((x|mask)*d+y,(x|second)*d+(y|first)))
        for target,source in arrows:
            for col in range(n):
                matrix[target][col]^=multiply(p,matrix[source][col])
        for target,source in arrows:
            for row in range(n):
                matrix[row][source]^=multiply(p,matrix[row][target])
    degree=max(v.bit_length()-1 for row in matrix for v in row)
    rowcuts=[[0]*n for _ in range(degree+1)]
    colcuts=[[0]*n for _ in range(degree+1)]
    for x,y,v,u in itertools.product(range(d),repeat=4):
        value=matrix[x*d+y][v*d+u]
        assert (value.bit_count()%2)==int(x==v and y==u)
        while value:
            bit=value&-value
            value^=bit
            k=bit.bit_length()-1
            rowcuts[k][x*d+y]|=1<<(v*d+u)
            colcuts[k][x*d+u]|=1<<(y*d+v)
    return [rank(a) for a in rowcuts],[rank(a) for a in colcuts]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    independent=0
    for p,order,t in itertools.product([2,3,6],["chain","pairs","layers"],[2,3]):
        direct=full_profile(t,order,p)
        reduced=profile(t,order,p)
        assert direct==(reduced["row_ranks"],reduced["column_ranks"]),(p,order,t)
        independent+=1
    profiles=[]
    for p,order,t in itertools.product([2,3,6],["chain","pairs","layers"],[2,4,6]):
        result=profile(t,order,p)
        products=[r*c for r,c in zip(result["row_ranks"],result["column_ranks"])]
        lower=sum(math.isqrt(q) for q in products)
        upper=sum(math.isqrt(q)+(math.isqrt(q)**2!=q) for q in products)
        assert lower>1<<(2*t)
        result["integer_rank_sum_lower"]=lower
        result["integer_rank_sum_upper"]=upper
        profiles.append(result)
    expected={(p,order,t) for p,order,t in itertools.product([2,3,6],["chain","pairs","layers"],[2,4,6])}
    assert len(profiles)==27 and {(r["polynomial_mask"],r["order"],r["t"]) for r in profiles}==expected
    payload=[{k:r[k] for k in ["t","order","polynomial_mask","degree","row_ranks","column_ranks"]} for r in profiles]
    report={"status":"passed","independent_full_polynomial_profiles":independent,
            "complete_screen_profiles":len(profiles),"all_rank_sums_exceed_label_dimension_squared":True,
            "exact_profiles_sha256":hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()).hexdigest(),
            "profiles":profiles,"scope":"Finite polynomial-shear controls; no improving seed or asymptotic bound."}
    temporary=args.output.with_suffix(args.output.suffix+".tmp")
    with temporary.open("w") as out:
        json.dump(report,out,indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    temporary.replace(args.output)
    persisted=json.loads(args.output.read_text())
    assert len(persisted["profiles"])==27
    print(json.dumps({k:v for k,v in report.items() if k!="profiles"},indent=2))


if __name__=="__main__":
    main()
