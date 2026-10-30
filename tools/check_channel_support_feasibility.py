#!/usr/bin/env python3
"""Exact finite-carry feasibility and independent arithmetic checks."""
import argparse
import collections
import hashlib
import itertools
import json
import math
import os
import random
from pathlib import Path


def data(a,b):
    assert len(a)==len(b) and a
    m=len(a[0])
    assert all(len(row)==m for row in a)
    bounds=tuple(max(abs(rhs),sum(abs(x) for x in row)) for row,rhs in zip(a,b))
    digits=list(itertools.product((0,1),repeat=m))
    increments=[tuple(sum(x*y for x,y in zip(row,u)) for row in a) for u in digits]
    return m,bounds,digits,increments


def successor(c,increment):
    numer=tuple(x+y for x,y in zip(c,increment))
    if any(x%2 for x in numer):
        return None
    return tuple(x//2 for x in numer)


def solve(a,b):
    m,bounds,digits,increments=data(a,b)
    start=tuple(-x for x in b)
    zero=(0,)*len(b)
    predecessor={start:None}
    queue=collections.deque([start])
    while queue:
        c=queue.popleft()
        if c==zero:
            path=[c]
            labels=[]
            while predecessor[c] is not None:
                previous,digit=predecessor[c]
                labels.append(digit)
                path.append(previous)
                c=previous
            path.reverse()
            labels.reverse()
            witness=[sum((1<<j)*u[i] for j,u in enumerate(labels)) for i in range(m)]
            return {"feasible":True,"witness":witness,"carry_path":path,
                    "digit_path":labels,"bounds":bounds,
                    "visited_count":len(predecessor)}
        for u,increment in zip(digits,increments):
            next_c=successor(c,increment)
            if next_c is None:
                continue
            assert all(abs(x)<=bound for x,bound in zip(next_c,bounds))
            if next_c not in predecessor:
                predecessor[next_c]=(c,u)
                queue.append(next_c)
    return {"feasible":False,"closed_carries":sorted(predecessor),
            "bounds":bounds,"visited_count":len(predecessor)}


def verify(a,b,result):
    m,bounds,digits,increments=data(a,b)
    assert tuple(result["bounds"])==bounds
    start=tuple(-x for x in b)
    zero=(0,)*len(b)
    if result["feasible"]:
        x=result["witness"]
        assert len(x)==m and all(isinstance(v,int) and v>=0 for v in x)
        assert [sum(v*w for v,w in zip(row,x)) for row in a]==b
        path=[tuple(c) for c in result["carry_path"]]
        labels=[tuple(u) for u in result["digit_path"]]
        assert path[0]==start and path[-1]==zero and len(path)==len(labels)+1
        assert len(labels)<=math.prod(2*v+1 for v in bounds)-1
        for c,next_c,u in zip(path,path[1:],labels):
            assert len(u)==m and all(v in (0,1) for v in u)
            assert all(2*q-p==sum(v*w for v,w in zip(row,u))
                       for p,q,row in zip(c,next_c,a))
        return
    closed={tuple(c) for c in result["closed_carries"]}
    assert start in closed and zero not in closed
    for c in closed:
        assert len(c)==len(b) and all(abs(x)<=v for x,v in zip(c,bounds))
        for increment in increments:
            next_c=successor(c,increment)
            assert next_c is None or next_c in closed


def extended_gcd(a,b):
    old_r,r=a,b
    old_s,s=1,0
    old_t,t=0,1
    while r:
        q=old_r//r
        old_r,r=r,old_r-q*r
        old_s,s=s,old_s-q*s
        old_t,t=t,old_t-q*t
    return old_r,old_s,old_t


def two_variable_oracle(a,b,c):
    """Solve a*x+b*y=c independently, by the general gcd solution."""
    if not a and not b:
        return c==0
    if not a:
        return c%b==0 and c//b>=0
    if not b:
        return c%a==0 and c//a>=0
    g,u,v=extended_gcd(abs(a),abs(b))
    if c%g:
        return False
    x0=u*(1 if a>0 else -1)*(c//g)
    y0=v*(1 if b>0 else -1)*(c//g)
    low=high=None
    for base,step in [(x0,b//g),(y0,-a//g)]:
        if step>0:
            bound=-(base//step)  # ceil((-base)/step)
            low=bound if low is None else max(low,bound)
        elif step<0:
            bound=base//(-step)
            high=bound if high is None else min(high,bound)
        elif base<0:
            return False
    return low is None or high is None or low<=high


def bounded_digits(a,b,length):
    _,_,_,increments=data(a,b)
    states={tuple(-x for x in b)}
    for _ in range(length):
        states={next_c for c in states for increment in increments
                if (next_c:=successor(c,increment)) is not None}
    return (0,)*len(b) in states


def support_system(row_edges,column_edges,row_ends,column_ends):
    """All selected edge counts are 1+x; return the system for x."""
    row_vertices=sorted(set(v for e in row_edges for v in e[:2])|set(row_ends))
    col_vertices=sorted(set(v for e in column_edges for v in e[:2])|set(column_ends))
    labels=sorted({e[2] for e in row_edges+column_edges})
    a=[]
    rhs=[]
    nr=len(row_edges)
    nc=len(column_edges)
    for v in row_vertices:
        a.append([int(s==v)-int(t==v) for s,t,_ in row_edges]+[0]*nc)
        rhs.append(int(v==row_ends[0])-int(v==row_ends[1]))
    for v in col_vertices:
        a.append([0]*nr+[int(s==v)-int(t==v) for s,t,_ in column_edges])
        rhs.append(int(v==column_ends[0])-int(v==column_ends[1]))
    for label in labels:
        a.append([int(l==label) for _,_,l in row_edges]
                 +[-int(l==label) for _,_,l in column_edges])
        rhs.append(0)
    return a,[v-sum(row) for row,v in zip(a,rhs)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    digest=hashlib.sha256()
    cases=positive=negative=0
    largest_closed=0
    for a,b,c in itertools.product(range(-3,4),range(-3,4),range(-6,7)):
        matrix=[[a,b]]
        result=solve(matrix,[c])
        verify(matrix,[c],result)
        assert result["feasible"]==two_variable_oracle(a,b,c),(a,b,c,result)
        cases+=1
        positive+=result["feasible"]
        negative+=not result["feasible"]
        largest_closed=max(largest_closed,result["visited_count"])
        digest.update(json.dumps([matrix,[c],result],sort_keys=True).encode())
    assert cases==637
    rng=random.Random(20261030)
    bounded_cases=200
    bounded_feasible=0
    for _ in range(bounded_cases):
        matrix=[[rng.randrange(-2,3) for _ in range(3)] for _ in range(2)]
        rhs=[rng.randrange(-4,5) for _ in range(2)]
        length=3
        brute=any([sum(v*w for v,w in zip(row,x)) for row in matrix]==rhs
                  for x in itertools.product(range(1<<length),repeat=3))
        assert bounded_digits(matrix,rhs,length)==brute
        result=solve(matrix,rhs)
        verify(matrix,rhs,result)
        assert not brute or result["feasible"]
        bounded_feasible+=brute
        largest_closed=max(largest_closed,result["visited_count"])
        digest.update(json.dumps([matrix,rhs,result],sort_keys=True).encode())
    odd_row=[(0,2,"A"),(2,0,"A"),(0,1,"A")]
    even_col=[(0,1,"A"),(1,0,"A")]
    a,b=support_system(odd_row,even_col,(0,1),(0,0))
    parity=solve(a,b)
    verify(a,b,parity)
    assert not parity["feasible"]
    bridge_row=[(0,0,"A"),(0,0,"A"),(1,1,"B"),(1,1,"B"),(0,1,"B")]
    bridge_col=[(0,1,"A"),(0,1,"A"),(1,0,"B"),(1,0,"B")]
    a2,b2=support_system(bridge_row,bridge_col,(0,1),(1,0))
    bridge=solve(a2,b2)
    verify(a2,b2,bridge)
    assert bridge["feasible"]
    counts=[1+x for x in bridge["witness"]]
    assert sum(counts[:5])==sum(counts[5:])
    report={"status":"passed","signed_two_variable_cases":cases,
            "signed_two_variable_feasible":positive,
            "signed_two_variable_infeasible":negative,
            "bounded_digit_multiequation_cases":bounded_cases,
            "bounded_digit_length":3,"bounded_digit_feasible":bounded_feasible,
            "largest_visited_count":largest_closed,"case_sha256":digest.hexdigest(),
            "mixed_endpoint_parity_example":{"matrix":a,"rhs":b,"certificate":parity},
            "two_label_bridge_example":{"matrix":a2,"rhs":b2,
                "positive_edge_counts":counts,"certificate":bridge}}
    if args.output:
        temporary=args.output.with_suffix(args.output.suffix+".tmp")
        with temporary.open("w") as out:
            json.dump(report,out,indent=2)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        temporary.replace(args.output)
    print(json.dumps({k:v for k,v in report.items() if "example" not in k},indent=2))


if __name__=="__main__":
    main()
