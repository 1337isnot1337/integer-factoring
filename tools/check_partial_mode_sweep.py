#!/usr/bin/env python3
"""Compare the unique-path implementation with direct polynomial matrix products."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic


def pmul(a,b):
    c=0
    while b:
        bit=b&-b;c^=a<<(bit.bit_length()-1);b^=bit
    return c


def rank(rows):
    basis={}
    for a in rows:
        while a:
            p=a.bit_length()-1
            if p not in basis:
                basis[p]=a;break
            a^=basis[p]
    return len(basis)


def labels_at(h):
    weights=[];w=1
    while w<=h:
        weights.append(w);w*=3
    labels=[];a=[0]*len(weights)
    def rec(j,n):
        if j==0:
            a[0]=n;labels.append(tuple(a));return
        for k in range(n//weights[j]+1):
            a[j]=k;rec(j-1,n-k*weights[j])
    rec(len(weights)-1,h)
    return weights,labels


def direct_profile(h,order):
    weights,labels=labels_at(h);d=len(labels);degree=h*len(order)
    row_groups=defaultdict(list);column_groups=defaultdict(list)
    for a,alpha in enumerate(labels):
        for y,second in enumerate(labels):
            row_groups[tuple(x-z for x,z in zip(alpha,second))].append(a*d+y)
            column_groups[tuple(x+z for x,z in zip(alpha,second))].append(a*d+y)
    column_groups=list(column_groups.values());location={code:(bi,i) for bi,g in enumerate(column_groups) for i,code in enumerate(g)}
    column_polys=[[[0]*len(g) for _ in g] for g in column_groups]
    binom=[]
    for n in range(h+1):
        binom.append(sum(1<<k for k in range(n+1) if k&~n==0))
    mode_cache={}
    def mode(j,a,y,beta):
        key=(j,a,y,beta)
        if key not in mode_cache:
            p=0
            for k in range(min(a,y)+1):
                v=a-k
                if k&~y==0 and 0<=v<=beta and v&~beta==0:
                    p^=binom[weights[j]*k]
            mode_cache[key]=p
        return mode_cache[key]
    row_ranks=[0]*(degree+1);direct_multiplications=0;entries=0
    for pairs in row_groups.values():
        n=len(pairs);current=[[int(i==j) for j in range(n)] for i in range(n)]
        for edge in order:
            gate=[[0]*n for _ in range(n)]
            for i,rc in enumerate(pairs):
                a,y=divmod(rc,d);alpha,second=labels[a],labels[y]
                for j,cc in enumerate(pairs):
                    beta,u=divmod(cc,d);first_in,second_in=labels[beta],labels[u]
                    if any(alpha[k]!=first_in[k] or second[k]!=second_in[k] for k in range(len(weights)) if k not in (edge,edge+1)):
                        continue
                    gate[i][j]=pmul(mode(edge,alpha[edge],second[edge],first_in[edge]),mode(edge+1,alpha[edge+1],second[edge+1],first_in[edge+1]))
            output=[[0]*n for _ in range(n)]
            for i in range(n):
                for k in range(n):
                    if gate[i][k]:
                        for j in range(n):
                            if current[k][j]:
                                output[i][j]^=pmul(gate[i][k],current[k][j]);direct_multiplications+=1
            current=output
        row_matrices=[[0]*n for _ in range(degree+1)]
        for i,rc in enumerate(pairs):
            a,y=divmod(rc,d)
            for j,cc in enumerate(pairs):
                beta,u=divmod(cc,d);p=current[i][j]
                assert p.bit_count()%2==(i==j)
                bi,r=location[a*d+u];bj,c=location[y*d+beta];assert bi==bj
                column_polys[bi][r][c]=p
                while p:
                    bit=p&-p;k=bit.bit_length()-1
                    assert k<=degree
                    row_matrices[k][i]|=1<<j;p^=bit
                entries+=1
        for k,matrix in enumerate(row_matrices):row_ranks[k]+=rank(matrix)
    column_ranks=[0]*(degree+1)
    for matrix in column_polys:
        n=len(matrix);planes=[[0]*n for _ in range(degree+1)]
        for i in range(n):
            for j in range(n):
                p=matrix[i][j]
                while p:
                    bit=p&-p;k=bit.bit_length()-1;planes[k][i]|=1<<j;p^=bit
        for k,plane in enumerate(planes):column_ranks[k]+=rank(plane)
    return {"energy":h,"labels":d,"application_order":order,"row_ranks":row_ranks,"column_ranks":column_ranks,"direct_nonzero_polynomial_products":direct_multiplications,"supported_entries":entries}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output");args=parser.parse_args();started=monotonic()
    here=Path(__file__).resolve().parent;source=here/"explore_partial_mode_sweep.cpp"
    reports=[]
    with tempfile.TemporaryDirectory(prefix="partial-sweep-") as temp:
        exe=Path(temp)/"sweep"
        subprocess.run(["g++","-O2","-std=c++17",str(source),"-o",str(exe)],check=True)
        for h,order in [(9,[0,1]),(27,[0,1,2]),(27,[0,2,1])]:
            direct=direct_profile(h,order)
            native=json.loads(subprocess.run([str(exe),str(h),"".join(map(str,order))],check=True,capture_output=True,text=True).stdout)
            assert native["row_ranks"]==direct["row_ranks"]
            assert native["column_ranks"]==direct["column_ranks"]
            direct["native_endpoint_checks"]=native["endpoint_checks"]
            direct["rank_sum"]=native["rank_sum"]
            direct["finite_target"]=native["finite_target"]
            reports.append(direct)
    result={"profiles":reports,"source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),source)},"seconds":round(monotonic()-started,3)}
    text=json.dumps(result,indent=2)+"\n"
    if args.output:Path(args.output).write_text(text)
    print(json.dumps({"profiles":[{k:v for k,v in p.items() if k not in ("row_ranks","column_ranks")} for p in reports],"seconds":result["seconds"]},indent=2))


if __name__=="__main__":main()
