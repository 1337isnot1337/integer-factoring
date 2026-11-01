#!/usr/bin/env python3
"""Exact powers of a constant partial-mode sweep, retaining both cut ranks."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from time import monotonic

from check_mismatched_weyl_spectrum import product,rank
from check_partial_mode_sweep import labels_at


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--energy",type=int,default=63)
    parser.add_argument("--order",default="012");parser.add_argument("--maximum-power",type=int,default=256)
    parser.add_argument("--output");args=parser.parse_args();started=monotonic()
    h=args.energy;weights,labels=labels_at(h);d=len(labels);order=list(map(int,args.order))
    assert sorted(order)==list(range(len(weights)-1))
    row_groups=defaultdict(list);col_groups=defaultdict(list)
    for a,alpha in enumerate(labels):
        for y,second in enumerate(labels):
            row_groups[tuple(x-z for x,z in zip(alpha,second))].append(a*d+y)
            col_groups[tuple(x+z for x,z in zip(alpha,second))].append(a*d+y)
    row_groups=list(row_groups.values());col_groups=list(col_groups.values())
    locations={code:(bi,j) for bi,group in enumerate(col_groups) for j,code in enumerate(group)}
    blocks=[]
    for pairs in row_groups:
        n=len(pairs);current=[1<<i for i in range(n)]
        for edge in order:
            gate=[]
            for rc in pairs:
                a,y=divmod(rc,d);alpha,second=labels[a],labels[y];row=0
                for j,cc in enumerate(pairs):
                    beta,u=divmod(cc,d);first_in,second_in=labels[beta],labels[u]
                    if any(alpha[k]!=first_in[k] or second[k]!=second_in[k] for k in range(len(weights)) if k not in (edge,edge+1)):
                        continue
                    if all(alpha[k]&~(first_in[k]+second[k])==0 for k in (edge,edge+1)):
                        row|=1<<j
                gate.append(row)
            current=product(gate,current)
        blocks.append(current)
    powers=[];power=1;previous=None
    while power<=args.maximum_power:
        row_rank=sum(map(rank,blocks));col_matrices=[[0]*len(g) for g in col_groups]
        for pairs,matrix in zip(row_groups,blocks):
            for i,row in enumerate(matrix):
                a,y=divmod(pairs[i],d)
                while row:
                    bit=row&-row;j=bit.bit_length()-1;beta,u=divmod(pairs[j],d)
                    bi,r=locations[a*d+u];bj,c=locations[y*d+beta];assert bi==bj
                    col_matrices[bi][r]^=1<<c;row^=bit
        col_rank=sum(map(rank,col_matrices))
        powers.append({"power":power,"row_rank":row_rank,"column_rank":col_rank,"rank_product":row_rank*col_rank})
        if previous==row_rank:break
        previous=row_rank;blocks=[product(A,A) for A in blocks];power*=2
    if h==63 and order==[0,1,2]:assert powers[0]["row_rank"]==5195 and powers[0]["column_rank"]==2400
    result={"energy":h,"labels":d,"application_order":order,"powers":powers,
            "row_rank_stabilized":len(powers)>1 and powers[-1]["row_rank"]==powers[-2]["row_rank"],
            "finite_rank_sum_target":d**(2**.5),"seconds":round(monotonic()-started,3),
            "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(__file__).with_name("check_partial_mode_sweep.py"),Path(__file__).with_name("check_mismatched_weyl_spectrum.py"))}}
    text=json.dumps(result,indent=2)+"\n"
    if args.output:Path(args.output).write_text(text)
    print(text,end="")


if __name__=="__main__":main()
