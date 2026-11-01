#!/usr/bin/env python3
"""Build the small mu N Delta core from exact weighted-Weyl block dumps."""
import argparse
from collections import defaultdict
import json
import struct
from time import monotonic

from witt_exchange_jets import multiply
from weighted_energy_pencil import rank_two


def poly_degree(f):
    return f.bit_length()-1


def poly_mod(a,b):
    while a and poly_degree(a)>=poly_degree(b):
        a^=b<<(poly_degree(a)-poly_degree(b))
    return a


def poly_gcd(a,b):
    while b:a,b=b,poly_mod(a,b)
    return a


def poly_multiply(a,b):
    result=0
    while b:
        bit=b&-b;result^=a<<(bit.bit_length()-1);b^=bit
    return result


def poly_divide(a,b):
    result=0
    while a and poly_degree(a)>=poly_degree(b):
        shift=poly_degree(a)-poly_degree(b);result^=1<<shift;a^=b<<shift
    assert a==0
    return result


def minimal_polynomial(A):
    n=len(A);basis={};power=[1<<i for i in range(n)]
    for k in range(n+1):
        value=sum(row<<(i*n) for i,row in enumerate(power));witness=1<<k
        while value:
            pivot=value.bit_length()-1
            if pivot not in basis:
                basis[pivot]=(value,witness);break
            vector,mask=basis[pivot];value^=vector;witness^=mask
        if not value:return witness
        power=multiply(power,A)
    raise AssertionError("minimal polynomial exceeded dimension")


def read_header(f):
    assert f.read(8)==b"WBIN001\0"
    h,d,blocks,pw,nw,cut=struct.unpack("<6I",f.read(24));assert pw==1
    weights=struct.unpack("<%dI"%nw,f.read(4*nw))
    labels=[struct.unpack("<%dI"%nw,f.read(4*nw)) for _ in range(d)]
    return h,d,blocks,weights,labels


def main():
    parser=argparse.ArgumentParser();parser.add_argument("prefix");parser.add_argument("--output");args=parser.parse_args()
    started=monotonic()
    with open(args.prefix+".row.bin","rb") as f:
        h,d,blocks,weights,labels=read_header(f)
        sums={};sum_labels=[];pair_sum=[-1]*(d*d);disjoint=[False]*(d*d)
        for a,aa in enumerate(labels):
            for b,bb in enumerate(labels):
                if any(x&y for x,y in zip(aa,bb)):continue
                disjoint[a*d+b]=True;sigma=tuple(x+y for x,y in zip(aa,bb))
                if sigma not in sums:sums[sigma]=len(sums);sum_labels.append(sigma)
        for a,aa in enumerate(labels):
            for b,bb in enumerate(labels):pair_sum[a*d+b]=sums.get(tuple(x+y for x,y in zip(aa,bb)),-1)
        core=[0]*len(sums)
        for _ in range(blocks):
            n=struct.unpack("<I",f.read(4))[0];pairs=struct.unpack("<%dI"%n,f.read(4*n));polys=struct.unpack("<%dQ"%(n*n),f.read(8*n*n))
            for i,code_i in enumerate(pairs):
                sigma=pair_sum[code_i]
                if sigma<0:continue
                for j,code_j in enumerate(pairs):
                    if disjoint[code_j] and polys[i*n+j]&1:core[sigma]^=1<<pair_sum[code_j]
    parity=defaultdict(list)
    for i,sigma in enumerate(sum_labels):parity[tuple(x&1 for x in sigma)].append(i)
    report=[];annihilator=1
    for group in parity.values():
        indices={i:j for j,i in enumerate(group)}
        small=[]
        for i in group:
            row=0;value=core[i]
            while value:
                bit=value&-value;j=bit.bit_length()-1;assert j in indices
                row|=1<<indices[j];value^=bit
            small.append(row)
        polynomial=minimal_polynomial(small)
        annihilator=poly_multiply(poly_divide(annihilator,poly_gcd(annihilator,polynomial)),polynomial)
        report.append({"dimension":len(group),"rank":rank_two(small),"minimal_polynomial":hex(polynomial),"minimal_degree":poly_degree(polynomial)})
    result={"energy":h,"labels":d,"core_dimension":len(sums),"core_rank":rank_two(core),
            "core_minimal_polynomial":hex(annihilator),"core_minimal_degree":poly_degree(annihilator),
            "central_annihilator":hex(annihilator<<1),"word_algebra_dimension_upper":4*(poly_degree(annihilator)+1),
            "parity_blocks":report,"row_rank_for_nonzero_identity_coefficient_lower":d*d-2*len(sums),
            "seconds":round(monotonic()-started,3)}
    if args.output:open(args.output,"w").write(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
