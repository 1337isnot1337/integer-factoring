#!/usr/bin/env python3
"""Exact all-field formal crossing cutoff73 from compatible cores and meshes."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib
import json
from time import monotonic
from check_compatible_compression_budget import verify_cores,compatible_complement
from check_rational_cut_mesh_budget import mesh_weights
from check_interior_jacobian_budget import fraction_record


def kernel_complement(d,h):
    assert 1<=h<=d
    return (2*d-h)**2//4


def formal_table(d,descriptors):
    # Index1 is class-dependent and must never be read from this table.
    F=[d*d,None]+[compatible_complement(d,m,descriptors)[0] for m in range(2,2*d)]
    assert F[2]==d*d//4 and F[-1]==0
    for m in range(2*d-2,1,-1):F[m]=max(F[m],F[m+1])
    return F


def budgets_without_one(d,F,H):
    result=[]
    for j in range(H,d*H+1):
        nums=[None,None]+[H*H*F[m]+j*j*m for m in range(2,2*d)]
        for m in range(2*d-2,1,-1):nums[m]=min(nums[m],nums[m+1])
        nums[1]=nums[2]
        nums[0]=min(H*H*d*d,nums[2])
        result.append(nums)
    return result


def intervals(d,k,b):
    low=(k-1)**2*b+1;high=min(k*k*b,d*d)
    if b==1:
        return [(low,high-1,k-1,'below_square'),(high,high,k,'at_square')]
    return [(low,high,k,'general')]


def parameters(d,b,low,high,h,F):
    kernel=kernel_complement(d,h)
    f=kernel if b==1 else (F[b] if b<len(F) else 0)
    return max(d*d-high,f),d*d-b,max(0,low+f-d*d),kernel


def scaled_class_bound(d,k,b,low,high,h,F,G,H,scale):
    R,C,e,kernel=parameters(d,b,low,high,h,F)
    base,weights=mesh_weights(k,H)
    terms=[base*(2*d*d+e)]
    for j,(s,w) in enumerate(weights,H):
        ss=s*s
        g0=min(Q(G[j-H][0],H*H),kernel+ss)
        gb=Q(G[j-H][b],H*H) if b<len(F) else ss*b
        if b==1:gb=min(gb,kernel+ss)
        if low<=ss*b:gb=g0
        mu=(ss*k*k-1)/(k**4-1);nu=(k*k-ss)/(k**4-1)
        U=max(gb,min(low,ss*b)+mu*R+nu*C)
        V=max(g0,b+nu*R+mu*C)
        terms.append(w*(U+V))
    return sum((x.numerator*scale)//x.denominator for x in terms),len(terms)


def dimension_certificate(d,F,floor):
    scale=10**9;Hlist=(1,4,32)
    Gs={H:budgets_without_one(d,F,H) for H in Hlist}
    assert 3363**2-2*2378**2==1
    assert floor.numerator**2378>d**3363*floor.denominator**2378
    total=base_count=1;counts={H:0 for H in Hlist};refined=[]
    minimum_high=None;below_count=square_count=0
    for k in range(2,d+1):
        for b in range(1,(d*d-1)//((k-1)**2)+1):
            parts=intervals(d,k,b)
            assert parts[0][0]==(k-1)**2*b+1 and parts[-1][1]==min(k*k*b,d*d)
            for z in range(1,len(parts)):assert parts[z-1][1]+1==parts[z][0]
            for low,high,h,label in parts:
                assert low<=high
                total+=1;below_count+=label=='below_square';square_count+=label=='at_square'
                *_,e,_kernel=parameters(d,b,low,high,h,F)
                if Q(2*d*d+e,2*k-1)>=floor:base_count+=1;continue
                history=[]
                for H in Hlist:
                    value,terms=scaled_class_bound(d,k,b,low,high,h,F,Gs[H],H,scale)
                    history.append(dict(mesh_denominator=H,scaled_lower_bound=value,terms=terms))
                    if Q(value,scale)>=floor:
                        counts[H]+=1
                        if H>1:
                            refined.append(dict(k=k,b=b,low=low,high=high,h=h,subcase=label,bounds=history))
                            if minimum_high is None or value<minimum_high[0]:minimum_high=(value,k,b,label,H)
                        break
                else:raise AssertionError((d,k,b,label,history[-1]))
    expected=1+sum((d*d-1)//((k-1)**2) for k in range(2,d+1))+(d-1)
    assert total==expected and total==base_count+sum(counts.values())
    assert below_count==square_count==d-1
    return dict(d=d,scalar_subclasses=total,base_term_classes=base_count,
                mesh_classes=counts,rank_one_below_square=below_count,rank_one_at_square=square_count,
                certified_common_floor=fraction_record(floor),upper_Pell_exponent=[3363,2378],
                downward_rounding_scale=scale,
                minimum_among_refined=dict(numerator=minimum_high[0],denominator=scale,
                  k=minimum_high[1],b=minimum_high[2],subcase=minimum_high[3],mesh_denominator=minimum_high[4]),
                refined_classes=refined,formal_complement_table_without_one=F)


def kernel_controls():
    maxima=products=0
    for d in range(2,74):
        for h in range(1,d+1):
            assert kernel_complement(d,h)==max(s*(d-max(0,s-d+h)) for s in range(d+1))
            maxima+=1
        for a in range(1,d*d+1):
            for p in range(1,d+1):
                if a%p==0 and a//p<=d:
                    assert min(p,a//p)<=isqrt(a)
                    products+=1
    return dict(integer_kernel_maxima=maxima,factor_rank_product_controls=products)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--cores',type=Path,default=Path(__file__).with_name('compatible-compression-cores.json'))
    args=parser.parse_args();started=monotonic()
    descriptors,core_records=verify_cores(args.cores)
    controls=kernel_controls()
    records=[dimension_certificate(d,formal_table(d,descriptors),floor)
             for d,floor in ((72,Q(424)),(73,Q(21583,50)))]
    result=dict(status='PASS',scope='All-field supplied formal crossing limits through73; generic Laurent cut ranks; no factoring exponent improvement',
                previous_formal_cutoff=71,all_field_formal_cutoff=73,
                exact_fitting_or_half_rank_bounds_imported=False,
                small_source_tensor_powers_not_excluded=True,
                core_cross_blocks_verified=sum(len(x['pairs']) for x in core_records),
                kernel_controls=controls,dimensions=records,
                scalar_subclasses=sum(x['scalar_subclasses'] for x in records),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                literal_core_data_sha256=hashlib.sha256(args.cores.read_bytes()).hexdigest(),
                seconds=monotonic()-started)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='dimensions'},indent=2))
    for r in records:print(json.dumps({k:v for k,v in r.items() if k not in ('refined_classes','formal_complement_table_without_one')},indent=2))

if __name__=='__main__':main()
