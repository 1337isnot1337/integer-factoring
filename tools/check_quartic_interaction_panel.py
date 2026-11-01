#!/usr/bin/env python3
"""Exact fixed non-Gaussian interaction panel and all-scale domination."""
import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

from cube_homogeneous_gaussian import (ZERO,monomial,gaussian_tensor,
                                      multiply,axis_matrix,serializable)
from check_block_smith_crossings import smith_profile
from check_smith_profile_transfer import toeplitz_rank
from check_fixed_cube_projection_floor import partition_inequality

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/quartic-interaction-panel.json'


def baseline_kernel():
    x2=monomial(2,0,0);y2=monomial(0,2,0);xy=monomial(1,1,0)
    steps=[ZERO,y2 ^ xy,x2 ^ xy,xy]
    return [[steps[p ^ q] if i==j else ZERO
             for q in range(4) for j in range(2)]
            for p in range(4) for i in range(2)]


def integer_polynomial(polynomial,variable):
    result=0
    for e in polynomial:result ^= 1 << e[variable]
    return result


def complete_profiles(tensor):
    profiles=[];hinges=0
    for physical_axis,variable in ((2,0),(1,1)):
        matrix=axis_matrix(tensor,physical_axis)
        values=[[integer_polynomial(p,variable) for p in row] for row in matrix]
        degree=max(p.bit_length()-1 for row in values for p in row)
        orders,precision=smith_profile(values,degree)
        assert len(orders)==16 and max(orders)<precision
        checks=[]
        for cutoff in range(1,max(orders)+2):
            rank=toeplitz_rank(values,cutoff)
            predicted=sum(max(cutoff-n,0) for n in orders)
            assert rank==predicted
            checks.append(dict(cutoff=cutoff,rank=rank));hinges+=1
        profiles.append(dict(physical_axis=physical_axis,variable=variable,
                             orders=orders,precision=precision,degree_bound=degree,
                             toeplitz_checks=checks))
    return profiles,hinges


def report():
    kernel=baseline_kernel();seed=gaussian_tensor(kernel)
    modifiers=[('x4+y4',monomial(4,0,0) ^ monomial(0,4,0)),
               ('x4+x2y2',monomial(4,0,0) ^ monomial(2,2,0)),
               ('y4+x2y2',monomial(0,4,0) ^ monomial(2,2,0))]
    panel=[('Gaussian baseline',0,ZERO)]+[
        (f'S{S}:{name}',S,h) for S in (153,149) for name,h in modifiers]
    records=[];endpoint_checks=0;hinges=0
    for name,support,h in panel:
        tensor=tuple(p ^ multiply(h,seed[mask ^ support])
                     if support and mask & support==support else p
                     for mask,p in enumerate(seed))
        for mask,p in enumerate(tensor):
            expected=int(all(((mask >> (2*q+i))&1)==((mask >> (2*(q ^ 3)+i))&1)
                             for q in range(4) for i in range(2)))
            assert len(p)%2==expected
            assert all(sum(e)==mask.bit_count() for e in p)
            endpoint_checks+=1
        assert tensor[0]==seed[0]
        assert all(tensor[(1 << i) | (1 << j)]==seed[(1 << i) | (1 << j)]
                   for i,j in combinations(range(8),2))
        if support:
            assert tensor[support] ^ seed[support]==h and h
            assert len(h)%2==0 and all(sum(e)==4 for e in h)
        profiles,count=complete_profiles(tensor);hinges+=count
        paired=sorted(a+b for a in profiles[0]['orders'] for b in profiles[1]['orders'])
        dominance=None
        if records:
            baseline=records[0]['paired_orders']
            c=b=0;gaps=[]
            for candidate,original in zip(paired,baseline):
                c+=candidate;b+=original;gaps.append(b-c)
            assert min(gaps)>=0
            dominance=dict(baseline_prefix_minus_candidate_prefix=gaps,
                           all_prefixes_nonnegative=True,
                           scope='Candidate Zrow*Zcolumn is >= baseline product at every 0<t<1')
        records.append(dict(name=name,quartic_support=support,interaction=h,
                            tensor_entries=tensor,non_Gaussian=bool(support),
                            profiles=profiles,paired_orders=paired,dominance=dominance))
    expected=sorted(sum(values) for values in product((0,1),(0,1),(0,3),(0,3)))
    assert all(list(p['orders'])==expected for p in records[0]['profiles'])
    sources=['tools/check_quartic_interaction_panel.py','tools/cube_homogeneous_gaussian.py',
             'tools/check_block_smith_crossings.py','tools/check_smith_profile_transfer.py',
             'tools/check_fixed_cube_projection_floor.py']
    return serializable(dict(
        scope='Six fixed homogeneous non-Gaussian four-fermion interactions; complete profile products are all-scale dominated by the Gaussian baseline, so common-type envelope rho>=3/4. No universal non-Gaussian or factoring lower bound.',
        cases=records,partition_certificate=partition_inequality(),
        counts=dict(perturbations=6,endpoint_homogeneity_entries=endpoint_checks,
                    complete_profiles=14,Toeplitz_hinge_checks=hinges,
                    paired_prefix_comparisons=6*256),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],counts=result['counts'])))


if __name__=='__main__':main()
