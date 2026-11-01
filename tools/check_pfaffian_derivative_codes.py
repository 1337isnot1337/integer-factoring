#!/usr/bin/env python3
"""Exact maps, matching/closed blocks, Smith and independent literal controls."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import smith_profile,binary_rank
from check_smith_profile_transfer import toeplitz_rank
from pfaffian_derivative_codes import (blocks,closed_block,code,direct_block,edges,
                                       literal_tensor,matchings,reverse,subset_tensor)
from probe_determinant_inverse_codes import fit

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/pfaffian-derivative-code-controls.json'


def report():
    cases=[];controls=Counter();computed={}
    for n in (2,4,6):
        data=code(n)
        for i,(s,vector) in enumerate(zip(data['subsets'],data['embedding'])):
            assert len(vector)%2==1
            for j,other in enumerate(data['embedding']):
                assert len(set(vector) & set(other))%2==int(i==j);controls['orthogonal_pairings']+=1
            for edge,(u,v) in enumerate(edges(n)):
                actual={w ^ (1 << edge) for w in vector if w >> edge & 1}
                expected=set(matchings(n,s ^ (1 << u) ^ (1 << v))) if (s >> u & 1) and (s >> v & 1) else set()
                assert actual==expected;controls['edge_deletions']+=1
        orders=[];records=[]
        for a,delta,multiplicity in blocks(n):
            key=(a,delta)
            if key not in computed:
                matrix=closed_block(a,delta);direct=direct_block(a,delta)
                assert direct==matrix;controls['direct_block_entries']+=len(matrix)**2
                assert all((p.bit_count()%2)==int(i==j) for i,row in enumerate(matrix) for j,p in enumerate(row))
                degree=max(p.bit_length()-1 for row in matrix for p in row)
                profile,precision=smith_profile(matrix,degree)
                hinges=[]
                for cutoff in range(1,max(profile)+2):
                    rank=toeplitz_rank(matrix,cutoff)
                    assert rank==sum(max(cutoff-v,0) for v in profile)
                    hinges.append(dict(cutoff=cutoff,rank=rank));controls['toeplitz_hinges']+=1
                if a>=2:assert matrix==[[p << 4 for p in row] for row in closed_block(a-2,delta)]
                if a==0 and delta:
                    rows=[sum(((p >> (delta//2))&1) << j for j,p in enumerate(row)) for row in matrix]
                    assert binary_rank(rows)==delta-1;controls['leading_rank_formulas']+=1
                computed[key]=dict(a=a,delta=delta,size=len(matrix),matrix=matrix,
                    orders=profile,precision=precision,independent_toeplitz=hinges)
            record=dict(computed[key],multiplicity=multiplicity);records.append(record)
            orders.extend(list(record['orders'])*multiplicity)
        orders.sort();assert len(orders)==data['q']**2
        numeric=fit(dict(cuts=[dict(retained_orders_t=orders)]*2),data['heights'])
        literal=[]
        if n<=4:
            _,tensor=literal_tensor(n);q=data['q'];pairs=list(product(range(q),repeat=2))
            assert subset_tensor(n)==tensor
            controls['uniform_subset_tensor_entries']+=len(tensor)
            for labels,value in tensor.items():
                x,y,u,v=labels
                assert value.bit_count()%2==int(x==v and y==u)
                assert value.bit_length()-1<=sum(data['heights'][i] for i in labels)
                controls['literal_endpoint_entries']+=1
            for cut in ('zero','infinity'):
                if cut=='zero':matrix=[[tensor[x,y,u,v] for u,v in pairs] for x,y in pairs]
                else:
                    matrix=[[reverse(tensor[x,y,u,v],sum(data['heights'][i] for i in (x,y,u,v)))
                             for y,v in pairs] for x,u in pairs]
                degree=max(p.bit_length()-1 for row in matrix for p in row)
                actual,precision=smith_profile(matrix,degree)
                assert list(actual)==orders
                literal.append(dict(cut=cut,degree=degree,orders=actual,precision=precision))
                controls['complete_literal_profiles']+=1
        cases.append(dict(n=n,q=data['q'],code=data,blocks=records,
            orders=orders,opposite_profile='Derived equal by self-dual coherent compression; literal both-cut controls for n=2,4.',
            numerical_fit=numeric,literal_full_profiles=literal))
    paths=['tools/pfaffian_derivative_codes.py','tools/check_pfaffian_derivative_codes.py',
           'tools/check_block_smith_crossings.py','tools/check_smith_profile_transfer.py',
           'tools/check_cancellation_module_automata.py','tools/probe_determinant_inverse_codes.py']
    return json.loads(json.dumps(dict(scope='Uniform self-dual Pfaffian-code block construction; finite complete profiles at n2,4,6 with direct matching and Toeplitz controls, numerical natural pressure only. No stronger score.',
        factoring_exponent_improved=False,cases=cases,controls=dict(controls),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=report()
    if args.write:RECEIPT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(RECEIPT.read_text())==data
    print(json.dumps(dict(status='PASS',controls=data['controls'],
        sampled_rhos=[c['numerical_fit']['rho'] for c in data['cases']],factoring_exponent_improved=False)))


if __name__=='__main__':main()
