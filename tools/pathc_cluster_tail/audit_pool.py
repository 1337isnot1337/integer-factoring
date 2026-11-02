#!/usr/bin/env python3
"""Fixed Path C reproduction; see research/clustered-hits-and-pooled-rank.md."""
from fractions import Fraction
from math import gcd,isqrt
from pathlib import Path
import json,lzma,sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from prime_form_factor_base import build_prime_form_factor_base
from dynamic_leading_norm_relations import leading_norm_word
from quadratic_ideal_hnf import compose_forms_hnf
from binary_quadratic_forms import inverse_form,ambiguous_factor
from arithmetic_two_group import power
from weighted_norm_relations import least_power_above
from relation_two_saturation import first_useful_relation_top

OUT=ROOT/'evidence/clustered-power-relations'
orbit=json.loads((OUT/'receipt.json').read_text())
audit=json.loads((OUT/'row-audit.json').read_text())
used=audit['total_composition_upper_bound'];before=used;LIMIT=250000

def compose(a,b):
    global used
    if used>=LIMIT: raise RuntimeError('Predeclared composition cap reached')
    used+=1
    return compose_forms_hnf(a,b)

def extend(basis,row):
    v=list(map(Fraction,row))
    for j,w in sorted(basis.items()):
        if v[j]:
            c=v[j];v=[a-c*b for a,b in zip(v,w)]
    if not any(v):return False
    j=next(i for i,a in enumerate(v) if a);c=v[j]
    basis[j]=[a/c for a in v];return True

N=1000036000099;Y=19;delta=-N;identity=(1,1,(1-delta)//4)
base=build_prime_form_factor_base(delta,Y)
records={r.prime:(r.root,tuple(r.positive)) for r in base['records']}
primes=sorted(records);H=isqrt(N//3)
five=next(c for c in audit['results'] if c['prime']==5)
h=five['period'];j7=next(c['discrete_index'] for c in orbit['cases'] if c['N']==N and c['prime']==7)
combined=[x+y for x,y in zip(five['selected_rows'][1],five['selected_rows'][2])]
assert combined[2]==1
j19=(-combined[0]-combined[1]*j7)%h
assert power(records[5][1],j19,identity,compose)==records[19][1]
base_flags=lzma.decompress((ROOT/next(c['flags_file'] for c in orbit['certificates'] if c['N']==N)).read_bytes())
jumps={5:1,7:j7,19:j19};K={p:least_power_above(p,H) for p in primes}
cache={(c['prime'],r['k']):r['row'] for c in audit['results'] for r in c['rows']}

def row_at(p,k):
    if (p,k) in cache:return cache[p,k]
    form=power(records[p][1],k,identity,compose)
    a=form[0];parts=[]
    for q in primes:
        e=0
        while a%q==0:a//=q;e+=1
        if e:parts.append((q,e))
    assert a==1
    word,new=leading_norm_word(N,form,parts,records);assert not new
    v=[0]*len(primes);v[primes.index(p)]=k
    represented=identity
    for q,e in word.items():
        v[primes.index(q)]-=e
        g=records[q][1] if e>=0 else inverse_form(records[q][1])
        represented=compose(represented,power(g,abs(e),identity,compose))
    assert represented==form
    cache[p,k]=v;return v

events=[];firsts={}
for p in primes:
    local=bytes(base_flags[(j*jumps[p])%h] for j in range(h))
    firsts[p]=next(j for j in range(h) if local[(K[p]+j)%h])
    for offset in range(h-K[p]):
        if local[K[p]+offset]:events.append((offset,p,K[p]+offset))
events.sort()
basis={};selected=[];rank_steps=[];stop=None
for offset,p,k in events:
    row=row_at(p,k)
    if extend(basis,row):
        selected.append(row)
        rank_steps.append(dict(prime=p,exponent=k,offset=offset,rank=len(basis),
                               round_robin_norm_queries=len(primes)*offset+primes.index(p)+1))
        if len(basis)==len(primes):stop=rank_steps[-1];break
assert stop
_,factor,cert=first_useful_relation_top([records[p][1] for p in primes],selected,identity,
                   compose,inverse_form,lambda f:ambiguous_factor(f,N),verify_relations=True)
assert factor and 1<factor<N and N%factor==0
baseline=[];one_each=[]
for p in primes:
    k=K[p]+firsts[p];one_each.append(row_at(p,k))
    baseline.append(dict(prime=p,K=K[p],first_offset=firsts[p],norm_queries=firsts[p]+1))
rank={}
for row in one_each:extend(rank,row)
assert len(rank)==len(primes)
result=dict(scope='Finite derived-orbit pooled-rank diagnostic. Norm-query counts describe schedules on the complete frozen flags, not a new timed native selector. No uniform rank-window bound follows.',
            N=N,cutoff=Y,period=h,primes=primes,derived_generator_indices=jumps,
            index19_witness_relation=combined,round_robin_rank_steps=rank_steps,selected_rows=selected,
            pooled_round_robin_norm_queries=stop['round_robin_norm_queries'],factor=factor,extraction=cert,
            one_first_row_per_target=baseline,one_per_target_total_norm_queries=sum(x['norm_queries'] for x in baseline),
            one_per_target_rows=one_each,additional_audit_compositions=used-before,
            total_composition_upper_bound=used,composition_cap=LIMIT,factoring_exponent='1/2; unchanged')
(OUT/'pool-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['round_robin_rank_steps','pooled_round_robin_norm_queries','one_per_target_total_norm_queries','factor','total_composition_upper_bound']}))
