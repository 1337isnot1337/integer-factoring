#!/usr/bin/env python3
"""Fixed Path C reproduction; see research/clustered-hits-and-pooled-rank.md."""
from collections import Counter
from fractions import Fraction
from math import gcd,isqrt,prod
from pathlib import Path
import json,lzma,sys,time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from prime_form_factor_base import build_prime_form_factor_base
from leading_power_smooth_relations import compose_with_prime_form
from dynamic_leading_norm_relations import leading_norm_word
from binary_quadratic_forms import inverse_form,ambiguous_factor
from quadratic_ideal_hnf import compose_forms_hnf
from arithmetic_two_group import power
from relation_two_saturation import first_useful_relation_top

OUT=ROOT/'evidence/clustered-power-relations'
old=json.loads((OUT/'receipt.json').read_text())
LIMIT=250000
used=old['base_walk_compositions']+1000  # conservative prior phase-audit allowance
audit_start=used
started=time.monotonic()


def charge():
    global used
    if used>=LIMIT: raise RuntimeError('Predeclared total composition cap reached')
    used+=1


def compose(a,b):
    charge()
    return compose_forms_hnf(a,b)


def extend_basis(basis,row):
    work=list(map(Fraction,row))
    for pivot,vector in sorted(basis.items()):
        scalar=work[pivot]
        if scalar:
            work=[a-scalar*b for a,b in zip(work,vector)]
    if not any(work): return False
    pivot=next(i for i,x in enumerate(work) if x)
    scalar=work[pivot]
    basis[pivot]=[x/scalar for x in work]
    return True


N=1000036000099;Y=19;delta=-N
identity=(1,1,(1-delta)//4)
height=isqrt((-delta)//3)
base=build_prime_form_factor_base(delta,Y)
records={r.prime:(r.root,tuple(r.positive)) for r in base['records']}
primes=sorted(records);indices={p:i for i,p in enumerate(primes)}
assert primes==[5,7,19]
cert=next(x for x in old['certificates'] if x['N']==N)
base_flags=lzma.decompress((ROOT/cert['flags_file']).read_bytes())
results=[]
for case in (x for x in old['cases'] if x['N']==N):
    p=case['prime'];period=case['period'];jump=case['discrete_index'];K=case['K']
    flags=bytes(base_flags[(j*jump)%len(base_flags)] for j in range(period))
    phases=[k for k in range(K,period) if flags[k]]
    basis={};selected=[];selected_witnesses=[];occurrences={};rows=[];verified=set()
    previous_k=None;previous_form=None;previous_row=None
    adjacent=duplicates=0;first_full_rank=None;rank_steps=[]
    for k in phases:
        if previous_k is not None and k==previous_k+1:
            charge()
            form=compose_with_prime_form(previous_form,p,*records[p])
        else:
            form=power(records[p][1],k,identity,compose)
        remaining=form[0];parts=[]
        for q in primes:
            exponent=0
            while remaining%q==0:
                remaining//=q;exponent+=1
            if exponent: parts.append((q,exponent))
        assert remaining==1 and gcd(form[0],N)==1
        word,new=leading_norm_word(N,form,parts,records)
        assert not new and abs(word.get(p,0))<=K-1
        row=[0]*len(primes);row[indices[p]]=k
        for q,e in word.items(): row[indices[q]]-=e
        row=tuple(row)
        if row not in verified:
            represented=identity
            for q,e in word.items():
                generator=records[q][1] if e>=0 else inverse_form(records[q][1])
                represented=compose(represented,power(generator,abs(e),identity,compose))
            assert represented==form
            verified.add(row)
        occurrences.setdefault(row,[]).append(k)
        if previous_k is not None and k==previous_k+1:
            adjacent+=1;duplicates+=row==previous_row
        if extend_basis(basis,row):
            selected.append(list(row));selected_witnesses.append(dict(k=k,form=list(form),word=word))
            rank_steps.append(dict(exponent=k,queries=k-K+1,rank=len(basis)))
            if len(basis)==len(primes) and first_full_rank is None: first_full_rank=k
        rows.append(dict(k=k,row=list(row),leading=form[0]))
        previous_k,previous_form,previous_row=k,form,row
    assert all(len(ks)<=2*K-1 and max(ks)-min(ks)<=2*K-2 for ks in occurrences.values())
    factor=None;extraction=None
    if len(basis)==len(primes):
        _,factor,extraction=first_useful_relation_top([records[q][1] for q in primes],selected,identity,
                    compose,inverse_form,lambda form:ambiguous_factor(form,N),verify_relations=True)
        if factor is not None: assert 1<factor<N and N%factor==0
    result=dict(N=N,cutoff=Y,prime=p,K=K,period=period,smooth_rows=len(rows),
                distinct_rows=len(verified),adjacent_smooth_pairs=adjacent,adjacent_duplicate_rows=duplicates,
                maximum_row_multiplicity=max(map(len,occurrences.values())),multiplicity_bound=2*K-1,
                rank=len(basis),base_size=len(primes),rank_steps=rank_steps,first_full_rank_exponent=first_full_rank,
                selected_rows=selected,selected_witnesses=selected_witnesses,factor=factor,extraction=extraction,rows=rows)
    results.append(result)
    print('row_audit',p,'rows',len(rows),'distinct',len(verified),'adjacentdup',duplicates,'/',adjacent,
          'rank',len(basis),'firstfull',first_full_rank,'factor',factor,'charged',used,flush=True)
receipt=dict(scope='Audit only the already measured smooth phases; skipped phases are supplied by the frozen full-cycle data. No fast selector or uniform rank-window theorem is claimed.',
             prior_walk_compositions=old['base_walk_compositions'],prior_phase_audit_allowance=1000,
             new_audit_compositions=used-audit_start,total_composition_upper_bound=used,composition_cap=LIMIT,
             seconds=time.monotonic()-started,primes=primes,results=results,factoring_exponent='1/2; unchanged')
(OUT/'row-audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('complete',used,'seconds',receipt['seconds'],flush=True)
