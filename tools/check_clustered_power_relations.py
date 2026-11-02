#!/usr/bin/env python3
"""Replay frozen Path C tail data and exact relation/order certificates.

The full producers also check every arithmetic phase. This faster replay
checks all frozen statistical quantities, arithmetic samples and hit rows.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd,isqrt,prod
from pathlib import Path
import json,lzma

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form,ambiguous_factor
from check_carry_fourier_elimination import prime_factors
from check_cofactor_orbit_moments import cyclic_window_counts,sufficient_moment_window,maximum_zero_run
from dynamic_leading_norm_relations import leading_norm_word
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'evidence/clustered-power-relations'


def extend(basis,row):
    work=list(map(Fraction,row))
    for pivot,vector in sorted(basis.items()):
        if work[pivot]:
            c=work[pivot];work=[x-c*y for x,y in zip(work,vector)]
    if not any(work): return False
    pivot=next(i for i,x in enumerate(work) if x);c=work[pivot]
    basis[pivot]=[x/c for x in work]
    return True


def model(N,Y):
    delta=-N if N%4==3 else -4*N
    identity=(1,delta%2,(delta%2-delta)//4)
    base=build_prime_form_factor_base(delta,Y)
    records={r.prime:(r.root,tuple(r.positive)) for r in base['records']}
    return identity,records,isqrt((-delta)//3)


def word_value(word,records,identity):
    result=identity
    for p,e in word.items():
        g=records[p][1] if e>=0 else inverse_form(records[p][1])
        result=compose_forms_hnf(result,power(g,abs(e),identity,compose_forms_hnf))
    return result


def check_statistics(stored,flags,N,K):
    h,s=len(flags),sum(flags)
    assert stored['period']==h and stored['good']==s
    assert stored['first_from_K']==next(j for j in range(h) if flags[(K+j)%h])
    assert stored['maximum_zero_run']==maximum_zero_run(flags)
    assert stored['lag_pair_counts']==[sum(flags[j]*flags[(j+lag)%h] for j in range(h)) for lag in range(17)]
    for tail in stored['tails']:
        W=(tail['target_mean']*h+s-1)//s
        histogram=Counter(cyclic_window_counts(flags,W))
        assert (tail['window'],tail['mean_numerator'],tail['mean_denominator'])==(W,W*s,h)
        assert tail['empty_phases']==histogram.get(0,0) and tail['phases']==h
    q0=(N.bit_length()+1).bit_length()
    assert [(m['q'],m['constant']) for m in stored['moments']]==[(q,C) for q in (q0,q0*q0) for C in (1,2)]
    for m in stored['moments']:
        q,C=m['q'],m['constant'];D=(h+s-1)//s
        W=sufficient_moment_window(2*N,D,q,C)
        histogram=Counter(cyclic_window_counts(flags,W))
        numerator=sum(n*(s*W-h*a)**(2*q) for a,n in histogram.items() if h*a<s*W)
        unit=h**(q+1)*(q*W*s)**q
        assert (m['window'],m['density_denominator'],m['period_bound'])==(W,D,2*N)
        assert (m['lower_tail_scaled_sum'],m['constant_one_scaled_bound'])==(numerator,unit)
        assert m['window_covers_period']==(W>=h) and m['minimum_hits']==min(histogram)
        required=m['required_integer_constant']
        assert required>=1 and numerator<=required**q*unit
        assert required==1 or numerator>(required-1)**q*unit
        assert m['holds']==(numerator<=C**q*unit)
        if m['holds']: assert min(histogram)>=1


def run():
    orbit=json.loads((DATA/'receipt.json').read_text())
    rows=json.loads((DATA/'row-audit.json').read_text())
    pool=json.loads((DATA/'pool-audit.json').read_text())
    assert orbit['base_walk_compositions']==sum(c['period'] for c in orbit['certificates'])==219306
    assert rows['prior_walk_compositions']==orbit['base_walk_compositions']
    assert rows['total_composition_upper_bound']==(rows['prior_walk_compositions']+
                rows['prior_phase_audit_allowance']+rows['new_audit_compositions'])
    assert pool['total_composition_upper_bound']==rows['total_composition_upper_bound']+pool['additional_audit_compositions']
    for name,digest in orbit['source_sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest()==digest
    flags_by_N={};models={}
    for certificate in orbit['certificates']:
        N,Y,p=certificate['N'],certificate['cutoff'],certificate['base_prime']
        identity,records,H=model(N,Y);models[N]=(identity,records,H)
        path=ROOT/certificate['flags_file'];compressed=path.read_bytes();flags=lzma.decompress(compressed)
        assert sha256(compressed).hexdigest()==certificate['archive_sha256']
        assert sha256(flags).hexdigest()==certificate['flags_sha256']
        h=certificate['period'];assert len(flags)==h and flags[0] and h<=2*N
        g=records[p][1];assert power(g,h,identity,compose_forms_hnf)==identity
        assert all(power(g,h//q,identity,compose_forms_hnf)!=identity for q in prime_factors(h))
        assert all(flags[j]==flags[-j] for j in range(h))
        flags_by_N[N]=flags
    local_flags={}
    for case in orbit['cases']:
        N,p=case['N'],case['prime'];identity,records,H=models[N]
        base=flags_by_N[N];h=case['period'];jump=case['discrete_index']
        assert h==len(base)//gcd(len(base),jump)
        assert power(records[case['base_prime']][1],jump,identity,compose_forms_hnf)==records[p][1]
        flags=bytes(base[(j*jump)%len(base)] for j in range(h));local_flags[N,p]=flags
        entries=bytes(flags[j]*(1-flags[j-1]) for j in range(h))
        assert (case['smooth_classes'],case['clusters'])==(sum(flags),sum(entries))
        assert case['maximum_hit_run']==maximum_zero_run(bytes(1-x for x in flags))
        K=least_power_above(p,H);assert K==case['K']
        check_statistics(case['raw'],flags,N,K);check_statistics(case['entry'],entries,N,K)
        P=prod(records)
        for sample in case['direct_arithmetic_samples']:
            form=power(records[p][1],sample['phase'],identity,compose_forms_hnf)
            divisor=gcd(form[0],N)
            good=int(1<divisor<N or pow(P,form[0].bit_length(),form[0])==0)
            assert list(form)==sample['form'] and good==sample['good']==flags[sample['phase']]
    row_cache={};adjacent_total=duplicate_total=0
    for case in rows['results']:
        N,p,K=case['N'],case['prime'],case['K'];identity,records,H=models[N]
        primes=sorted(records);flags=local_flags[N,p];basis={};selected=[];rank_steps=[];occurrences={}
        assert [r['k'] for r in case['rows']]==[k for k in range(K,len(flags)) if flags[k]]
        adjacent=duplicates=0;previous=None
        for item in case['rows']:
            k=item['k'];row=item['row'];form=power(records[p][1],k,identity,compose_forms_hnf)
            word={q:(k if q==p else 0)-r for q,r in zip(primes,row)}
            assert prod(q**abs(e) for q,e in word.items())==item['leading']==form[0]<=H
            assert abs(word[p])<=K-1 and word_value(word,records,identity)==form
            key=tuple(row);occurrences.setdefault(key,[]).append(k);row_cache[p,k]=row
            if previous and k==previous['k']+1:
                adjacent+=1;duplicates+=row==previous['row']
            if extend(basis,row):
                selected.append(row);rank_steps.append(dict(exponent=k,queries=k-K+1,rank=len(basis)))
            previous=item
        assert rank_steps==case['rank_steps'] and selected==case['selected_rows']
        assert (adjacent,duplicates,len(occurrences))==(case['adjacent_smooth_pairs'],case['adjacent_duplicate_rows'],case['distinct_rows'])
        assert all(len(ks)<=2*K-1 and max(ks)-min(ks)<=2*K-2 for ks in occurrences.values())
        assert case['maximum_row_multiplicity']==max(map(len,occurrences.values()))
        _,factor,_=first_useful_relation_top([records[q][1] for q in primes],selected,identity,
                       compose_forms_hnf,inverse_form,lambda f:ambiguous_factor(f,N),verify_relations=True)
        assert factor==case['factor'] and N%factor==0 and 1<factor<N
        adjacent_total+=adjacent;duplicate_total+=duplicates
    N=pool['N'];identity,records,H=models[N];primes=sorted(records);base=flags_by_N[N]
    jumps={int(p):j for p,j in pool['derived_generator_indices'].items()}
    assert all(power(records[5][1],j,identity,compose_forms_hnf)==records[p][1] for p,j in jumps.items())
    events=[];firsts={}
    for p in primes:
        K=least_power_above(p,H);flags=bytes(base[(j*jumps[p])%len(base)] for j in range(len(base)))
        firsts[p]=next(j for j in range(len(base)) if flags[(K+j)%len(base)])
        events.extend((j,p,K+j) for j in range(len(base)-K) if flags[K+j])
    basis={};selected=[];steps=[]
    for offset,p,k in sorted(events):
        if (p,k) in row_cache: row=row_cache[p,k]
        else:
            form=power(records[p][1],k,identity,compose_forms_hnf);a=form[0];parts=[]
            for q in primes:
                e=0
                while a%q==0:a//=q;e+=1
                if e:parts.append((q,e))
            assert a==1
            word,new=leading_norm_word(N,form,parts,records);assert not new
            row=[(k if q==p else 0)-word.get(q,0) for q in primes]
            assert word_value(word,records,identity)==form
        if extend(basis,row):
            selected.append(row);steps.append(dict(prime=p,exponent=k,offset=offset,rank=len(basis),
                                                  round_robin_norm_queries=len(primes)*offset+primes.index(p)+1))
            if len(basis)==len(primes):break
    assert steps==pool['round_robin_rank_steps'] and selected==pool['selected_rows']
    assert pool['one_per_target_total_norm_queries']==sum(firsts[p]+1 for p in primes)
    assert pool['pooled_round_robin_norm_queries']==steps[-1]['round_robin_norm_queries']
    _,factor,_=first_useful_relation_top([records[q][1] for q in primes],selected,identity,
                       compose_forms_hnf,inverse_form,lambda f:ambiguous_factor(f,N),verify_relations=True)
    assert factor==pool['factor'] and 1<factor<N and N%factor==0
    assert pool['total_composition_upper_bound']<=pool['composition_cap']==250000
    return dict(scope='Frozen full-pattern statistics, minimal period certificates, all 670 hit-row identities and pooled factor certificate; full producers separately replay every arithmetic phase.',
                base_periods=[len(f) for f in flags_by_N.values()],adjacent_smooth_pairs=adjacent_total,
                duplicate_rows=duplicate_total,pooled_norm_queries=pool['pooled_round_robin_norm_queries'],
                one_per_target_norm_queries=pool['one_per_target_total_norm_queries'],factor=factor,
                factoring_exponent='1/2; unchanged')


if __name__=='__main__':
    print(json.dumps(run()))
