#!/usr/bin/env python3
"""Fixed Path C reproduction; see research/clustered-hits-and-pooled-rank.md."""
from collections import Counter
from hashlib import sha256
from math import gcd,isqrt,prod
from pathlib import Path
import json,lzma,sys,time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from prime_form_factor_base import build_prime_form_factor_base
from leading_power_smooth_relations import compose_with_prime_form
from binary_quadratic_forms import reduce_form
from quadratic_ideal_hnf import compose_forms_hnf
from arithmetic_two_group import power
from weighted_norm_relations import least_power_above
from check_cofactor_orbit_moments import cyclic_window_counts,ceiling_root,sufficient_moment_window,maximum_zero_run

OUT=ROOT/'evidence/clustered-power-relations';OUT.mkdir(exist_ok=True)
LIMIT=250000
used=0
started=time.monotonic()


def summarize(flags,N,p,K,kind):
    h,s=len(flags),sum(flags)
    assert 0<s<h
    D=(h+s-1)//s
    b=N.bit_length();q0=(b+1).bit_length()
    moments=[]
    for q in (q0,q0*q0):
        for C in (1,2):
            W=sufficient_moment_window(2*N,D,q,C)
            hist=Counter(cyclic_window_counts(flags,W))
            numerator=sum(number*(s*W-h*a)**(2*q) for a,number in hist.items() if h*a<s*W)
            unit=h**(q+1)*(q*W*s)**q
            required=max(1,ceiling_root(max(1,(numerator+unit-1)//unit),q))
            holds=numerator<=C**q*unit
            if holds: assert min(hist)>=1
            moments.append(dict(q=q,constant=C,window=W,density_denominator=D,
                                period_bound=2*N,window_covers_period=W>=h,minimum_hits=min(hist),
                                lower_tail_scaled_sum=numerator,constant_one_scaled_bound=unit,
                                required_integer_constant=required,holds=holds))
    tails=[]
    for mean in (1,2,4,8):
        W=(mean*h+s-1)//s
        hist=Counter(cyclic_window_counts(flags,W))
        tails.append(dict(target_mean=mean,window=W,mean_numerator=W*s,mean_denominator=h,
                          empty_phases=hist.get(0,0),phases=h))
    first=next(j for j in range(h) if flags[(K+j)%h])
    pair=[sum(flags[j]*flags[(j+lag)%h] for j in range(h)) for lag in range(17)]
    return dict(kind=kind,period=h,good=s,density_denominator=D,first_from_K=first,
                maximum_zero_run=maximum_zero_run(flags),lag_pair_counts=pair,
                tails=tails,moments=moments)


def base_cycle(N,Y,p0,targets,expected):
    global used
    delta=-N if N%4==3 else -4*N
    identity=(1,delta%2,(delta%2-delta)//4)
    height=isqrt((-delta)//3)
    base=build_prime_form_factor_base(delta,Y)
    records={r.prime:(r.root,tuple(r.positive)) for r in base['records']}
    assert N%4==3 and p0 in records
    P=prod(records)
    current=identity;forms=[];index={}
    while True:
        if used>=LIMIT: raise RuntimeError('Predeclared arithmetic composition cap reached')
        assert current not in index
        index[current]=len(forms);forms.append(current)
        current=compose_with_prime_form(current,p0,*records[p0]);used+=1
        if current==identity: break
    h=len(forms)
    assert h==expected and h<=2*N
    assert index[reduce_form(records[p0][1])]==1
    def good(a):
        g=gcd(a,N)
        return int(1<g<N or pow(P,a.bit_length(),a)==0)
    flags=bytes(good(f[0]) for f in forms)
    assert flags[0] and all(flags[j]==flags[-j] for j in range(h))
    archive=OUT/f'N{N}-Y{Y}-g{p0}.flags.xz';archive.write_bytes(lzma.compress(flags))
    result=[]
    for p in targets:
        jump=index[reduce_form(records[p][1])]
        period=h//gcd(h,jump)
        local=bytes(flags[(j*jump)%h] for j in range(period))
        entries=bytes(local[j]*(1-local[j-1]) for j in range(period))
        assert 0<sum(local)<period and sum(entries)>0
        assert all(local[j]==local[-j] for j in range(period))
        K=least_power_above(p,height)
        raw=summarize(local,N,p,K,'smooth_or_gcd')
        entry=summarize(entries,N,p,K,'cluster_entry')
        samples=[]
        phases=sorted({0,1,2,K%period,(K+raw['first_from_K'])%period,period//2,period-1})
        for phase in phases:
            direct=power(records[p][1],phase,identity,compose_forms_hnf)
            assert direct==forms[(phase*jump)%h] and good(direct[0])==local[phase]
            samples.append(dict(phase=phase,form=list(direct),good=local[phase]))
        result.append(dict(N=N,input_bits=N.bit_length(),cutoff=Y,prime=p,K=K,
                           base_prime=p0,base_period=h,discrete_index=jump,period=period,
                           smooth_classes=sum(local),clusters=sum(entries),
                           mean_cluster_numerator=sum(local),mean_cluster_denominator=sum(entries),
                           maximum_hit_run=maximum_zero_run(bytes(1-x for x in local)),
                           raw=raw,entry=entry,direct_arithmetic_samples=samples))
        print('case',N,Y,p,'period',period,'smooth',sum(local),'clusters',sum(entries),
              'raw_C',[x['required_integer_constant'] for x in raw['moments']],
              'entry_C',[x['required_integer_constant'] for x in entry['moments']],flush=True)
    certificate=dict(N=N,cutoff=Y,base_prime=p0,period=h,base_form=list(records[p0][1]),
                     all_classes_distinct=True,returned_identity=True,flags_file=str(archive.relative_to(ROOT)),
                     flags_sha256=sha256(flags).hexdigest(),archive_sha256=sha256(archive.read_bytes()).hexdigest())
    return certificate,result


certificates=[];cases=[]
for args in [(450296711,23,2,[2],20456),(1000036000099,19,5,[5,7],198850)]:
    cert,rows=base_cycle(*args);certificates.append(cert);cases.extend(rows)
names=['tools/pathc_cluster_tail/generate_orbits.py','tools/pathc_cluster_tail/audit_rows.py','tools/pathc_cluster_tail/audit_pool.py','tools/prime_form_factor_base.py','tools/leading_power_smooth_relations.py',
       'tools/quadratic_ideal_hnf.py','tools/check_cofactor_orbit_moments.py']
receipt=dict(scope='Predeclared two-cycle full-smooth-or-proper-gcd experiment; generator-seven derived in the exact generator-five cyclic subgroup. Entry lower-tail moments are finite diagnostics, not uniform density or window bounds.',
             arithmetic_composition_cap=LIMIT,base_walk_compositions=used,seconds=time.monotonic()-started,
             certificates=certificates,cases=cases,
             source_sha256={n:sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
             producer_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),factoring_exponent='1/2; unchanged')
(OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('complete',used,'seconds',receipt['seconds'],flush=True)
