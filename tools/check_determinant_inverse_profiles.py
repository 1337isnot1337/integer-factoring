#!/usr/bin/env python3
"""Reconstruct the registered determinant panel and its independent controls."""
from hashlib import sha256
import json
import math
from pathlib import Path

from coherent_hot_compounds_f8 import partial_profiles
from determinant_inverse_codes import code,serializable
from probe_determinant_inverse_codes import entry_controls,lse

ROOT=Path(__file__).resolve().parents[1]


def without_seconds(value):
    return {k:v for k,v in value.items() if k!='seconds'}


def main():
    d=json.loads((ROOT/'evidence/determinant-code-panel.json').read_text())
    assert d['source_commit']=='4fcd81c666c56a998e1437d3ed57e1940ebd6b8a'
    assert d['source_tree']=='b21bc8ba7d4cac31c3a56ad3de1bdcf5ff91e840'
    assert d['status']=='complete' and d['pre_status']==d['post_status']==''
    assert not d['factoring_exponent_improved']
    for path,digest in d['source_sha256'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    assert [(r['k'],r['kind']) for r in d['cases']]==[
        (k,kind) for k in (2,3) for kind in ('diagonal','antidiagonal','cauchy')]
    entries=nonbinary=orders=0
    for r in d['cases']:
        acquired=code(r['k'],r['kind']);assert serializable(acquired)==r['maps']
        expected=r['complete_profile'];assert expected['all_orders_resolved']
        assert len(r['attempts'])==1 and expected['cutoff']==127
        actual=partial_profiles(acquired['embedding'],acquired['retraction'],r['k']**2,127)
        assert serializable(without_seconds(actual))==without_seconds(expected)
        controls=entry_controls(acquired)
        assert serializable(controls)==r['independent_controls']
        entries+=controls['entries'];nonbinary+=controls['nonbinary_entries']
        for cut in expected['cuts']:
            assert cut['tail_count']==0 and len(cut['retained_orders_t'])==acquired['q']**2
            orders+=len(cut['retained_orders_t'])
        fit=r['numerical_fit'];rho=fit['rho'];scale=fit['scale']
        hs=[x['height'] for x in acquired['labels']]
        raw=sum(lse([-scale*v for v in cut['retained_orders_t']]) for cut in expected['cuts'])
        z=lse([-scale*h/rho for h in hs])
        assert abs(raw-4*rho*z)<1e-11
        assert all(abs(p-math.exp(-scale*h/rho-z))<1e-13
                   for p,h in zip(fit['physical_probabilities'],hs))
        assert rho>1/math.sqrt(2) and not fit['strictly_below_critical_numerically']
        assert fit['global_minimum_certified'] is False
    assert entries==3006 and nonbinary==179 and orders==2616
    print(json.dumps(dict(status='PASS',complete_cases=6,complete_orders=orders,
                          independent_polynomial_entries=entries,nonbinary_entries=nonbinary,
                          factoring_exponent_improved=False)))
    print('Numerical pressure-point consistency is checked; no global scale optimum or new score is claimed.')


if __name__=='__main__':main()
