#!/usr/bin/env python3
"""Controls for the abstract prefix-to-partition uniform-axis proof.

Source crossing orders are replayed from the prior independent receipt;
they are not recomputed here. Abstract vector cases are explicitly separate
from realized tensors. The exponentially tiny universal window is proved
analytically, not sampled in floating point.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
from random import Random

import numpy as np
from scipy.special import logsumexp, softmax


def check_case(r,s,e):
    n=len(e)
    e=np.asarray(e,dtype=float)
    r,s=np.asarray(r,dtype=float),np.asarray(s,dtype=float)
    lam=r+s[0]-e[0]
    e=e-e[0]
    delta=float(s[1]-s[0])
    assert np.all(np.cumsum(lam)<=np.cumsum(e)+1e-12)
    assert sum(lam[:2])<=sum(e[:2])-delta+1e-12
    checked=0
    for scale in (.1,1.,3.):
        ee,ll,dd=scale*e,scale*lam,scale*delta
        phi=float(logsumexp(-ee))
        entropy=phi+float(softmax(-ee)@ee)
        x=float(ee[1])
        assert entropy<=math.log(n)+1e-12
        assert entropy<=(n-1)*(1+x)*math.exp(-x)+1e-10
        for ratio in (Fraction(1,2),Fraction(1,3),Fraction(1,10),Fraction(1,100)):
            aa,bb=(1-ratio*ratio)/(1+ratio*ratio),2*ratio/(1+ratio*ratio)
            a,b=float(aa),float(bb)
            z=math.exp(-x/a)*dd*dd/(8*n*a*a)
            actual_defect=float(logsumexp(-ll/a)-logsumexp(-ee/a))
            assert actual_defect+1e-9>=math.log1p(z)
            weak=float(logsumexp(-scale*(s-s[0])/b))
            margin=a*float(logsumexp(-ll/a))+b*weak-phi
            lower=a*math.log1p(z)+b*math.log1p(math.exp(-dd/b))-(1-a)*entropy
            assert margin+1e-9>=lower
            checked+=1
    return checked


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    source=root/'evidence/cold-gap-bound-controls.json'
    data=json.loads(source.read_text())
    assert data['status']=='pass'
    source_checks=0
    for case in data['cases']:
        gx,gy=case['physical_energies']
        e=sorted(x+y for x in gx for y in gy)
        for strong in (0,1):
            source_checks+=check_case(case['profiles'][strong],case['profiles'][1-strong],e)
    # Exact rational Abel controls and separate unrealized abstract spectra.
    rng=Random(0xA815C01D)
    abstract_checks,abel_checks=0,0
    for _ in range(512):
        n=rng.randrange(2,17)
        e=sorted([0]+[rng.randrange(0,21) for _ in range(n-1)])
        lam=e[:]
        for _ in range(rng.randrange(1,6)):
            i=rng.randrange(n-1)
            j=rng.randrange(i+1,n)
            amount=rng.randrange(1,5)
            lam[i]-=amount
            lam[j]+=amount
            lam.sort()
        assert all(sum(lam[:k])<=sum(e[:k]) for k in range(1,n+1))
        deficit=sum(e[:2])-sum(lam[:2])
        delta=Fraction(deficit,rng.randrange(1,5))
        s=[Fraction(0),delta]+[delta+rng.randrange(0,9) for _ in range(n-2)]
        s.sort()
        w=sorted([Fraction(rng.randrange(1,100),100) for _ in range(n)],reverse=True)
        assert sum(wi*(li-ei) for wi,li,ei in zip(w,lam,e))<=0
        abel_checks+=1
        abstract_checks+=check_case(lam,s,e)
    # Exact coarse constants are polynomial inequalities, not tiny-float tests.
    constants=[]
    for n in (2,4,6,9,15,25,144,9216):
        X=1024*n*n
        assert X*X>48*n*(n-1)*(1+X)
        assert X>384*n*(n-1)
        assert X>3*n and X>2*(n-1)
        assert X*X>n
        constants.append({'n':n,'X':X,'log_universal_window':-2*X})
    names=['check_uniform_axis_window.py','check_cold_gap_bound.py','cold_gap_bound.py']
    report={'status':'pass','source_receipt_sha256':sha256(source.read_bytes()).hexdigest(),
            'source_independent_orders_reused_not_recomputed':data['independent_untruncated_orders'],
            'source_partition_gap_controls':source_checks,'abstract_vector_cases_not_tensors':512,
            'abstract_partition_gap_controls':abstract_checks,'exact_rational_Abel_controls':abel_checks,
            'coarse_constant_controls':constants,
            'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
            'scope':'Analytic proof establishes the uniform window. These finite floating '
                    'checks test ingredient formulas, not the exponentially small window. '
                    'Abstract spectra are not tensors or favorable construction candidates. '
                    'No factoring improvement or independent mathematical review.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','source_partition_gap_controls',
                     'abstract_partition_gap_controls','exact_rational_Abel_controls')}))


if __name__=='__main__':
    main()
