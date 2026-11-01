#!/usr/bin/env python3
"""Bounded numerical diagnostic of necessary binary profiles, never tensors."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar


def lower_profile(D,X,Y,u):
    C=4*(D+X+Y)
    v=min((C-u)/2,3*D+2*(X+Y)-u)
    return [0.,float(u),float(v),float(C-u-v)]


def evaluate(point):
    D,X,Y,y=map(float,point)
    tangent=math.exp(y)
    a,b=1/math.sqrt(1+tangent*tangent),tangent/math.sqrt(1+tangent*tangent)
    C,K=4*(D+X+Y),2*(D+min(X,Y))
    U=min(C/3,(3*D+2*(X+Y))/2,2*D+2*min(X,Y))
    lo,hi=max(0.,K-U),min(U,K)
    assert lo<=hi+1e-10
    def term(profile,weight):
        return weight*math.log1p(sum(math.exp(-v/weight) for v in profile[1:]))
    def value(u):
        return term(lower_profile(D,X,Y,u),a)+term(lower_profile(D,X,Y,K-u),b)
    candidates=[(value(lo),lo),(value(hi),hi)]
    if hi-lo>1e-11:
        opt=minimize_scalar(value,bounds=(lo,hi),method='bounded',options={'xatol':1e-11})
        candidates.append((float(opt.fun),float(opt.x)))
    result,u=min(candidates)
    result+=D-math.log1p(math.exp(-2*X))-math.log1p(math.exp(-2*Y))
    r,s=lower_profile(D,X,Y,u),lower_profile(D,X,Y,K-u)
    for p in (r,s):
        assert min(np.diff(p))>=-1e-8
        assert p[1]<=2*D+2*min(X,Y)+1e-8
        assert p[1]+p[2]<=3*D+2*(X+Y)+1e-8
        assert abs(sum(p)-C)<1e-8
    return result,{'parameters':list(map(float,point)),'weights':[a,b],
                   'zero_infinity_normalized_profiles':[r,s],
                   'physical_ground_slack':D,'physical_excitation_energies':[2*X,2*Y],
                   'allocation_interval':[lo,hi],'allocated_second_order':u,
                   'margin':result,'scope':'Necessary real profiles only; not a tensor or a bound.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    starts=[list(base)+[angle] for base in ((0,0,0),(1,1,1),(1,2,1),(4,1,1)) for angle in (-2,0,2)]
    report={'status':'running_incomplete','registered_starts':starts,'attempts':[],
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Necessary profile relaxation, no realized tensor, global optimum or factoring improvement.'}
    begun=time.monotonic()
    best=None
    for initial in starts:
        if time.monotonic()-begun>=180:
            break
        started=time.monotonic()
        queries=0
        trial_best,last=None,None
        def objective(point):
            nonlocal queries,trial_best,last,best
            if time.monotonic()-started>=15 or time.monotonic()-begun>=180:
                raise TimeoutError('registered time cap between complete relaxation queries')
            value,details=evaluate(point)
            queries+=1
            last=details
            if trial_best is None or value<trial_best['margin']:
                trial_best=details
            if best is None or value<best['margin']:
                best=details
            return value
        result=None
        try:
            result=minimize(objective,initial,method='SLSQP',bounds=[(0,16)]*3+[(-4,4)],
                            options={'maxiter':500,'ftol':1e-12})
            status='optimizer_returned'
            message=str(result.message)
        except TimeoutError as error:
            status='budget_reached_partial'
            message=str(error)
        record={'initial_parameters':initial,'queries':queries,'best':trial_best,'last':last,
                'status':status,'message':message,'optimizer_success':bool(result.success) if result is not None else False,
                'seconds':time.monotonic()-started}
        report['attempts'].append(record)
        report['best']=best
        report['seconds']=time.monotonic()-begun
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'initial':initial,'queries':queries,'status':status,'message':message,
                          'best_margin':trial_best['margin'] if trial_best else None}),flush=True)
    report['status']='complete_registered_attempts' if len(report['attempts'])==len(starts) else 'budget_reached_partial'
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
