#!/usr/bin/env python3
"""Exploratory F4 rectangular pressure search with independent paired types.

Requires NumPy and SciPy and the optional C++17 profile engine. Arithmetic
in the local series is exact, but all weights, pivot comparisons, Gibbs
roots, and optimizations are floating point. Every saved state is a search
start, not a proof or a global optimum. Exact endpoint and profile evidence
must be supplied separately before using a state in a factoring theorem.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize

from native_rectangular_smith_f4 import RectangularOracle


def logsumexp(values):
    top=max(values)
    return top+math.log(sum(math.exp(v-top) for v in values))


def paired_value_gradient(oracle, params, cut=1):
    """Use one native profile pair; differentiate the minimized type value."""
    x=np.asarray(params,dtype=float)
    uniform_x=[1/oracle.A]*oracle.A
    uniform_y=[1/oracle.B]*oracle.B
    _,gradient,orders=oracle.evaluate(x,uniform_x,uniform_y,cut,True)
    gx=np.zeros(oracle.A);gy=np.zeros(oracle.B)
    for side in range(2):
        for leg in range(4):
            g=gx if leg in (0,3) else gy
            start=side*oracle.stride+oracle.offsets[leg]
            g[1:]+=x[start:start+len(g)-1]
    minimum_x=float(min(gx));minimum_y=float(min(gy))
    constant=(minimum_x+minimum_y+logsumexp([-v for v in orders[:oracle.n]])
              +logsumexp([-v for v in orders[oracle.n:]]))
    def equation(rate):
        return constant-rate*(logsumexp([-(v-minimum_x)/rate for v in gx])
                              +logsumexp([-(v-minimum_y)/rate for v in gy]))
    low=1e-10;high=8.
    if equation(low)<=0:
        raise ArithmeticError('Nonpositive paired Gibbs root needs exact investigation.')
    while equation(high)>0:
        high*=2
        if high>1e8:raise ArithmeticError('Paired root not bracketed.')
    rate=brentq(equation,low,high,xtol=2e-13,rtol=2e-13)
    px=np.exp(-gx/rate-logsumexp(list(-gx/rate)))
    py=np.exp(-gy/rate-logsumexp(list(-gy/rate)))
    entropy=-sum(p*math.log(p) for p in list(px)+list(py) if p>0)
    if entropy<1e-10:raise ArithmeticError('Type entropy is too small for stable optimization.')
    derivative=np.asarray(gradient,dtype=float)
    for side in range(2):
        for leg in range(4):
            p=px if leg in (0,3) else py
            start=side*oracle.stride+oracle.offsets[leg]
            derivative[start:start+len(p)-1]+=p[1:]-1/len(p)
    derivative/=2*entropy
    return rate/2,derivative,px.tolist(),py.tolist()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed',type=Path)
    parser.add_argument('--starts',type=Path)
    parser.add_argument('--cut',type=int,choices=(0,1),default=1)
    parser.add_argument('--restarts',type=int,default=12)
    parser.add_argument('--maxfun',type=int,default=2500)
    parser.add_argument('--seconds',type=float,default=300.)
    parser.add_argument('--random-seed',type=int,default=7351)
    parser.add_argument('--height-bound',type=float,default=60.)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    seed=json.loads(args.seed.read_text());a,b=seed['dimensions'][:2]
    oracle=RectangularOracle(seed['values'],a,b,seed.get('denominator_degree',0),
                             seed.get('strict_endpoint',True))
    rng=np.random.default_rng(args.random_seed)
    starts=[]
    if args.starts:
        source=json.loads(args.starts.read_text())
        for record in source.get('records',[source]):
            if record.get('dimensions',seed['dimensions'])==seed['dimensions']:
                vector=record.get('parameters',record.get('best_parameters'))
                if vector is not None and len(vector)==oracle.parameters:starts.append(np.array(vector,dtype=float))
    starts.append(np.r_[np.zeros(oracle.parameters-1),.15])
    initial=time.monotonic();deadline=initial+args.seconds
    best=None;attempts=0;failures=0;history=[]
    coefficient_hash=hashlib.sha256(json.dumps(seed['values'],separators=(',',':')).encode()).hexdigest()
    report={'status':'exploratory floating search; no certified exponent or global optimum',
            'seed_coefficient_sha256':coefficient_hash,
            'constructor':{k:seed[k] for k in ('indices','ambient','alpha','beta','delta','repair_kind') if k in seed},
            'dimensions':seed['dimensions'],'denominator_degree':oracle.denominator_degree,
            'cut':args.cut,'random_seed':args.random_seed,'records':[]}

    def save():
        report.update(profile_queries_attempted=attempts,failed_profiles=failures,
                      elapsed_seconds=time.monotonic()-initial,restarts=history)
        if best is not None:report['records']=[best]
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')

    class Deadline(Exception):pass

    def objective(x):
        nonlocal attempts,failures,best
        if time.monotonic()>deadline:raise Deadline()
        attempts+=1
        try:value,gradient,px,py=paired_value_gradient(oracle,x,args.cut)
        except (ArithmeticError,ValueError,OverflowError):
            failures+=1
            return 100.,np.zeros(oracle.parameters)
        if best is None or value<best['rho_floating']-1e-10:
            best={'dimensions':seed['dimensions'],'parameters':list(map(float,x)),
                  'rho_floating':value,'type_x_v':px,'type_y_u':py,'cut':args.cut,
                  'query':attempts}
            save()
        return value,gradient

    try:
        for restart in range(args.restarts):
            if restart<len(starts):start=starts[restart].copy()
            elif best is not None:
                start=np.array(best['parameters'])
                magnitude=(.03,.15,.6,2.)[(restart-len(starts))%4]
                start[:-1]+=rng.normal(0,magnitude,size=oracle.parameters-1)
                start[-1]*=math.exp(rng.normal(0,.2))
            else:start=np.r_[rng.normal(0,1,oracle.parameters-1),.15]
            start[:-1]=np.clip(start[:-1],-args.height_bound,args.height_bound)
            start[-1]=min(10.,max(.002,start[-1]))
            result=minimize(objective,start,method='L-BFGS-B',jac=True,
                            bounds=[(-args.height_bound,args.height_bound)]*(oracle.parameters-1)+[(.002,10.)],
                            options={'maxfun':args.maxfun,'maxiter':args.maxfun,'ftol':1e-13,
                                     'gtol':1e-8,'maxls':35,'maxcor':20})
            history.append({'restart':restart,'optimizer_success':bool(result.success),
                            'message':str(result.message),'evaluations':int(result.nfev),
                            'final_rho':float(result.fun)})
            save()
            print(json.dumps({'restart':restart,'best':None if best is None else best['rho_floating'],
                              'queries':attempts,'seconds':time.monotonic()-initial}),flush=True)
    except Deadline:
        report['time_budget_reached']=True
    save()
    print(json.dumps({'best':None if best is None else best['rho_floating'],
                      'queries':attempts,'failed_profiles':failures,
                      'seconds':time.monotonic()-initial}),flush=True)


if __name__=='__main__':main()
