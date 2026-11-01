#!/usr/bin/env python3
"""One registered bounded weighted-type diagnostic, never a proof certificate.

Native polynomial arithmetic is exact; weights, pivot comparisons and the
optimizer are floating point. A strict numerical witness requires a separate
exact rational replay and a complete construction before any score claim.
"""
import argparse
import datetime
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import platform
import random
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import minimize

from native_weighted_smith_f4 import Oracle

ROOT=Path(__file__).resolve().parents[1]
SEED=0x51554152544943
RATE=2*math.sqrt(2)
LIMIT=50000
BOUNDS=[(-6.,6.)]*24+[(.05,4.)]


class BudgetStop(Exception):pass
class WitnessFound(Exception):pass
class NativeRejected(Exception):pass


def logsumexp(values):
    top=max(values)
    return top+math.log(sum(math.exp(v-top) for v in values))


def starts():
    rng=random.Random(SEED);result=[]
    def natural(scale,a=1.,b=1.):
        return [scale*h for _ in range(4) for h in (a,b,a+b)]+[0.]*12+[scale]
    for scale in (.25,.5,1.,2.):
        result.append(dict(kind='natural occupation',params=natural(scale)))
    for a,b in ((.5,1.5),(1.5,.5),(.75,1.25),(1.25,.75)):
        result.append(dict(kind='unequal copy occupation',params=natural(.7,a,b)))
    for _ in range(8):
        scale=rng.uniform(.2,2.)
        point=natural(scale)
        point[:-1]=[max(-6.,min(6.,x+rng.uniform(-.5,.5))) for x in point[:-1]]
        result.append(dict(kind='prescribed seed perturbation',params=point))
    assert len(result)==16
    return result


def type_root(raw,weights):
    def equation(rate):return raw-rate*logsumexp([-w/rate for w in weights])
    low,high=1e-12,8.
    if equation(low)<=0:return None
    while equation(high)>0:
        high*=2
        if high>1e12:return None
    for _ in range(100):
        middle=(low+high)/2
        if equation(middle)>0:low=middle
        else:high=middle
    return (low+high)/8


def inputs():
    path=ROOT/'evidence/quartic-interaction-panel.json'
    evidence=json.loads(path.read_text())
    result=[]
    for name in ('Gaussian baseline','S149:x4+y4'):
        case=next(c for c in evidence['cases'] if c['name']==name)
        values=[]
        for entry in case['tensor_entries']:
            polynomial=0
            for x,y,z in entry:polynomial ^= 1 << x
            values.append([polynomial,0])
        assert len(values)==256
        result.append(dict(name=name,values=values))
    return path,result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    _,seeds=inputs();initial=starts()
    sources=['tools/search_quartic_weighted_type.py','tools/native_weighted_smith_f4.py',
             'tools/native_weighted_smith_f4.cpp','evidence/quartic-interaction-panel.json']
    def git(*parts):return subprocess.check_output(['git',*parts],cwd=ROOT,text=True).strip()
    data=dict(scope='Registered bounded floating search, not a pressure proof or factoring improvement.',
        started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_commit=git('rev-parse','HEAD'),
        source_tree=git('rev-parse','HEAD^{tree}'),pre_status=git('status','--porcelain'),
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        seed=SEED,starts=initial,bounds=BOUNDS,max_iterations=250,query_budget=LIMIT,
        witness_threshold=-1e-7,critical_rate=1/math.sqrt(2),input_tensors=seeds,
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        status='running',attempted_queries=0,completed_queries=0,rejected_queries=0,
        gradient_controls=[],runs=[],best=None)
    assert not data['pre_status'],'Run this source in a clean frozen checkout'
    log_path=args.output.with_suffix('.queries.jsonl.gz')
    digest=sha256();begun=time.monotonic();context={};local_best=None
    def save():
        data['seconds']=time.monotonic()-begun
        temporary=args.output.with_suffix('.checkpoint.tmp')
        temporary.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
        temporary.replace(args.output)
    save()
    unexpected=None
    with gzip.open(log_path,'wt',compresslevel=3) as log:
        def append_log(record):
            line=json.dumps(record,sort_keys=True,separators=(',',':'))+'\n'
            log.write(line);digest.update(line.encode())
        def query(point):
            nonlocal local_best
            if data['completed_queries']>=LIMIT:raise BudgetStop()
            x=list(map(float,point));data['attempted_queries']+=1
            assert len(x)==25 and all(lo-1e-9<=v<=hi+1e-9 for v,(lo,hi) in zip(x,BOUNDS))
            try:raw,gradient,spectra=oracle.raw_evaluate(x,context['cut'],True)
            except ArithmeticError as error:
                data['rejected_queries']+=1
                append_log(dict(event='native_rejected',attempt=data['attempted_queries'],
                                **context,params=x,message=str(error)))
                raise NativeRejected(str(error)) from error
            weights=[0.]+[sum(x[12*end+3*leg+a-1] for end in range(2) for leg in range(4))
                           for a in range(1,4)]
            normalizer=logsumexp([-w/RATE for w in weights])
            probabilities=[math.exp(-w/RATE-normalizer) for w in weights]
            margin=raw-RATE*normalizer
            gradient=list(gradient)
            for end in range(2):
                for leg in range(4):
                    for a in range(1,4):gradient[12*end+3*leg+a-1]+=probabilities[a]
            assert all(math.isfinite(v) for v in [margin,raw]+gradient+spectra+probabilities)
            data['completed_queries']+=1
            record=dict(query=data['completed_queries'],**context,params=x,margin=margin,
                        raw_pressure=raw,weights=weights,probabilities=probabilities,
                        spectra=spectra,gradient=gradient)
            append_log(record)
            if local_best is None or margin<local_best['margin']:local_best=record
            if data['best'] is None or margin<data['best']['margin']:
                data['best']=record
                data['best']['type_root_rho']=type_root(raw,weights)
            if margin < -1e-7:
                save()
                raise WitnessFound()
            return margin,np.asarray(gradient)
        try:
            for seed in seeds:
                oracle=Oracle(seed['values'],4,endpoint_check=True)
                for cut in (0,1):
                    context=dict(tensor=seed['name'],cut=cut,start='gradient control')
                    local_best=None
                    point=np.asarray(initial[8]['params']);cost,gradient=query(point)
                    errors=[]
                    for i in range(25):
                        plus=point.copy();minus=point.copy();plus[i]+=1e-6;minus[i]-=1e-6
                        finite=(query(plus)[0]-query(minus)[0])/2e-6
                        errors.append(abs(finite-gradient[i]))
                    control=dict(tensor=seed['name'],cut=cut,maximum_error=max(errors),errors=errors)
                    data['gradient_controls'].append(control)
                    if max(errors)>1e-4:
                        data['status']='gradient_control_failed'
                        raise NativeRejected('Critical-margin gradient failed its independent finite difference control')
                    for index,start in enumerate(initial):
                        context=dict(tensor=seed['name'],cut=cut,start=index)
                        local_best=None;before=data['completed_queries'];step=time.monotonic()
                        try:
                            result=minimize(query,np.asarray(start['params']),method='L-BFGS-B',jac=True,
                                bounds=BOUNDS,options=dict(maxiter=250,maxls=30,ftol=1e-12,gtol=1e-8))
                            run=dict(**context,kind=start['kind'],status=int(result.status),
                                     message=str(result.message),success=bool(result.success),
                                     iterations=int(result.nit),terminal_params=result.x.tolist(),
                                     terminal_margin=float(result.fun),best=local_best)
                        except NativeRejected as error:
                            run=dict(**context,kind=start['kind'],status='native_rejected',
                                     message=str(error),best=local_best)
                        run.update(completed_queries=data['completed_queries']-before,
                                   seconds=time.monotonic()-step)
                        data['runs'].append(run);save()
                        print(json.dumps(dict(**context,queries=run['completed_queries'],
                            best=None if local_best is None else local_best['margin'],
                            status=run['status'])),flush=True)
            data['status']='complete_no_witness'
        except BudgetStop:data['status']='budget_complete_no_witness'
        except WitnessFound:data['status']='numerical_witness_requires_exact_replay'
        except NativeRejected as error:
            data['status']='failed_native_or_gradient_control';data['failure']=str(error)
        except Exception as error:
            data['status']='unexpected_failure';data['failure']=repr(error);unexpected=error
        finally:
            data['last_context']=dict(context)
            data['last_context_best']=local_best
            log.flush()
    data.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                post_status=git('status','--porcelain'),query_log=str(log_path),
                uncompressed_query_sha256=digest.hexdigest(),
                compressed_query_sha256=sha256(log_path.read_bytes()).hexdigest())
    save()
    assert not data['post_status']
    print(json.dumps(dict(status=data['status'],completed=data['completed_queries'],
        rejected=data['rejected_queries'],runs=len(data['runs']),
        best=None if data['best'] is None else data['best']['margin'],seconds=data['seconds'])),flush=True)
    if unexpected is not None:raise unexpected


if __name__=='__main__':main()
