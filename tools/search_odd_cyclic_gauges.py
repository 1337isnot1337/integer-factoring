#!/usr/bin/env python3
"""Explore physical gauges of odd-cyclic matrix crossings over F4.

Finite searches only: exact seed endpoints and local field arithmetic, but
floating weights and optimization. No numerical result is a proof.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize

from native_rectangular_smith_f4 import RectangularOracle
from odd_cyclic_crossing import tensor_entry
from search_rectangular_f4_pressure import paired_value_gradient

MUL = ((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))


def seed(q, corrections):
    degree = max([1]+[max(a,b) for a,b,_ in corrections])
    matrices = [[[0]*q for _ in range(q)] for _ in range(degree+1)]
    for i in range(q):
        matrices[0][i][(-i)%q] = 1
        matrices[1][i][i] ^= 1
        matrices[1][i][(-i)%q] ^= 1
    for a,b,matrix in corrections:
        for i in range(q):
            for j in range(q):
                matrices[a][i][j] ^= matrix[i][j]
                matrices[b][i][j] ^= matrix[i][j]
    values = []
    for v in range(q):
        for u in range(q):
            for y in range(q):
                for x in range(q):
                    value = [tensor_entry(matrix,x,y,u,v) for matrix in matrices]
                    total = 0
                    for coefficient in value:
                        total ^= coefficient
                    assert total == int(x == v and y == u)
                    values.append([sum((a&1)<<d for d,a in enumerate(value)),
                                   sum((a>>1)<<d for d,a in enumerate(value))])
    return {'dimensions':[q]*4,'strict_endpoint':True,'values':values,
            'matrix_coefficients':matrices,'degree':degree}


def catalogue(q, random_count, random_seed):
    yield {'kind':'reflection_pencil'}, seed(q,[])
    # Simultaneous nonzero dilation leaves the reflection reference fixed.
    positions = [(0,0),(0,1)]+[(1,j) for j in range(q)]
    for powers in ((1,2),(1,3),(2,3),(1,4)):
        for i,j in positions:
            for coefficient in (1,2,3):
                matrix = [[0]*q for _ in range(q)]
                matrix[i][j] = coefficient
                description = {'kind':'single_entry','powers':powers,'row':i,
                               'column':j,'coefficient':coefficient}
                yield description, seed(q,[(powers[0],powers[1],matrix)])
    rng = random.Random(random_seed)
    for trial in range(random_count):
        corrections = []
        for _ in range(1+trial%3):
            left = [rng.randrange(4) for _ in range(q)]
            right = [rng.randrange(4) for _ in range(q)]
            matrix = [[MUL[x][y] for y in right]for x in left]
            a,b = sorted(rng.sample(range(1,5),2))
            corrections.append((a,b,matrix))
        yield {'kind':'sampled_rank_one_corrections','trial':trial}, seed(q,corrections)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alphabet',type=int,default=5)
    parser.add_argument('--random-count',type=int,default=32)
    parser.add_argument('--random-seed',type=int,default=311026)
    parser.add_argument('--restarts',type=int,default=6)
    parser.add_argument('--maxfun',type=int,default=400)
    parser.add_argument('--seconds',type=float,default=1800)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.alphabet<3 or args.alphabet%2==0:
        raise ValueError('physical alphabet must be odd and at least three')
    started=time.monotonic();deadline=started+args.seconds
    rng=np.random.default_rng(args.random_seed)
    report={'scope':'finite floating physical-gauge search with exact F4 endpoints; no exhaustive coefficient or continuous-parameter claim',
            'arguments':vars(args)|{'output':str(args.output)},'seeds':[],
            'best_records':[],'completed_catalogue':False}
    previous=[];queries=failures=0
    failure_reasons=Counter();failure_examples=[]
    seed_best=None

    def save():
        report.update(elapsed_seconds=time.monotonic()-started,
                      profile_queries=queries,failed_queries=failures,
                      failure_reasons=dict(failure_reasons),
                      failure_examples=failure_examples)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        temporary=args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(report,indent=2)+'\n');temporary.replace(args.output)

    class Deadline(Exception):pass
    try:
        for seed_index,(description,data) in enumerate(catalogue(args.alphabet,args.random_count,args.random_seed)):
            if time.monotonic()>deadline:raise Deadline()
            oracle=RectangularOracle(data['values'],args.alphabet,args.alphabet)
            seed_best=None;seed_queries=queries;histories=[]
            for cut in (0,1):
                def objective(parameters):
                    nonlocal queries,failures,seed_best
                    if time.monotonic()>deadline:raise Deadline()
                    queries+=1
                    try:rho,gradient,px,py=paired_value_gradient(oracle,parameters,cut)
                    except (ArithmeticError,ValueError,OverflowError) as error:
                        reason=type(error).__name__+': '+str(error)
                        failure_reasons[reason]+=1
                        if len(failure_examples)<8 and not any(e['reason']==reason for e in failure_examples):
                            failure_examples.append({'reason':reason,'seed':seed_index,
                                                     'cut':cut,'parameters':parameters.tolist()})
                        failures+=1;return 100.,np.zeros(oracle.parameters)
                    if seed_best is None or rho<seed_best['rho_floating']-1e-10:
                        seed_best={'rho_floating':float(rho),'parameters':parameters.tolist(),
                                   'type_x_v':px,'type_y_u':py,'cut':cut,'query':queries}
                    return rho,gradient
                for restart in range(args.restarts):
                    if restart==0:
                        parameters=np.r_[np.zeros(oracle.parameters-1),.4]
                    elif restart==1 and previous:
                        parameters=np.array(previous[0]['parameters'])
                    else:
                        parameters=np.array(seed_best['parameters'])
                        parameters[:-1]+=rng.normal(0,(.08,.3,1.2,3.)[restart%4],oracle.parameters-1)
                        parameters[-1]*=math.exp(rng.normal(0,.2))
                    parameters[:-1]=np.clip(parameters[:-1],-30,30)
                    parameters[-1]=min(8,max(.002,parameters[-1]))
                    result=minimize(objective,parameters,method='L-BFGS-B',jac=True,
                                    bounds=[(-30,30)]*(oracle.parameters-1)+[(.002,8)],
                                    options={'maxfun':args.maxfun,'maxiter':args.maxfun,
                                             'ftol':1e-12,'gtol':1e-8,'maxls':30})
                    histories.append({'cut':cut,'restart':restart,'nfev':int(result.nfev),
                                      'success':bool(result.success),'message':str(result.message)})
            if seed_best is None:raise AssertionError('no valid pressure for an endpoint seed')
            previous=[seed_best]
            item={'index':seed_index,'description':description,'degree':data['degree'],
                  'coefficient_sha256':hashlib.sha256(json.dumps(data['values'],separators=(',',':')).encode()).hexdigest(),
                  'queries':queries-seed_queries,'optimizer_runs':histories,**seed_best}
            report['seeds'].append(item)
            kept={**item,'matrix_coefficients':data['matrix_coefficients'],'values':data['values'],
                  'dimensions':data['dimensions']}
            _,_,orders=oracle.evaluate(seed_best['parameters'],seed_best['type_x_v'],
                                       seed_best['type_y_u'],seed_best['cut'],True)
            kept['scaled_profiles']=[orders[:oracle.n],orders[oracle.n:]]
            report['best_records']=sorted(report['best_records']+[kept],key=lambda r:r['rho_floating'])[:8]
            save()
            print(json.dumps({'seed':seed_index,'rho':seed_best['rho_floating'],
                              'best':report['best_records'][0]['rho_floating'],'queries':queries,
                              'seconds':time.monotonic()-started}),flush=True)
        report['completed_catalogue']=True
    except Deadline:
        report['time_budget_reached']=True
        if seed_best is not None:
            report['interrupted_seed']={'index':seed_index,'description':description,
                                        'matrix_coefficients':data['matrix_coefficients'],
                                        'values':data['values'],**seed_best}
    save()


if __name__=='__main__':main()
