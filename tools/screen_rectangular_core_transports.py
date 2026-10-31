#!/usr/bin/env python3
"""Search coupled constant core maps with a certified unchanged endpoint.

Every candidate is an exact admissible constructor. Only the requested finite
fixed starts and local refinements are searched; no exclusion is inferred.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from check_rectangular_identity_transports import IDENTITY,conjugation_spec,typed_conjugation_spec
from check_rectangular_refined_repairs import height
from rectangular_coordinate_cache import cached_coordinates
from rectangular_coordinate_transforms import transform_coordinate_matrix
from search_rectangular_core_pressure import CoreDefectPressure

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--M',type=int,default=4);parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2)
    parser.add_argument('--target',choices=('partial','identity'),default='partial')
    parser.add_argument('--starts',type=Path,action='append',required=True);parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--candidates',type=Path);parser.add_argument('--equal-height-only',action='store_true')
    parser.add_argument('--seconds',type=float,default=3600);parser.add_argument('--refine',type=int,default=6);parser.add_argument('--maxfun',type=int,default=160)
    parser.add_argument('--height-bound',type=float,default=32);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic();starts=[]
    for filename in args.starts:
        raw=json.loads(filename.read_text());points=raw if isinstance(raw,list) else [raw['best']['parameters']]
        for p in points:
            if len(p)==19:p=list(p[:11])+[0.]*3+list(p[11:18])+[0.]*3+[p[-1]]
            if len(p)!=25:raise ValueError('A full core start has25 parameters.')
            starts.append(np.array(p))
    if args.candidates:
        candidates=json.loads(args.candidates.read_text())
    else:
        directions=[(r,c) for r in range(12) for c in range(12) if r!=c and (not args.equal_height_only or height(r)==height(c))]
        directions.sort(key=lambda rc:(height(rc[0])!=height(rc[1]),abs(height(rc[0])-height(rc[1])),rc))
        candidates=[{'name':'identity','V':[row[:] for row in IDENTITY]}]
        for row,column in directions:
            for coefficient in (1,2,3):
                V=[line[:] for line in IDENTITY];V[row][column]=coefficient
                candidates.append({'name':f'E_{row}_{column}_{coefficient}','row':row,'column':column,'coefficient':coefficient,'V':V})
    repair='partial' if args.target=='partial' else 'inverse'
    data,base,cache=cached_coordinates(args.cache,args.M,args.alpha,args.beta,full_core=True,repair_kind=repair)
    records=[];best=None;queries=0;failures=0
    names=['screen_rectangular_core_transports.py','rectangular_coordinate_cache.py','rectangular_coordinate_transforms.py','check_rectangular_identity_transports.py','search_rectangular_core_pressure.py']
    sources={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    def save():
        result={'status':'exploratory finite fixed-start transport screen and local refinement; no global exclusion or exact rate','M':args.M,'alpha':args.alpha,'beta':args.beta,'target':args.target,'full_core':True,'candidate_count':len(candidates),'completed_candidates':len(records),'queries':queries,'failures':failures,'elapsed_seconds':time.monotonic()-started,'best':best,'records':records,'starts':[v.tolist() for v in starts],'cache':cache,'source_sha256':sources}
        temporary=args.output.with_suffix(args.output.suffix+'.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    def oracle_for(candidate):
        V=candidate['V']
        spec=typed_conjugation_spec(args.alpha,args.beta,V) if args.target=='partial' else conjugation_spec(args.alpha,args.beta,V)
        compiled=[transform_coordinate_matrix(matrix,d,V,spec[1]) for matrix,d in zip(base,data)]
        return CoreDefectPressure(args.M,args.alpha,args.beta,full_core=True,compiled=compiled,source_data=data),spec
    def query(oracle,spec,candidate,theta):
        nonlocal best,queries,failures
        if time.monotonic()-started>args.seconds:raise TimeoutError
        queries+=1
        try:rho,gradient,detail=oracle.evaluate(theta)
        except ArithmeticError:
            failures+=1;return 100.,np.zeros(25),None
        if best is None or rho<best['rho']-1e-10:
            best={'rho':rho,'parameters':theta.tolist(),'candidate':candidate,'transport_spec':spec,**detail};save()
            print(json.dumps({'query':queries,'candidate':candidate['name'],'rho':rho,'seconds':time.monotonic()-started}),flush=True)
        return rho,gradient,detail
    try:
        for index,candidate in enumerate(candidates):
            if time.monotonic()-started>args.seconds:raise TimeoutError
            try:oracle,spec=oracle_for(candidate)
            except ArithmeticError as error:
                failures+=1;records.append({'candidate':candidate,'status':'unresolved coordinate window','error':str(error)});save();continue
            values=[]
            for theta in starts:
                rho,_,detail=query(oracle,spec,candidate,theta)
                if detail is not None:values.append((rho,theta.tolist()))
            if values:
                rho,theta=min(values,key=lambda v:v[0]);records.append({'candidate':candidate,'rho':rho,'parameters':theta,'status':'screened'})
            else:records.append({'candidate':candidate,'status':'all queries unresolved'})
            save();print(json.dumps({'screened':index+1,'total':len(candidates),'best_rho':None if best is None else best['rho'],'seconds':time.monotonic()-started}),flush=True)
        selected=sorted([record for record in records if 'rho' in record],key=lambda record:record['rho'])[:args.refine]
        for record in selected:
            candidate=record['candidate'];oracle,spec=oracle_for(candidate)
            def objective(x):
                theta=np.r_[x[:-1]*x[-1],x[-1]];rho,gradient,_=query(oracle,spec,candidate,theta)
                return rho,np.r_[gradient[:-1]*x[-1],gradient[-1]+gradient[:-1]@x[:-1]]
            p=np.array(record['parameters']);x=np.r_[p[:-1]/p[-1],p[-1]]
            opt=minimize(objective,x,jac=True,method='L-BFGS-B',bounds=[(-args.height_bound,args.height_bound)]*24+[(.002/args.M,8/args.M)],options={'maxfun':args.maxfun,'maxiter':args.maxfun,'ftol':1e-13,'gtol':1e-7,'maxls':20})
            record['refinement']={'rho':float(opt.fun),'evaluations':int(opt.nfev),'success':bool(opt.success),'message':str(opt.message)};save()
    except TimeoutError:records.append({'status':'time budget reached'})
    finally:save()

if __name__=='__main__':main()
