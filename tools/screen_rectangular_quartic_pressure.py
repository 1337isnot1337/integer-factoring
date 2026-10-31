#!/usr/bin/env python3
"""Screen proved scalar quartic numerators, then refine promising candidates."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from search_rectangular_core_pressure import CoreDefectPressure

def canonical(terms):
    coefficients={}
    for coefficient,power in terms:coefficients[power]=coefficients.get(power,0)^coefficient
    return tuple((coefficient,power) for power,coefficient in sorted(coefficients.items()) if coefficient)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--M',type=int,default=8);parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2);parser.add_argument('--maximum-power',type=int,default=24);parser.add_argument('--perturbed',action='store_true');parser.add_argument('--full-core',action='store_true');parser.add_argument('--starts',type=Path,action='append',required=True);parser.add_argument('--seconds',type=float,default=3600);parser.add_argument('--refine',type=int,default=5);parser.add_argument('--maxfun',type=int,default=150);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic();starts=[];dimension=25 if args.full_core else 19
    for filename in args.starts:
        data=json.loads(filename.read_text());points=data if isinstance(data,list) else [data['best']['parameters']]
        for p in points:
            if len(p)==5:p=list(p[:4])+[0.]*(dimension-5)+[p[-1]]
            elif len(p)==19 and args.full_core:p=list(p[:11])+[0.]*3+list(p[11:18])+[0.]*3+[p[-1]]
            if len(p)!=dimension:raise ValueError('Start has the wrong dimension.')
            starts.append(np.array(p))
    powers=[0]+[signed for value in range(1,args.maximum_power+1) for signed in (value,-value)]
    candidates=[canonical([(1,power)]) for power in powers]
    if args.perturbed:
        candidates.extend(canonical([(1,0),(coefficient,power),(coefficient,power+1)]) for power in powers for coefficient in (1,2,3))
    candidates=list(dict.fromkeys(candidates));records=[];best=None;queries=0;failures=0
    sources=['rectangular_defect_profiles.py','rectangular_defect_coordinates.py','native_fock_f4.cpp','search_rectangular_core_pressure.py',Path(__file__).name]
    source_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources}
    def save():
        result={'status':'exploratory fixed-start scalar screen and local refinement; no global exclusion or exact rate','M':args.M,'alpha':args.alpha,'beta':args.beta,'full_core':args.full_core,'maximum_power':args.maximum_power,'perturbed':args.perturbed,'candidate_count':len(candidates),'completed_candidates':len(records),'queries':queries,'failures':failures,'elapsed_seconds':time.monotonic()-started,'best':best,'records':records,'starts':[v.tolist() for v in starts],'source_sha256':source_hashes}
        temporary=args.output.with_suffix(args.output.suffix+'.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    def query(oracle,terms,theta):
        nonlocal best,queries,failures
        if time.monotonic()-started>args.seconds:raise TimeoutError
        queries+=1
        try:value,gradient,details=oracle.evaluate(theta)
        except ArithmeticError:
            failures+=1;return 100.,np.zeros(dimension),None
        if best is None or value<best['rho']-1e-10:
            best={'rho':value,'quartic_terms':terms,'parameters':np.asarray(theta).tolist(),**details};save()
            print(json.dumps({'queries':queries,'rho':value,'quartic_terms':terms,'seconds':time.monotonic()-started}),flush=True)
        return value,gradient,details
    try:
        for index,terms in enumerate(candidates):
            if time.monotonic()-started>args.seconds:raise TimeoutError
            try:oracle=CoreDefectPressure(args.M,args.alpha,args.beta,list(terms),full_core=args.full_core)
            except ArithmeticError as error:
                failures+=1;records.append({'quartic_terms':terms,'status':'builder unresolved','error':str(error)});save();continue
            choices=[]
            for theta in starts:
                value,_,details=query(oracle,terms,theta)
                if details is not None:choices.append((value,theta.tolist()))
            if choices:
                value,theta=min(choices,key=lambda v:v[0]);records.append({'quartic_terms':terms,'rho':value,'parameters':theta,'status':'screened'})
            else:records.append({'quartic_terms':terms,'status':'all profiles unresolved'})
            save();print(json.dumps({'screened':index+1,'total':len(candidates),'best_rho':None if best is None else best['rho'],'seconds':time.monotonic()-started}),flush=True)
        selected=sorted([r for r in records if 'rho' in r],key=lambda r:r['rho'])[:args.refine]
        for record in selected:
            terms=record['quartic_terms'];oracle=CoreDefectPressure(args.M,args.alpha,args.beta,list(terms),full_core=args.full_core)
            def objective(x):
                theta=np.r_[x[:-1]*x[-1],x[-1]];value,gradient,_=query(oracle,terms,theta)
                return value,np.r_[gradient[:-1]*x[-1],gradient[-1]+gradient[:-1]@x[:-1]]
            p=np.array(record['parameters']);x=np.r_[p[:-1]/p[-1],p[-1]]
            opt=minimize(objective,x,jac=True,method='L-BFGS-B',bounds=[(-4*args.M,4*args.M)]*(dimension-1)+[(.002/args.M,8/args.M)],options={'maxfun':args.maxfun,'maxiter':args.maxfun,'ftol':1e-13,'gtol':1e-7,'maxls':20})
            record['refinement']={'rho':float(opt.fun),'evaluations':int(opt.nfev),'success':bool(opt.success),'message':str(opt.message)};save()
    except TimeoutError:records.append({'status':'time budget reached'})
    finally:save()

if __name__=='__main__':main()
