#!/usr/bin/env python3
"""Scout deliberately unbalanced bulk charges and refine distinct phase starts.

The bulk zero-temperature cancellation diamond is a guide to starts, not
a pressure certificate. All saved numerical values remain exploratory.
"""
import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from search_rectangular_core_pressure import CoreDefectPressure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--M',type=int,default=8)
    parser.add_argument('--full-alphabet-gamma',type=int,choices=(0,2,3))
    parser.add_argument('--starts',type=Path)
    parser.add_argument('--resume',type=Path)
    parser.add_argument('--include-edges',action='store_true')
    parser.add_argument('--seconds',type=float,default=3600)
    parser.add_argument('--maxfun',type=int,default=180)
    parser.add_argument('--refine',type=int,default=12)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();began=time.monotonic();deadline=began+args.seconds
    oracle=CoreDefectPressure(args.M,full_alphabet_gamma=args.full_alphabet_gamma)
    dimension=oracle.parameters;best=None;scouts=[];refinements=[];failures=[]
    core_starts=[('zero_core',np.zeros(dimension-1))]
    if args.starts:
        source=json.loads(args.starts.read_text());point=source.get('best',source)
        vector=point.get('parameters')
        if vector is not None and len(vector)==dimension:
            core_starts.append(('inherited_core',np.array(vector[:-1])/vector[-1]))
    phases=[(0.,0.),(1.,0.),(-1.,0.),(0.,1.),(0.,-1.),
            (1.5,0.),(-1.5,0.),(0.,1.5),(0.,-1.5)]
    phases.extend((sx*a,sy*b) for a,b in ((.5,.5),(.75,.75),(1.25,.25),(.25,1.25))
                  for sx in (-1,1) for sy in (-1,1))
    if args.include_edges:
        phases.extend([(2.,0.),(-2.,0.),(0.,2.),(0.,-2.),
                       (1.,1.),(1.,-1.),(-1.,1.),(-1.,-1.),
                       (1.5,.5),(1.5,-.5),(-1.5,.5),(-1.5,-.5),
                       (.5,1.5),(.5,-1.5),(-.5,1.5),(-.5,-1.5)])
    previous=None
    if args.resume:
        previous=json.loads(args.resume.read_text())
        if previous['M']!=args.M or previous.get('full_alphabet_gamma')!=args.full_alphabet_gamma:
            raise ValueError('Resume construction does not match.')
        scouts=previous['scouts'];refinements=previous.get('refinements',[]);best=previous.get('best')
        failures=previous.get('failures',[])
        vectors=previous.get('core_start_vectors')
        if vectors is None:
            vectors={}
            for record in scouts:vectors.setdefault(record['core'],record['normalized_parameters'][:-1])
        core_starts=[(name,np.array(vector,dtype=float)) for name,vector in vectors.items()]
    report={'status':'exploratory floating bulk-charge phase screen; no exact rate or exclusion',
            'M':args.M,'full_alphabet_gamma':args.full_alphabet_gamma,
            'phases':[list(p) for p in phases],'tau_starts':[3.,5.],
            'core_start_names':[name for name,_ in core_starts],
            'core_start_vectors':{name:vector.tolist() for name,vector in core_starts},
            'resume_source':str(args.resume) if args.resume else None,
            'previous_queries':previous.get('queries',0) if previous else 0}

    def save():
        report.update(best=best,scouts=scouts,refinements=refinements,failures=failures,
                      queries=oracle.queries,elapsed_seconds=time.monotonic()-began)
        temporary=args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(report,indent=2)+'\n');temporary.replace(args.output)

    def query(normalized,context):
        nonlocal best
        if time.monotonic()>deadline:raise TimeoutError
        theta=np.r_[normalized[:-1]*normalized[-1],normalized[-1]]
        try:rho,gradient,details=oracle.evaluate(theta)
        except (ArithmeticError,ValueError,OverflowError) as error:
            failures.append({'context':context,'error':repr(error)});return None
        if best is None or rho<best['rho']-1e-10:
            best={'rho':rho,'parameters':theta.tolist(),'normalized_parameters':normalized.tolist(),
                  'context':context,**details};save()
            print(json.dumps({'query':oracle.queries,'rho':rho,'context':context,'seconds':time.monotonic()-began}),flush=True)
        derivative=np.r_[gradient[:-1]*normalized[-1],gradient[-1]+gradient[:-1]@normalized[:-1]]
        return rho,derivative,details

    try:
        completed={(tuple(record['phase']),record['tau'],record['core']) for record in scouts}
        for a,b in phases:
            for tau in (3.,5.):
                for core_name,core in core_starts:
                    if ((a,b),tau,core_name) in completed:continue
                    x=np.r_[core.copy(),tau/args.M]
                    # Keep the inherited differences while replacing each sum.
                    d0=x[0]-x[1];d1=x[2]-x[3]
                    x[0]=(a*args.M+d0)/2;x[1]=(a*args.M-d0)/2
                    x[2]=(b*args.M+d1)/2;x[3]=(b*args.M-d1)/2
                    context={'phase':[a,b],'tau':tau,'core':core_name}
                    result=query(x,context)
                    if result is not None:
                        rho,_,details=result
                        scouts.append({'rho':rho,'normalized_parameters':x.tolist(),**context,
                                       'type_x_core':details['type_x_core'],'type_y_core':details['type_y_core']})
                    save()
        # Refine distinct interior phases. Exact edges lose a physical half
        # of the thermal term, so they are optional finite controls only.
        by_phase={}
        for record in sorted(scouts,key=lambda r:r['rho']):by_phase.setdefault(tuple(record['phase']),record)
        selected=[record for record in sorted(by_phase.values(),key=lambda r:r['rho'])
                  if sum(abs(a) for a in record['phase'])<2]
        selected=selected[:args.refine]
        for index,start in enumerate(selected):
            context={'refinement':index,'initial_phase':start['phase'],'initial_core':start['core']}
            def objective(x):
                result=query(x,context)
                return (100.,np.zeros(dimension)) if result is None else result[:2]
            result=minimize(objective,np.array(start['normalized_parameters']),jac=True,
                            method='L-BFGS-B',bounds=[(-8*args.M,8*args.M)]*(dimension-1)+[(.01/args.M,16/args.M)],
                            options={'maxfun':args.maxfun,'maxiter':args.maxfun,'ftol':1e-13,'gtol':1e-8,'maxls':25})
            refinements.append({'initial_phase':start['phase'],'initial_rho':start['rho'],
                                'final_rho':float(result.fun),'evaluations':int(result.nfev),
                                'success':bool(result.success),'message':str(result.message),
                                'final_normalized_parameters':result.x.tolist()})
            save()
    except TimeoutError:report['time_budget_reached']=True
    finally:save()
    print(json.dumps({'best_rho':None if best is None else best['rho'],'queries':oracle.queries,
                      'scouts':len(scouts),'refinements':len(refinements),'failures':len(failures)}),flush=True)


if __name__=='__main__':main()
