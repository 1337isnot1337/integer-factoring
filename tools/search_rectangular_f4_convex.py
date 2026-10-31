#!/usr/bin/env python3
"""Fixed-type convex gauge search with attaining-minor separation.

This uses floating-point optimization and the optional native F4 oracle.
It supplies reproducible candidates and numerical feasibility diagnostics,
not exact optimality certificates. NumPy, SciPy, and C++17 are required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize
from scipy.special import logsumexp, softmax

from native_rectangular_smith_f4 import RectangularOracle
from search_rectangular_f4_pressure import paired_value_gradient


class RectangularConvexSearch:
    def __init__(self,seed,cut=1,height_bound=80.,normalized_bound=40.):
        a,b=seed['dimensions'][:2]
        self.oracle=RectangularOracle(seed['values'],a,b,seed.get('denominator_degree',0),
                                     seed.get('strict_endpoint',True))
        self.cut=cut;self.height_bound=height_bound
        self.p=self.oracle.parameters;self.n=self.oracle.n
        self.dimension=self.p+2*self.n
        self.constraints={};self.equalities={}
        # Every nonzero entry gives an exact first-minor inequality. Including
        # them prevents an initial relaxation from hiding a very small pivot.
        o=self.oracle
        self.normalized_bound=normalized_bound
        for coordinate in range(self.p-1):
            for sign in (-1,1):
                row=np.zeros(self.dimension);row[coordinate]=sign;row[self.p-1]=normalized_bound
                self.constraints['normalized_bound',coordinate,sign]=row
        for index,pair in enumerate(seed['values']):
            bits=pair[0]|pair[1]
            if not bits:continue
            z=index;labels=[]
            for q in (a,b,b,a):labels.append(z%q);z//=q
            for side in range(2):
                row=np.zeros(self.dimension)
                for leg,label in enumerate(labels):
                    if label:row[side*o.stride+o.offsets[leg]+label-1]+=1
                valuation=-(bits.bit_length()-1)+o.denominator_degree if side==0 else (bits & -bits).bit_length()-1
                row[self.p-1]=valuation;row[self.p+side*self.n]=-1
                self.constraints['entry',side,index]=row
        for side in range(2):
            for k in range(self.n-1):
                row=np.zeros(self.dimension)
                row[self.p+side*self.n+k+1]=1
                row[self.p+side*self.n+k]=-1
                self.constraints['sorted',side,k]=row

    def separate(self,parameters,candidate=None,tolerance=2e-8):
        o=self.oracle;profile=o.prefixes(parameters,self.cut)
        additions=0;worst=0.
        for side in range(2):
            selected=self.cut if side==0 else 1-self.cut
            features=np.zeros(self.dimension);rows=[];cols=[]
            for k in range(self.n):
                i=side*self.n+k;r=profile['rows'][i];c=profile['cols'][i]
                labels=(r//o.B,r%o.B,c//o.A,c%o.A) if selected==0 else (r//o.B,c//o.A,r%o.B,c%o.A)
                for leg,a in enumerate(labels):
                    if a:features[side*o.stride+o.offsets[leg]+a-1]+=1
                features[self.p-1]=profile['base_prefix'][i]
                features[self.p+side*self.n+k]=-1
                rows.append(r);cols.append(c)
                key=('minor',side,k+1,tuple(sorted(rows)),tuple(sorted(cols)))
                slack=0. if candidate is None else float(features@candidate)
                worst=min(worst,slack)
                if k==self.n-1:self.equalities[side]=features.copy()
                elif candidate is None or slack < -tolerance:
                    if key not in self.constraints:
                        self.constraints[key]=features.copy();additions+=1
        return profile,additions,worst

    def mean_vector(self,px,py):
        o=self.oracle;linear=np.zeros(self.dimension)
        for side in range(2):
            for leg in range(4):
                p=px if leg in (0,3) else py
                start=side*o.stride+o.offsets[leg]
                linear[start:start+len(p)-1]=p[1:]
        return linear

    def solve(self,px,py,parameters,rounds=60,maxiter=250,deadline=None):
        linear=self.mean_vector(px,py)
        profile,_,_=self.separate(parameters)
        x=np.r_[parameters,profile['orders']]
        bound=self.height_bound
        lower=np.r_[np.full(self.p-1,-bound),.002,np.full(2*self.n,-10000.)]
        upper=np.r_[np.full(self.p-1,bound),10.,np.full(2*self.n,10000.)]
        x=np.clip(x,lower,upper)
        entropy=-sum(p*math.log(p) for p in list(px)+list(py) if p>0)
        initial_cost=self.oracle.evaluate(parameters,px,py,self.cut)[0]
        best={'parameters':list(map(float,parameters)),'rho_fixed_type':initial_cost/(2*entropy),
              'type_x_v':list(map(float,px)),'type_y_u':list(map(float,py))};history=[]

        def objective(z):
            gradient=linear.copy();value=float(linear@z)
            for side in range(2):
                sl=slice(self.p+side*self.n,self.p+(side+1)*self.n)
                value+=float(logsumexp(-z[sl]));gradient[sl]=-softmax(-z[sl])
            return value,gradient

        class Deadline(Exception):pass

        def callback(z):
            if deadline is not None and time.monotonic()>deadline:raise Deadline()

        for round_index in range(rounds):
            matrix=np.array(list(self.constraints.values()))
            equality=np.array(list(self.equalities.values()))
            started=time.monotonic()
            try:
                result=minimize(objective,x,method='SLSQP',jac=True,bounds=Bounds(lower,upper),
                                constraints=[LinearConstraint(matrix,0,np.inf),LinearConstraint(equality,0,0)],
                                callback=callback,options={'maxiter':maxiter,'ftol':2e-11})
            except Deadline:break
            x=result.x
            try:profile,added,worst=self.separate(x[:self.p],x)
            except ArithmeticError:
                history.append({'round':round_index,'native_profile_failed':True});break
            actual=np.r_[x[:self.p],profile['orders']]
            true_value=objective(actual)[0]
            rho=true_value/(2*entropy)
            if best is None or rho<best['rho_fixed_type']:
                best={'parameters':x[:self.p].tolist(),'rho_fixed_type':rho,
                      'type_x_v':list(map(float,px)),'type_y_u':list(map(float,py))}
            listed=float(min(matrix@x));eq=float(max(abs(equality@x)))
            entry={'round':round_index,'rho_fixed_type':rho,'relaxed_rho':float(result.fun/(2*entropy)),
                   'new_constraints':added,'separation_minimum_slack':worst,
                   'listed_minimum_slack':listed,'equality_error':eq,
                   'optimizer_success':bool(result.success),'optimizer_message':str(result.message),
                   'iterations':int(result.nit),'seconds':time.monotonic()-started,
                   'height_box_hit':bool(np.any(abs(x[:self.p-1])>=bound-1e-5))}
            history.append(entry);print(json.dumps(entry),flush=True)
            if not added and worst>=-2e-7 and listed>=-2e-7 and eq<=2e-7:break
            # Restore the actual feasible profile before the next solve.
            x=actual
        return best,history


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed',type=Path)
    parser.add_argument('--start',type=Path)
    parser.add_argument('--start-index',type=int,default=0)
    parser.add_argument('--cut',type=int,choices=(0,1),default=1)
    parser.add_argument('--iterations',type=int,default=12)
    parser.add_argument('--rounds',type=int,default=60)
    parser.add_argument('--maxiter',type=int,default=250)
    parser.add_argument('--height-bound',type=float,default=80.)
    parser.add_argument('--normalized-bound',type=float,default=40.)
    parser.add_argument('--seconds',type=float,default=900.)
    parser.add_argument('--interior',type=float,default=0.,help='Mix each Gibbs type with this uniform proportion.')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not 0<=args.interior<=1:raise ValueError('--interior must lie between zero and one.')
    seed=json.loads(args.seed.read_text());search=RectangularConvexSearch(seed,args.cut,args.height_bound,args.normalized_bound)
    o=search.oracle
    if args.start:
        source=json.loads(args.start.read_text());records=source.get('records',[source])
        records=[r for r in records if r.get('dimensions',seed['dimensions'])==seed['dimensions']]
        params=np.array(records[args.start_index]['parameters'])
    else:params=np.r_[np.zeros(o.parameters-1),.1]
    start=time.monotonic();deadline=start+args.seconds
    report={'status':'exploratory convex separation and Gibbs iteration; not an exact certificate',
            'dimensions':seed['dimensions'],'cut':args.cut,'denominator_degree':o.denominator_degree,
            'height_bound':args.height_bound,'normalized_height_bound':args.normalized_bound,
            'seed_coefficient_sha256':hashlib.sha256(json.dumps(seed['values'],separators=(',',':')).encode()).hexdigest(),
            'iterations':[],'records':[]}
    for iteration in range(args.iterations):
        _,_,px,py=paired_value_gradient(o,params,args.cut)
        px=(1-args.interior)*np.array(px)+args.interior/o.A
        py=(1-args.interior)*np.array(py)+args.interior/o.B
        best,history=search.solve(px,py,params,args.rounds,args.maxiter,deadline)
        if best is None:break
        params=np.array(best['parameters'])
        rho,_,qx,qy=paired_value_gradient(o,params,args.cut)
        record={'dimensions':seed['dimensions'],'parameters':params.tolist(),'rho_floating':rho,
                'type_x_v':qx,'type_y_u':qy,'cut':args.cut,'iteration':iteration}
        report['iterations'].append({'iteration':iteration,'rounds':history,'rho_floating':rho})
        if not report['records'] or rho<report['records'][0]['rho_floating']:report['records']=[record]
        report['elapsed_seconds']=time.monotonic()-start
        report['constraints']=len(search.constraints)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'iteration':iteration,'paired_rho':rho,'seconds':time.monotonic()-start}),flush=True)
        if time.monotonic()>deadline:break


if __name__=='__main__':main()
