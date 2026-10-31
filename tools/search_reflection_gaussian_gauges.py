#!/usr/bin/env python3
"""Optimize all mode gauges of the compressed Gaussian in its regular cone.

Only internal Gaussian coefficients constrained to be regular are
removed. Both complete one-particle profiles and both physical-type
partition functions are charged. Floating optimization is exploratory.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.optimize import brentq,minimize,LinearConstraint
from scipy.special import expit
from reflection_dual_rail_gaussian import laurent_kernel,logical_heights
from native_matrix_smith_f4 import MatrixOracle


class GaussianGaugePressure:
    def __init__(self,interval):
        self.interval=interval;self.modes=interval//2;self.size=2*interval
        self.kernel=laurent_kernel(interval);self.cuts=[];constraints=[]
        for side,legs in enumerate(((0,1),(0,2))):
            row=[leg*self.modes+i for leg in legs for i in range(self.modes)]
            col=[leg*self.modes+i for leg in range(4) if leg not in legs for i in range(self.modes)]
            matrix=[];shifts=[]
            for i in row:
                values=[];powers=[]
                for j in col:
                    entry=self.kernel[i][j]
                    first=min(entry) if entry else 0
                    values.append((sum(1<<(p-first) for p in entry),0));powers.append(first)
                matrix.append(values);shifts.append(powers)
            precision=max(512,1<<(4*interval-1).bit_length())
            self.cuts.append((MatrixOracle(matrix,interval,shifts,precision),row,col))
            for group in (row,col):
                for a,i in enumerate(group):
                    for j in group[a+1:]:
                        if not self.kernel[i][j]:continue
                        order=min(self.kernel[i][j]) if side==0 else -max(self.kernel[i][j])
                        coefficient=np.zeros(2*self.size+1)
                        coefficient[side*self.size+i]=coefficient[side*self.size+j]=1
                        coefficient[-1]=order;constraints.append(coefficient)
        self.constraints=np.array(constraints)
        self.last=None;self.queries=0

    def evaluate(self,parameters,gradient=False):
        parameters=np.asarray(parameters,dtype=float);ell=parameters[-1]
        if ell<=0:raise ValueError('A positive profile scale is required.')
        if np.min(self.constraints@parameters)<-2e-7:raise ArithmeticError('An internal Gaussian coefficient is not regular.')
        log_b=0.;derivatives=np.zeros_like(parameters);records=[]
        for side,(oracle,row,col) in enumerate(self.cuts):
            heights=parameters[side*self.size:(side+1)*self.size]
            result=oracle.profile(heights[row],heights[col],ell,infinity=bool(side))
            orders=np.array(result['orders']);log_b+=float(np.logaddexp(0,-orders).sum())
            probabilities=expit(-orders);base=np.diff([0]+result['base_prefix'])
            for p,r,c,b in zip(probabilities,result['rows'],result['cols'],base):
                derivatives[side*self.size+row[r]]-=p
                derivatives[side*self.size+col[c]]-=p
                derivatives[-1]-=p*b
            records.append(result)
        physical=np.zeros((2,self.modes));height_array=parameters[:-1].reshape(2,4,self.modes)
        physical[0]=height_array[:,0,:].sum(axis=0)+height_array[:,3,:].sum(axis=0)
        physical[1]=height_array[:,1,:].sum(axis=0)+height_array[:,2,:].sum(axis=0)
        def equation(rate):return log_b-rate*float(np.logaddexp(0,-physical/rate).sum())
        if equation(1e-9)<=0:raise ArithmeticError('The positive Gibbs root is unresolved near zero entropy.')
        high=4.
        while equation(high)>0:high*=2
        rate=brentq(equation,1e-9,high,xtol=1e-13)
        p=expit(-physical/rate)
        entropy=float((np.logaddexp(0,-physical/rate)+physical/rate*p).sum())
        if entropy<1e-9:raise ArithmeticError('The physical entropy is too small for a stable rate.')
        for side in range(2):
            for leg in range(4):
                derivatives[side*self.size+leg*self.modes:side*self.size+(leg+1)*self.modes]+=p[0 if leg in (0,3) else 1]
        self.queries+=1
        self.last={'rho':rate/2,'entropy':entropy,'parameters':parameters.tolist(),
                   'critical_pressure':equation(math.sqrt(2)),'one_particle_profiles':records,
                   'internal_regularity_margin':float(np.min(self.constraints@parameters))}
        return (rate/2,derivatives/(2*entropy)) if gradient else rate/2

    def natural(self,ell):
        return np.array([ell*h for place in ('zero','infinity') for h in logical_heights(self.interval,place)]+[ell])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval',type=int,default=16);parser.add_argument('--restarts',type=int,default=20)
    parser.add_argument('--iterations',type=int,default=1200);parser.add_argument('--seconds',type=float,default=3600)
    parser.add_argument('--seed',type=int,default=310750);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();pressure=GaussianGaugePressure(args.interval);rng=np.random.default_rng(args.seed)
    from scipy.optimize import minimize_scalar
    scale=minimize_scalar(lambda x:pressure.evaluate(pressure.natural(math.exp(x))),bounds=(math.log(.001),math.log(2)),method='bounded')
    ell=math.exp(scale.x);natural=pressure.natural(ell);pressure.evaluate(natural);baseline=pressure.last.copy()
    # An exact-profile, finite-difference gradient control inside the cone.
    control=natural.copy();control[:-1]+=.017+rng.uniform(-.003,.003,len(control)-1);value,analytic=pressure.evaluate(control,True)
    gradient_error=0.
    for j in range(min(len(control),17)):
        delta=np.zeros_like(control);delta[j]=1e-6
        finite=(pressure.evaluate(control+delta)-pressure.evaluate(control-delta))/(2e-6)
        gradient_error=max(gradient_error,abs(finite-analytic[j]))
    assert gradient_error<2e-5,gradient_error
    start=time.monotonic();deadline=start+args.seconds;best=baseline.copy()
    result={'status':'exploratory floating search inside the proved regular Gaussian gauge cone',
            'interval':args.interval,'baseline':baseline,'runs':[],'gradient_max_error':gradient_error,
            'random_seed':args.seed,'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                for name in (Path(__file__).name,'reflection_dual_rail_gaussian.py')}}
    def save():
        result['best']=best;result['queries']=pressure.queries;result['seconds']=time.monotonic()-start
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    constraints=LinearConstraint(pressure.constraints,0,np.inf)
    for trial in range(args.restarts):
        if time.monotonic()>=deadline:break
        initial=natural.copy()
        if trial:
            initial[:-1]+=rng.normal(size=len(initial)-1)*ell*(.2 if trial<args.restarts//2 else 1.)
            margin=float(np.min(pressure.constraints@initial))
            if margin<0:initial[:-1]+=-margin/2+.0001
        def objective(parameters):
            nonlocal best
            if time.monotonic()>=deadline:raise TimeoutError('Requested search budget exhausted.')
            margin=float(np.min(pressure.constraints@parameters))
            if margin<-2e-7:
                violation=np.minimum(pressure.constraints@parameters,0)
                return 2.+float(violation@violation),2*pressure.constraints.T@violation
            try:answer=pressure.evaluate(parameters,True)
            except ArithmeticError:return 2.,np.zeros_like(parameters)
            if answer[0]<best['rho']:
                best=pressure.last.copy();save()
            return answer
        try:
            fit=minimize(objective,initial,jac=True,method='SLSQP',constraints=[constraints],
                         bounds=[(-50,50)]*(len(initial)-1)+[(.0005,2.)],
                         options={'maxiter':args.iterations,'ftol':1e-12})
            result['runs'].append({'trial':trial,'rho':float(fit.fun),'success':bool(fit.success),
                                   'message':str(fit.message),'iterations':int(fit.nit)})
        except TimeoutError:result['budget_exhausted']=True;break
        save();print(json.dumps({'trial':trial,'best_rho':best['rho'],'queries':pressure.queries,'seconds':result['seconds']}),flush=True)
    save()

if __name__=='__main__':main()
