#!/usr/bin/env python3
"""Explore bounded interval defects with four proved real charge gauges.

Each profile uses a structurally proved generic rank and the native exact
F4 series engine. Real gauge selection and pressure optimization are
floating point; the output is a search record, never an exact bound.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit,logsumexp

from native_matrix_smith_f4 import MatrixOracle
from rectangular_defect_profiles import active_matrix

def entropy(p):
    p=np.asarray(p);positive=p>0
    return float(-np.sum(p[positive]*np.log(p[positive])))

class DefectPressure:
    def __init__(self,infinity_data,zero_data,alpha=1,beta=2,quartic_terms=None):
        if infinity_data['M']!=zero_data['M']:raise ValueError('Mismatched intervals.')
        self.M=infinity_data['M'];self.data=[infinity_data,zero_data]
        self.oracles=[];self.metadata=[];self.charges=[];self.spectator_signs=[]
        self.coefficients=[alpha,beta];self.quartic_terms=quartic_terms or [(1,0)]
        for side,data in enumerate(self.data):
            matrix,shifts,metadata=active_matrix(data,alpha,beta,self.quartic_terms)
            n=len(matrix);self.oracles.append(MatrixOracle(matrix,3*n//4,shifts,precision=metadata['profile_precision']));self.metadata.append(metadata)
            signs=[1 if side==0 or j<self.M else -1 for j in data['active_indices']]
            self.charges.append(np.array([sum(signs[j] for j in range(len(signs)) if mask>>j&1) for mask in range(n)]))
            inactive=[j for j in range(2*self.M) if j not in data['active_indices']]
            self.spectator_signs.append(np.array([1 if side==0 or j<self.M else -1 for j in inactive]))
        self.base_heights=4*np.arange(self.M)-2*self.M+2
        self.queries=0;self.failures=0

    def evaluate(self,parameters):
        theta=np.asarray(parameters,dtype=float)
        if theta.shape!=(5,) or theta[-1]<=0:raise ValueError('Four charge weights and positive log(s) are required.')
        self.queries+=1;logB=0.;dB=np.zeros(5);profiles=[]
        for side,(oracle,data,charges) in enumerate(zip(self.oracles,self.data,self.charges)):
            offset=2*side
            result=oracle.profile(charges*theta[offset],charges*theta[offset+1],scale=theta[4])
            if not self.metadata[side]['exact_polynomial_matrix']:
                known=self.metadata[side]['known_below']*theta[4]+float(min(charges*theta[offset]))+float(min(charges*theta[offset+1]))
                if max(result['orders'])>=known-1e-8:raise ArithmeticError('Weighted profile exceeds the known Fock precision.')
            orders=np.array(result['orders']);logZ=logsumexp(-orders);weights=np.exp(-orders-logZ);logB+=float(logZ)
            dB[offset]-=float(weights@charges[np.array(result['rows'])]);dB[offset+1]-=float(weights@charges[np.array(result['cols'])])
            bases=np.diff(np.r_[0,np.array(result['base_prefix'])]);dB[4]-=float(weights@bases)
            spectators=np.array(data['spectator_orders'],dtype=float);signs=self.spectator_signs[side]
            values=spectators*theta[4]+signs*(theta[offset]+theta[offset+1]);occ=expit(-values)
            logB+=float(np.sum(np.logaddexp(0,-values)))
            dB[offset]-=float(occ@signs);dB[offset+1]-=float(occ@signs);dB[4]-=float(occ@spectators)
            profiles.append(orders.tolist())
        gx=self.base_heights*theta[4]+theta[0]+theta[1]+theta[2]+theta[3]
        gy=self.base_heights*theta[4]+theta[0]+theta[1]-theta[2]-theta[3]
        def partitions(R):
            return float(np.sum(np.logaddexp(0,-gx/R))+np.sum(np.logaddexp(0,-gy[1:-1]/R))+logsumexp([0.,-gy[0]/R,-gy[-1]/R]))
        def equation(R):return logB-R*partitions(R)
        low,high=1e-9,4.
        if equation(low)<=0:raise ArithmeticError('Nonpositive root; inspect exact profile and precision.')
        while equation(high)>0:
            high*=2
            if high>1e6:raise ArithmeticError('No Gibbs root bracket.')
        for _ in range(70):
            middle=(low+high)/2
            if equation(middle)>0:low=middle
            else:high=middle
        R=(low+high)/2;px=expit(-gx/R);py=expit(-gy/R)
        core_log=np.array([0.,-gy[0]/R,-gy[-1]/R]);core=np.exp(core_log-logsumexp(core_log));py[0],py[-1]=core[1],core[2]
        H=sum(entropy([p,1-p]) for p in px)+sum(entropy([p,1-p]) for p in py[1:-1])+entropy(core)
        physical=np.array([np.sum(px)+np.sum(py)]*2+[np.sum(px)-np.sum(py)]*2+[self.base_heights@(px+py)])
        gradient=(dB+physical)/(2*H)
        if not np.all(np.isfinite(gradient)):raise ArithmeticError('Nonfinite derivative.')
        return R/2,gradient,{'type_x_occupations':px.tolist(),'type_y_occupations':py.tolist(),'type_y_core':core.tolist(),'entropy':H,'log_profile_product':logB,'active_orders_times_log_s':profiles}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path);parser.add_argument('--M',type=int,default=8)
    parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2)
    parser.add_argument('--quartic-power',type=int,default=0);parser.add_argument('--restarts',type=int,default=8)
    parser.add_argument('--maxfun',type=int,default=400);parser.add_argument('--seconds',type=float,default=1800)
    parser.add_argument('--seed',type=int,default=8914);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--starts',type=Path);args=parser.parse_args();started=time.monotonic()
    data=[json.loads((args.directory/f'defect-{place}-M{args.M}.json').read_text()) for place in ('infinity','zero')]
    oracle=DefectPressure(*data,args.alpha,args.beta,[(1,args.quartic_power)])
    rng=np.random.default_rng(args.seed);best=None;records=[]
    def save():
        result={'status':'exploratory floating search; no exact rate claimed','M':args.M,'alpha':args.alpha,'beta':args.beta,'quartic_terms':oracle.quartic_terms,'queries':oracle.queries,'failures':oracle.failures,'elapsed_seconds':time.monotonic()-started,'best':best,'restarts':records,'builder_metadata':oracle.metadata,'random_seed':args.seed,'height_bounds':[-12,12],'log_s_bounds':[.002/args.M,8/args.M]}
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    def objective(x):
        nonlocal best
        if time.monotonic()-started>args.seconds:raise TimeoutError
        try:value,gradient,details=oracle.evaluate(x)
        except ArithmeticError:
            oracle.failures+=1;return 100.,np.zeros(5)
        if best is None or value<best['rho']-1e-10:
            best={'rho':value,'parameters':np.asarray(x).tolist(),**details};save()
            print(json.dumps({'query':oracle.queries,'rho':value,'parameters':np.asarray(x).tolist(),'seconds':time.monotonic()-started}),flush=True)
        return value,gradient
    starts=[]
    if args.starts:
        supplied=json.loads(args.starts.read_text());starts=supplied if isinstance(supplied,list) else [supplied['best']['parameters']]
    starts.extend([[0.,0.,0.,0.,t/args.M] for t in (.2,.8,2.,4.)])
    try:
        for x in starts:objective(np.array(x))
        for iteration in range(args.restarts):
            if iteration==0:x=np.array(best['parameters'])
            else:
                x=np.array(best['parameters'])+rng.normal(0,.2 if iteration%3 else .7,5);x[-1]=max(.01/args.M,best['parameters'][-1]*math.exp(rng.normal(0,.5)))
            run=minimize(objective,x,jac=True,method='L-BFGS-B',bounds=[(-12,12)]*4+[(.002/args.M,8/args.M)],options={'maxfun':args.maxfun,'maxiter':args.maxfun,'ftol':1e-13,'gtol':1e-7,'maxls':20})
            records.append({'restart':iteration,'success':bool(run.success),'message':str(run.message),'evaluations':int(run.nfev),'rho':float(run.fun)});save()
    except TimeoutError:records.append({'status':'time budget reached'})
    finally:save()
    print(json.dumps({'best_rho':None if best is None else best['rho'],'queries':oracle.queries,'seconds':time.monotonic()-started}),flush=True)

if __name__=='__main__':main()
