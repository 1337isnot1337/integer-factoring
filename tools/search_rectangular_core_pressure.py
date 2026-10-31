#!/usr/bin/env python3
"""Explore all y/u/v core weights and four charges in bounded defect bases.

The regular coordinate bases preserve every valuation profile. Numerical
optimization is exploratory and requires exact replay for any strict bound.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit,logsumexp

from check_rectangular_defect_core_basis import data_for
from native_matrix_smith_f4 import MatrixOracle
from rectangular_defect_coordinates import coordinate_matrix
from search_rectangular_defect_pressure import entropy

class CoreDefectPressure:
    def __init__(self,M=8,alpha=1,beta=2,quartic_terms=None,repair_spec=None):
        self.M=M;self.data=[data_for(place,M) for place in ('infinity','zero')]
        self.oracles=[];self.metadata=[];self.row_features=[];self.column_features=[];self.spectator_signs=[]
        self.quartic_terms=quartic_terms or [(1,0)];self.queries=0;self.failures=0
        self.base_heights=4*np.arange(M)-2*M+2
        for side,data in enumerate(self.data):
            matrix,shifts,meta=coordinate_matrix(data,alpha,beta,self.quartic_terms,repair_spec);n=len(matrix)
            self.oracles.append(MatrixOracle(matrix,n,shifts));self.metadata.append(meta)
            rf=np.zeros((n,19));cf=np.zeros((n,19));rf[:,2*side]=meta['row_charges'];cf[:,2*side+1]=meta['column_charges'];o=4+7*side
            for i,state in enumerate(meta['row_core_states']):
                if side==0:
                    v,y=state&3,(state>>2)&3
                    if v:rf[i,o+v-1]=1
                else:y=state&3
                if y:rf[i,o+3+y-1]=1
            for j,state in enumerate(meta['column_core_states']):
                if side==0:u=state&3
                else:
                    v,u=state&3,(state>>2)&3
                    if v:cf[j,o+v-1]=1
                if u:cf[j,o+5+u-1]=1
            self.row_features.append(rf);self.column_features.append(cf)
            inactive=[j for j in range(2*M) if j not in data['active_indices']]
            self.spectator_signs.append(np.array([1 if side==0 or j<M else -1 for j in inactive]))

    def evaluate(self,parameters):
        theta=np.asarray(parameters,dtype=float)
        if theta.shape!=(19,) or theta[-1]<=0:raise ValueError('Eighteen gauge weights and positive log(s) are required.')
        self.queries+=1;logB=0.;dB=np.zeros(19);profiles=[]
        for side,(oracle,data,rf,cf) in enumerate(zip(self.oracles,self.data,self.row_features,self.column_features)):
            rw=rf@theta;cw=cf@theta;result=oracle.profile(rw,cw,scale=theta[-1])
            if not self.metadata[side]['exact_polynomial_matrix']:
                known=self.metadata[side]['known_below']*theta[-1]+float(min(rw))+float(min(cw))
                if max(result['orders'])>=known-1e-8:raise ArithmeticError('Weighted profile exceeds known Fock precision.')
            orders=np.array(result['orders']);logZ=logsumexp(-orders);weights=np.exp(-orders-logZ);logB+=float(logZ)
            dB-=weights@(rf[np.array(result['rows'])]+cf[np.array(result['cols'])])
            dB[-1]-=float(weights@np.diff(np.r_[0,np.array(result['base_prefix'])]))
            spectators=np.array(data['spectator_orders'],dtype=float);signs=self.spectator_signs[side]
            values=spectators*theta[-1]+signs*(theta[2*side]+theta[2*side+1]);occ=expit(-values)
            logB+=float(np.sum(np.logaddexp(0,-values)));dB[2*side]-=float(occ@signs);dB[2*side+1]-=float(occ@signs);dB[-1]-=float(occ@spectators)
            profiles.append(orders.tolist())
        etaX=theta[0]+theta[1]+theta[2]+theta[3];etaY=theta[0]+theta[1]-theta[2]-theta[3]
        gx=self.base_heights*theta[-1]+etaX;gy=self.base_heights*theta[-1]+etaY
        vx=theta[4:7]+theta[11:14];uy=theta[7:9]+theta[9:11]+theta[14:16]+theta[16:18]
        coreX=np.array([0.,gx[0]+vx[0],gx[-1]+vx[1],gx[0]+gx[-1]+vx[2]])
        coreY=np.array([0.,gy[0]+uy[0],gy[-1]+uy[1]])
        def partitions(R):
            return float(np.sum(np.logaddexp(0,-gx[1:-1]/R))+np.sum(np.logaddexp(0,-gy[1:-1]/R))+logsumexp(-coreX/R)+logsumexp(-coreY/R))
        def equation(R):return logB-R*partitions(R)
        low,high=1e-9,4.
        if equation(low)<=0:raise ArithmeticError('Nonpositive root; inspect profiles and precision.')
        while equation(high)>0:
            high*=2
            if high>1e6:raise ArithmeticError('No Gibbs root bracket.')
        for _ in range(70):
            middle=(low+high)/2
            if equation(middle)>0:low=middle
            else:high=middle
        R=(low+high)/2;px=expit(-gx/R);py=expit(-gy/R)
        cx=np.exp(-coreX/R-logsumexp(-coreX/R));cy=np.exp(-coreY/R-logsumexp(-coreY/R))
        px[0],px[-1]=cx[1]+cx[3],cx[2]+cx[3];py[0],py[-1]=cy[1],cy[2]
        H=sum(entropy([p,1-p]) for p in px[1:-1])+sum(entropy([p,1-p]) for p in py[1:-1])+entropy(cx)+entropy(cy)
        physical=np.zeros(19);physical[:2]=np.sum(px)+np.sum(py);physical[2:4]=np.sum(px)-np.sum(py);physical[-1]=self.base_heights@(px+py)
        for o in (4,11):physical[o:o+3]=cx[1:];physical[o+3:o+5]=cy[1:];physical[o+5:o+7]=cy[1:]
        gradient=(dB+physical)/(2*H)
        if not np.all(np.isfinite(gradient)):raise ArithmeticError('Nonfinite derivative.')
        return R/2,gradient,{'type_x_core':cx.tolist(),'type_y_core':cy.tolist(),'type_x_occupations':px.tolist(),'type_y_occupations':py.tolist(),'entropy':H,'log_profile_product':logB,'active_orders_times_log_s':profiles}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--M',type=int,default=8);parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2);parser.add_argument('--quartic-power',type=int,default=0);parser.add_argument('--restarts',type=int,default=12);parser.add_argument('--maxfun',type=int,default=300);parser.add_argument('--seconds',type=float,default=1800);parser.add_argument('--seed',type=int,default=91378);parser.add_argument('--height-bound',type=float,default=24);parser.add_argument('--starts',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic();oracle=CoreDefectPressure(args.M,args.alpha,args.beta,[(1,args.quartic_power)]);best=None;records=[];rng=np.random.default_rng(args.seed)
    def save():
        metadata=[{k:v for k,v in m.items() if k not in ('row_core_states','column_core_states','row_charges','column_charges')} for m in oracle.metadata]
        result={'status':'exploratory floating search; no exact rate claimed','M':args.M,'alpha':args.alpha,'beta':args.beta,'quartic_terms':oracle.quartic_terms,'queries':oracle.queries,'failures':oracle.failures,'elapsed_seconds':time.monotonic()-started,'best':best,'restarts':records,'builder_metadata':metadata,'random_seed':args.seed,'normalized_height_bound':args.height_bound,'parameter_layout':'charge_inf_row,charge_inf_col,charge_zero_row,charge_zero_col; inf v_low,v_top,v_both,y_low,y_top,u_low,u_top; zero same; log_s; weights are multiplied by log_s'}
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    def objective(x):
        nonlocal best
        if time.monotonic()-started>args.seconds:raise TimeoutError
        theta=np.r_[x[:-1]*x[-1],x[-1]]
        try:value,gradient,details=oracle.evaluate(theta)
        except ArithmeticError:
            oracle.failures+=1;return 100.,np.zeros(19)
        if best is None or value<best['rho']-1e-10:
            best={'rho':value,'parameters':theta.tolist(),'normalized_heights':x[:-1].tolist(),**details};save();print(json.dumps({'query':oracle.queries,'rho':value,'seconds':time.monotonic()-started}),flush=True)
        return value,np.r_[gradient[:-1]*x[-1],gradient[-1]+gradient[:-1]@x[:-1]]
    starts=[]
    if args.starts:
        supplied=json.loads(args.starts.read_text());supplied=supplied if isinstance(supplied,list) else [supplied['best']['parameters']]
        for p in supplied:
            if len(p)==5:p=list(p[:4])+[0.]*14+[p[-1]]
            if len(p)!=19:raise ValueError('Start must have 5 or 19 parameters.')
            starts.append(np.r_[np.array(p[:-1])/p[-1],p[-1]])
    starts.extend([np.r_[np.zeros(18),t/args.M] for t in (.3,1.,3.)])
    try:
        for x in starts:objective(x)
        for iteration in range(args.restarts):
            x=np.r_[best['normalized_heights'],best['parameters'][-1]]
            if iteration:
                x[:-1]+=rng.normal(0,.6 if iteration%3 else 2.,18);x[-1]*=math.exp(rng.normal(0,.4))
            run=minimize(objective,x,jac=True,method='L-BFGS-B',bounds=[(-args.height_bound,args.height_bound)]*18+[(.002/args.M,8/args.M)],options={'maxfun':args.maxfun,'maxiter':args.maxfun,'ftol':1e-13,'gtol':1e-7,'maxls':20})
            records.append({'restart':iteration,'success':bool(run.success),'message':str(run.message),'evaluations':int(run.nfev),'rho':float(run.fun)});save()
    except TimeoutError:records.append({'status':'time budget reached'})
    finally:save()
    print(json.dumps({'best_rho':None if best is None else best['rho'],'queries':oracle.queries,'seconds':time.monotonic()-started}),flush=True)

if __name__=='__main__':main()
