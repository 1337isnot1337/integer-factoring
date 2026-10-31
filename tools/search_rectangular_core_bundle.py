#!/usr/bin/env python3
"""Proximal cutting-plane search at fixed types, alternating with Gibbs types.

Every cut is a subgradient of the convex log profile product. Held-out
convexity inequalities are checked on every query. Numerical QP and Gibbs
steps are exploratory; no numerical stopping condition certifies a bound.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,minimize
from search_rectangular_core_pressure import CoreDefectPressure

def physical_vector(detail,M):
    px=np.array(detail['type_x_occupations']);py=np.array(detail['type_y_occupations'])
    cx=np.array(detail['type_x_core']);cy=np.array(detail['type_y_core']);stride=detail.get('core_stride',7);physical=np.zeros(2*stride+5)
    physical[:2]=sum(px)+sum(py);physical[2:4]=sum(px)-sum(py)
    physical[-1]=(4*np.arange(M)-2*M+2)@(px+py)
    for o in (4,4+stride):
        physical[o:o+3]=cx[1:];physical[o+3:o+5]=cy[1:];physical[o+5:o+7]=cy[1:]
        if stride==10:physical[o+7:o+10]=cx[1:]
    return physical

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--M',type=int,default=8);parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2);parser.add_argument('--quartic-power',type=int,default=0);parser.add_argument('--iterations',type=int,default=400);parser.add_argument('--restarts',type=int,default=4);parser.add_argument('--seconds',type=float,default=1800);parser.add_argument('--seed',type=int,default=11203);parser.add_argument('--height-bound',type=float,default=24);parser.add_argument('--starts',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--full-core',action='store_true')
    args=parser.parse_args();start=time.monotonic();oracle=CoreDefectPressure(args.M,args.alpha,args.beta,[(1,args.quartic_power)],full_core=args.full_core);dimension=oracle.parameters
    rng=np.random.default_rng(args.seed);cuts=[];best=None;records=[];minimum_slack=0.
    def save():
        result={'status':'exploratory proximal bundle search; no exact rate or global optimum claimed','M':args.M,'alpha':args.alpha,'beta':args.beta,'quartic_power':args.quartic_power,'full_core':args.full_core,'queries':oracle.queries,'failures':oracle.failures,'elapsed_seconds':time.monotonic()-start,'best':best,'iterations':records,'cuts':len(cuts),'minimum_heldout_convexity_slack':minimum_slack,'normalized_height_bound':args.height_bound,'random_seed':args.seed}
        temporary=args.output.with_suffix(args.output.suffix+'.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    def query(theta):
        nonlocal best,minimum_slack
        if time.monotonic()-start>args.seconds:raise TimeoutError
        rho,gradient,detail=oracle.evaluate(theta);physical=physical_vector(detail,args.M)
        f=detail['log_profile_product'];g=gradient*(2*detail['entropy'])-physical
        if cuts:
            slacks=[f-(old_f+old_g@(theta-old_theta)) for old_theta,old_f,old_g in cuts]
            slacks.extend(old_f-(f+g@(old_theta-theta)) for old_theta,old_f,old_g in cuts)
            minimum_slack=min(minimum_slack,min(slacks))
            if min(slacks)<-2e-7:raise ArithmeticError(f'Convex subgradient control failed: {min(slacks)}.')
        cuts.append((theta.copy(),f,g.copy()))
        if best is None or rho<best['rho']-1e-10:
            best={'rho':rho,'parameters':theta.tolist(),'normalized_heights':(theta[:-1]/theta[-1]).tolist(),**detail};save()
            print(json.dumps({'query':oracle.queries,'rho':rho,'seconds':time.monotonic()-start,'cuts':len(cuts)}),flush=True)
        return rho,f,g,detail,physical
    raw=json.loads(args.starts.read_text());theta=np.array(raw['best']['parameters'] if isinstance(raw,dict) else raw[0])
    if len(theta)==5:theta=np.r_[theta[:4],np.zeros(dimension-5),theta[-1]]
    elif len(theta)==19 and args.full_core:theta=np.r_[theta[:11],np.zeros(3),theta[11:18],np.zeros(3),theta[-1]]
    if len(theta)!=dimension:raise ValueError('Wrong number of start parameters.')
    normal=np.zeros((2*(dimension-1),dimension+1))
    for j in range(dimension-1):normal[2*j,j]=1;normal[2*j,dimension-1]=args.height_bound;normal[2*j+1,j]=-1;normal[2*j+1,dimension-1]=args.height_bound
    scale=np.r_[np.ones(dimension-1),args.M];bounds=Bounds(np.r_[np.full(dimension-1,-60.),.002/args.M,-1e7],np.r_[np.full(dimension-1,60.),8/args.M,1e7])
    try:
        for restart in range(args.restarts):
            if restart:
                theta=np.array(best['parameters']);theta[:-1]+=rng.normal(0,.15,dimension-1);theta[-1]*=np.exp(rng.normal(0,.15))
                theta[:-1]=np.clip(theta[:-1],-.95*args.height_bound*theta[-1],.95*args.height_bound*theta[-1])
            rho,f,g,detail,physical=query(theta);anchor=theta.copy();H=detail['entropy'];current=(f+physical@anchor)/(2*H);mu=.03;serious=0
            for iteration in range(args.iterations):
                A=np.array([np.r_[-cg,1.] for _,_,cg in cuts]);lower=np.array([cf-cg@cp for cp,cf,cg in cuts]);
                def objective(x):
                    delta=(x[:dimension]-anchor)*scale
                    return (x[dimension]+physical@x[:dimension])/(2*H)+.5*mu*(delta@delta)
                def derivative(x):return np.r_[physical/(2*H)+mu*(scale**2)*(x[:dimension]-anchor),1/(2*H)]
                initial=np.r_[anchor,max(cf+cg@(anchor-cp) for cp,cf,cg in cuts)]
                opt=minimize(objective,initial,jac=derivative,method='SLSQP',constraints=[LinearConstraint(A,lower,np.inf),LinearConstraint(normal,0,np.inf)],bounds=bounds,options={'ftol':1e-11,'maxiter':300})
                candidate=opt.x[:dimension];model=(opt.x[dimension]+physical@candidate)/(2*H);predicted=current-model
                if not np.all(np.isfinite(candidate)) or np.min(A@opt.x-lower)<-1e-6:
                    records.append({'restart':restart,'iteration':iteration,'status':'QP unresolved','message':str(opt.message)});save();break
                if predicted<2e-9 or np.max(abs(candidate-anchor))<1e-9:
                    updated=query(anchor);newphysical=updated[4];change=max(abs(newphysical-physical))
                    physical=newphysical;H=updated[3]['entropy'];current=(updated[1]+physical@anchor)/(2*H)
                    records.append({'restart':restart,'iteration':iteration,'status':'type update at bundle stop','predicted':predicted,'type_change':change});save()
                    if change<1e-7:break
                    mu=.03;continue
                try:result=query(candidate)
                except ArithmeticError as error:
                    if 'Convex subgradient' in str(error):raise
                    oracle.failures+=1;mu=min(100.,mu*4);records.append({'restart':restart,'iteration':iteration,'status':'profile precision rejection'});save();continue
                value=(result[1]+physical@candidate)/(2*H);actual=current-value;accepted=actual>=.2*predicted
                if accepted:
                    anchor=candidate;current=value;serious+=1;mu=max(.0001,mu*.9)
                    if serious%6==0:
                        physical=result[4];H=result[3]['entropy'];current=(result[1]+physical@anchor)/(2*H)
                else:mu=min(10.,mu*1.15)
                records.append({'restart':restart,'iteration':iteration,'rho':result[0],'fixed_type_rho':value,'predicted':predicted,'actual':actual,'serious_step':bool(accepted),'proximal_weight':mu,'qp_success':bool(opt.success),'qp_iterations':int(opt.nit)});save()
                if best['rho']<1/np.sqrt(2)-1e-9:raise StopIteration
    except TimeoutError:records.append({'status':'time budget reached'})
    except StopIteration:records.append({'status':'below-threshold candidate needs exact validation'})
    finally:save()
    print(json.dumps({'best_rho':None if best is None else best['rho'],'queries':oracle.queries,'elapsed_seconds':time.monotonic()-start}),flush=True)

if __name__=='__main__':main()
