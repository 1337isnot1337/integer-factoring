#!/usr/bin/env python3
"""Search critical pressure after exact frozen-row coefficient cancellations.

This is an exploratory search, not a global obstruction or a pressure
certificate. At each selected pair of balanced charge gauges, it solves
F4 linear equations canceling entry coefficients below a specified
Gaussian ground-order margin. It enumerates small affine row spaces and
samples larger ones. Complete profiles include all Gaussian spectators.

The finite matrices are rebuilt from the proved polynomial construction.
Only exact polynomial arrays enter the simultaneous unit normalization.
"""
from __future__ import annotations
import argparse,json,math,time,itertools
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from native_matrix_smith_f4 import MatrixOracle
from check_rectangular_f4_repair import fmul,finv,ZERO
from frozen_core_row_profiles import combine_components,compile_interval

class CompiledRows:
 def __init__(self,M,beta=1):
  self.M=M;self.beta=beta
  self.data,self.compiled=compile_interval(M,beta)
  self.signs=[];self.vals=[];self.baseline=[]
  for side,(d,(comp,sh,md)) in enumerate(zip(self.data,self.compiled)):
   self.signs.append(np.array([1 if side==0 or j<M else -1 for j in range(2*M) if j not in d['active_indices']]))
   vals=[]
   for src in comp:
    vals.append(np.array([[sh[i][j]+((p[0]|p[1])&-(p[0]|p[1])).bit_length()-1 if p!=ZERO else np.inf for j,p in enumerate(row)] for i,row in enumerate(src)]))
   self.vals.append(np.array(vals))
 def oracle(self,g):
  result=[]
  selected=[i for i,c in enumerate(g) if c]
  for side,(comp,sh,md) in enumerate(self.compiled):
   result.append((MatrixOracle(combine_components(comp,g),len(comp[0]),sh,precision=md['profile_precision']),np.min(self.vals[side][selected],axis=0)))
  return result
 def logB(self,oracles,theta):
  total=0.;profiles=[]
  for side,((oracle,retained),(d,(_,_,md))) in enumerate(zip(oracles,zip(self.data,self.compiled))):
   rw=np.array(md['row_charges'])*theta[2*side];cw=np.array(md['column_charges'])*theta[2*side+1];res=oracle.profile(rw,cw,scale=theta[-1]);orders=np.array(res['orders'])
   if not md['exact_polynomial_matrix']:
    bound=min(md['builder_known_below']*theta[-1]+min(rw)+min(cw),np.min((retained+md['output_window'])*theta[-1]+rw[:,None]+cw[None,:]))
    if max(orders)>=bound-1e-8:raise ArithmeticError('Retained component precision is insufficient')
   total+=float(logsumexp(-orders));sp=np.array(d['spectator_orders']);values=sp*theta[-1]+self.signs[side]*(theta[2*side]+theta[2*side+1]);total+=float(np.sum(np.logaddexp(0,-values)));profiles.append(orders.tolist())
  return total,profiles
 def pressure(self,oracles,theta):
  logB,pro=self.logB(oracles,theta);R=math.sqrt(2);ci=theta[0]+theta[1];c0=theta[2]+theta[3];base=4*np.arange(self.M)-2*self.M+2;gx=base*theta[-1]+ci+c0;gy=base*theta[-1]+ci-c0
  partition=float(np.sum(np.logaddexp(0,-gx[1:-1]/R))+np.sum(np.logaddexp(0,-gy[1:-1]/R))-(gx[0]+gy[0])/R)
  return logB-R*partition,logB,pro
 def gaussian_check(self):
  from collections import Counter
  g=[1]+[0]*11;oracles=self.oracle(g);records=[]
  for side,((o,_),d) in enumerate(zip(oracles,self.data)):
   M=self.M;factors=list(range(5-2*M,2*M-4,2)) if side==0 else list(range(1,2*M-4,2))*2;ct=Counter(factors);sp=Counter(d['spectator_orders']);assert not (sp-ct)
   active=list((ct-sp).elements());expected=[4-4*M if side==0 else 0]
   for a in active:expected+= [v+a for v in expected]
   for a in (0,M/2,-M):
    md=self.compiled[side][2];actual=o.profile(np.array(md['row_charges'])*a,-np.array(md['column_charges'])*a)['orders'];assert actual==sorted(expected),(M,side,a,actual,expected)
   records.append({'place':d['place'],'orders':sorted(expected),'active_factors':active,'charges_checked':[0,M/2,-M]})
  return oracles,records

class AffineSpace:
 def __init__(self):self.rows={};self.add([int(i in (0,10)) for i in range(12)],1)
 def add(self,vec,rhs=0):
  a=list(vec)+[rhs]
  for p,b in sorted(self.rows.items()):
   if a[p]:
    z=a[p];a=[x^fmul(z,y) for x,y in zip(a,b)]
  pivot=next((i for i in range(12) if a[i]),None)
  if pivot is None:
   if a[-1]:raise ArithmeticError('inconsistent affine constraints')
   return
  z=finv(a[pivot]);a=[fmul(z,x) for x in a];self.rows[pivot]=a
 def solution(self,free_values):
  free=[i for i in range(12) if i not in self.rows];x=[0]*12
  for i,c in zip(free,free_values):x[i]=c
  for p,a in sorted(self.rows.items(),reverse=True):
   x[p]=a[-1]
   for j in range(p+1,12):x[p]^=fmul(a[j],x[j])
  return x
 def candidates(self,rng,cap=256):
  dim=12-len(self.rows)
  if 4**dim<=cap:return [self.solution(c) for c in itertools.product(range(4),repeat=dim)]
  starts=[[0]*dim]
  for i in range(dim):
   for c in (1,2,3):v=[0]*dim;v[i]=c;starts.append(v)
  starts=starts[:cap]
  starts += rng.integers(0,4,size=(max(0,cap-len(starts)),dim)).tolist()
  return [self.solution(v) for v in starts]

def leading_space(bank,a,b,margin=.25):
 space=AffineSpace();seen=set();details=[]
 for side,gauge in enumerate((a,b)):
  comp,sh,md=bank.compiled[side];r=np.array(md['row_charges']);c=np.array(md['column_charges']);tags=gauge*bank.M*(r[:,None]-c[None,:]);baseline=float(np.min(bank.vals[side][0]+tags));threshold=baseline-margin*bank.M;count=0
  for i in range(len(comp[0])):
   for j in range(len(comp[0])):
    limit=math.ceil(threshold-sh[i][j]-tags[i,j]-1e-9)
    if limit<=0:continue
    ps=[component[i][j] for component in comp];union=0
    for p in ps:union|=p[0]|p[1]
    if limit<union.bit_length():union&=(1<<limit)-1
    while union:
     bit=union&-union;union^=bit;k=bit.bit_length()-1;vec=tuple(((p[0]>>k)&1)|(((p[1]>>k)&1)<<1) for p in ps)
     if vec not in seen:
      assert vec[0]==0;seen.add(vec);space.add(vec);count+=1
     if len(space.rows)==12:break
    if len(space.rows)==12:break
   if len(space.rows)==12:break
  details.append({'place':md['place'],'baseline_minimum':baseline,'threshold':threshold,'new_constraints':count})
  if len(space.rows)==12:break
 return space,details

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--M',type=int,required=True);p.add_argument('--beta',type=int,default=1);p.add_argument('--grid',type=str,default='-2,-1,-0.5,0,0.5,1,2');p.add_argument('--cap',type=int,default=256);p.add_argument('--margin',type=float,default=.25);p.add_argument('--tau',type=float);p.add_argument('--output',type=Path,required=True);args=p.parse_args();start=time.monotonic();bank=CompiledRows(args.M,args.beta);base,controls=bank.gaussian_check();rng=np.random.default_rng(310827);grid=[float(v) for v in args.grid.split(',')];tau=args.tau if args.tau is not None else math.sqrt(2)*math.log(args.M);ell=tau/args.M;
 if args.cap<1 or args.margin<0 or ell<=0:raise ValueError('Use positive cap and scale, and nonnegative margin.')
 baseline,logB0,_=bank.pressure(base,np.array([0.,0.,0.,0.,ell]));L=args.M-2;h=lambda x:np.logaddexp(0,-x);theory=4*sum(h((2*j-1)*ell) for j in range(1,L+1))-4*math.sqrt(2)*sum(h(math.sqrt(2)*(2*j-1)*ell) for j in range(1,L//2+1));assert abs(baseline-theory)<1e-9
 records=[];best=None;cache={};queries=failures=0
 for a,b in itertools.product(grid,repeat=2):
  space,details=leading_space(bank,a,b,args.margin);candidates=space.candidates(rng,args.cap);current=None
  for g in candidates:
   key=tuple(g)
   if key not in cache:cache[key]=bank.oracle(g)
   theta=np.array([a*tau,-a*tau,b*tau,-b*tau,ell]);queries+=1
   try:pressure,logB,_=bank.pressure(cache[key],theta)
   except ArithmeticError:failures+=1;continue
   delta=logB-logB0;entry={'g':g,'a':a,'b':b,'critical_pressure':pressure,'defect_pressure':delta,'defect_over_ell':delta/ell}
   if current is None or pressure<current['critical_pressure']:current=entry
   if best is None or pressure<best['critical_pressure']-1e-11:
    best=entry;print(json.dumps({'M':args.M,'queries':queries,'spaces':len(records),'free_dimension':12-len(space.rows),**best}),flush=True)
  records.append({'a':a,'b':b,'free_dimension':12-len(space.rows),'candidates':len(candidates),'possible_candidates':4**(12-len(space.rows)),'constraints':details,'best':current})
  if len(records)%7==0:print(json.dumps({'progress_spaces':len(records),'queries':queries,'failures':failures,'seconds':time.monotonic()-start}),flush=True)
 result={'status':'exploratory critical-pressure screen with exact field cancellation constraints','M':args.M,'beta':args.beta,'normalized_margin':args.margin,'candidate_cap':args.cap,'random_seed':310827,'gauge_grid':grid,'tau':tau,'ell':ell,'gaussian_critical_pressure':baseline,'gaussian_theory_error':baseline-theory,'gaussian_profile_controls':controls,'best':best,'queries':queries,'failures':failures,'unique_rows':len(cache),'spaces':records,'seconds':time.monotonic()-start,'metadata':[{k:v for k,v in md.items() if k not in ('row_charges','column_charges')} for _,_,md in bank.compiled]};args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'done':True,'M':args.M,'queries':queries,'failures':failures,'unique_rows':len(cache),'best':best,'seconds':result['seconds']}),flush=True)
if __name__=='__main__':main()
