#!/usr/bin/env python3
"""Explore exact monomial cancellation branches with convex inner optimization.

The finite-difference coefficient polynomials are exact over the prime
field. Distinct residue groups use their generic leading orders; rounded
integer replays combine all additional band ties before evaluation.
Optimization and pressure values are exploratory, not certificates.
"""
import argparse,itertools,math,json,time
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize,LinearConstraint,Bounds
from prime_state_filter_pressure import filter_pressure
from prime_state_filter_profiles import beta_polynomials


def partitions(n):
 def rec(a):
  if len(a)==n:yield tuple(a);return
  for g in range(max(a)+2):yield from rec(a+[g])
 yield from rec([0])


def branches(p,offsets=(-1,0,1)):
 for part in partitions(p):
  ng=max(part)+1;first=[part.index(g) for g in range(ng)]
  free=[i for i in range(p) if i not in first]
  for values in itertools.product(offsets,repeat=len(free)):
   shifts=[0]*p
   for i,m in zip(free,values):shifts[i]=m
   H=np.full((p,ng),np.inf)
   for k in range(p):
    polys=[{} for _ in range(ng)]
    for j in range(k+1):
     d=p-1-k+j;g=part[d];m=shifts[d]
     for h in range(j+1):
      r=h-m;polys[g][r]=(polys[g].get(r,0)+(-1)**(k-j+h)*math.comb(k,j)*math.comb(j,h))%p
    for g,pol in enumerate(polys):
     nz=[r for r,c in pol.items() if c]
     if nz:H[k,g]=min(nz)
   yield part,tuple(shifts),H


def solve_branch(p,n,J,t,part,shifts,H):
 ng=H.shape[1];nf=ng-1;k=np.arange(p);part=np.array(part);shifts=np.array(shifts)
 rows=[];bounds=[]
 for digit in range(p):
  for group in range(ng):
   if np.isfinite(H[digit,group]):
    row=np.zeros(nf+p);row[nf+digit]=1
    if group:row[group-1]=2
    rows.append(row);bounds.append(2*n*H[digit,group])
 constraint=LinearConstraint(np.array(rows),-np.inf,np.array(bounds))
 def fun(x):
  E=np.r_[0,x[:nf]];d=x[nf:]
  a=(-d+(1-J)*k)*t;b=-((1-2*n)*k+2*n*shifts+2*E[part])*t
  la=logsumexp(a);lb=logsumexp(b)
  grad=np.empty(nf+p);wb=np.exp(b-lb)
  for g in range(1,ng):grad[g-1]=-2*t*np.sum(wb[part==g])
  grad[nf:]=-t*np.exp(a-la)
  return la+lb,grad
 x0=np.r_[np.zeros(nf),2*n*np.min(H,axis=1)]
 lower=np.r_[np.full(nf,-4*n),np.full(p,-20*n)]
 upper=np.r_[np.full(nf,4*n),np.full(p,20*n)]
 opt=minimize(fun,x0,jac=True,method='SLSQP',constraints=[constraint],bounds=Bounds(lower,upper),options={'ftol':1e-10,'maxiter':100})
 E=np.r_[0,opt.x[:nf]]
 ds=np.min(2*n*H-2*E[None,:],axis=1)
 es=E[part]+n*shifts
 phi=logsumexp((-ds+(1-J)*k)*t)+logsumexp(-((1-2*n)*k+2*es)*t)
 return phi,es,ds,bool(opt.success)


def run_search(p,n,offset_values=(-1,0,1),refine=50,limit=None):
 start=time.monotonic();allbranches=list(branches(p,offset_values))
 if limit is not None:allbranches=allbranches[:limit]
 tau=math.log(2*n)/n
 grid=list(itertools.product([-n/2,-4.,0.,1.,4.,n/4,n/2],tau*np.array([.5,1.,2.])))
 records=[];failures=0;digits=np.arange(p)
 for index,(part,shifts,H) in enumerate(allbranches):
  current=None
  for J,t in grid:
   phi,e,d,ok=solve_branch(p,n,J,t,part,shifts,H);failures+=not ok
   infinity=(1-2*n)*digits+2*e
   value=filter_pressure(p,n,1,J,t,d,infinity)
   if current is None or value<current['pressure']:
    current={'partition':part,'offsets':shifts,'charge':float(J),'scale':float(t),'exponents_per_P':e.tolist(),'beta_orders_per_P':d.tolist(),'pressure':float(value),'optimizer_success':ok}
  records.append(current)
  if (index+1)%100==0:print(json.dumps({'branches':index+1,'optimizer_failures':failures,'seconds':time.monotonic()-start}),flush=True)
 lookup={(part,shifts):H for part,shifts,H in allbranches};replays=[]
 for entry in sorted(records,key=lambda value:value['pressure'])[:refine]:
  part=tuple(entry['partition']);shifts=tuple(entry['offsets']);H=lookup[part,shifts]
  def objective(value):
   J,t=float(value[0]),math.exp(float(value[1]))
   phi,e,d,ok=solve_branch(p,n,J,t,part,shifts,H)
   return filter_pressure(p,n,1,J,t,d,(1-2*n)*digits+2*e)
  bounds=[(-n/2,n/2),(math.log(.03),math.log(.6))]
  initial=[entry['charge'],max(bounds[1][0],min(bounds[1][1],math.log(entry['scale'])))]
  optimized=minimize(objective,initial,method='Nelder-Mead',bounds=bounds,options={'maxiter':120,'xatol':1e-6,'fatol':1e-9})
  J,t=round(float(optimized.x[0]),9),math.exp(float(optimized.x[1]));phi,e,d,ok=solve_branch(p,n,J,t,part,shifts,H)
  P=1000001
  e=[round(value/n)*n if abs(value/n-round(value/n))<1e-7 else value for value in e]
  exponents=tuple(round(P*value) for value in e)
  beta=tuple(min(value) for value in beta_polynomials(p,n,P,exponents))
  infinity=tuple((1-2*n)*P*k+2*value for k,value in enumerate(exponents))
  actual=filter_pressure(p,n,P,J,t,beta,infinity)
  replays.append({'partition':part,'offsets':shifts,'charge':J,'scale':t,'generic_branch_pressure':float(optimized.fun),'ramification':P,'exponents':exponents,'exact_beta_orders':beta,'exact_infinity_orders':infinity,'replay_pressure':float(actual),'optimizer_success':bool(optimized.success and ok)})
 return {'status':'Finite exploratory branch search; generic branch limits and exact integer replays are distinct.','prime':p,'n':n,'offset_values':offset_values,'branches':len(allbranches),'grid_pairs':len(grid),'grid':[(float(a),float(b)) for a,b in grid],'optimizer_failures':failures,'best_generic':min(records,key=lambda value:value['pressure']),'generic_records':records,'rounded_replays':replays,'best_replay':min(replays,key=lambda value:value['replay_pressure']) if replays else None,'seconds':time.monotonic()-start}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--prime',type=int,default=5);parser.add_argument('--n',type=int,default=32)
 parser.add_argument('--offsets',default='-1,0,1');parser.add_argument('--refine',type=int,default=50)
 parser.add_argument('--limit',type=int);parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();beta_polynomials(args.prime,args.n,1,(0,)*args.prime)
 if args.refine<0 or args.limit is not None and args.limit<1:raise ValueError('Use nonnegative refinement count and a positive optional limit.')
 result=run_search(args.prime,args.n,tuple(map(int,args.offsets.split(','))),args.refine,args.limit)
 args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({key:value for key,value in result.items() if key in ('prime','n','branches','optimizer_failures','best_replay','seconds')}))


if __name__=='__main__':main()
