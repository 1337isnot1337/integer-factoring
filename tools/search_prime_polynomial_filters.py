#!/usr/bin/env python3
"""Search polynomial prime-state filters by exact residue-space feasibility.

The degree grid is finite. Four priority orders select Pareto candidates;
this does not enumerate every feasible profile or every polynomial.
Witness construction and profile orders are exact. Pressure is numerical.
"""
import argparse,itertools,math,time,json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from polynomial_state_filter_synthesis import synthesize
from prime_state_filter_pressure import gaussian_pressure,filter_pressure
from prime_state_filter_profiles import beta_polynomials
from search_prime_state_filters import digit_log


def greedy(prime,m,A,permutation):
 G=tuple(4*m+k*m for k in range(prime))
 B=list(G if synthesize(prime,m,A,G,construct=False) else (0,)*prime)
 calls=1
 for k in permutation:
  lo=B[k];hi=max(A[prime-1-k+j]+m*j for j in range(k+1))+1
  while hi-lo>1:
   mid=(lo+hi)//2;trial=B[:];trial[k]=mid;calls+=1
   if synthesize(prime,m,A,trial,construct=False):lo=mid
   else:hi=mid
  B[k]=lo
 return synthesize(prime,m,A,B),calls


def evaluate(p,n,P,A,B,grid):
 k=np.arange(p);gg=np.array(grid);J=gg[:,0,None];t=gg[:,1,None]
 values=gg[:,2]+logsumexp(((1-J)*k-2*np.array(B)/P)*t,axis=1)+logsumexp(-((1-2*n)*k-2*np.array(A)/P)*t,axis=1)
 i=int(np.argmin(values));return {'J':float(gg[i,0]),'t':float(gg[i,1]),'pressure':float(values[i])}


def run_search(p,n,P,limit=None):
 start=time.monotonic();m=n*P;L=4*m
 levels=sorted(set((-m,-(m//2),-2,-1,0,1,2,m//2,m)))
 degrees=[(L,)+tuple(L+d for d in row) for row in itertools.product(levels,repeat=p-1)]
 degrees.sort(key=lambda row:(sum(abs(value-L) for value in row),row))
 if limit is not None:degrees=degrees[:limit]
 alternating=[]
 for k in range((p+1)//2):
  alternating.append(p-1-k)
  if k!=p-1-k:alternating.append(k)
 priorities=list(dict.fromkeys((tuple(range(p)),tuple(reversed(range(p))),tuple(sorted(range(p),key=lambda i:(abs(i-(p-1)/2),i))),tuple(alternating))))
 scales=math.log(2*n)/n*np.array((.5,.75,1.,1.25,1.5,2.))
 charges=(-n/2,-8.,-4.,0.,1.,4.,8.,n/2)
 grid=[(J,float(t),float(gaussian_pressure(p,n,J,t)-digit_log(p,(2*n-1+J)*t)-digit_log(p,(1-2*n)*t))) for J,t in itertools.product(charges,scales)]
 best=None;records=[];seen=set();calls=constructions=0
 for index,A in enumerate(degrees):
  for priority in priorities:
   witness,count=greedy(p,m,A,priority);calls+=count;constructions+=1
   if witness is None:raise ArithmeticError('A greedy sequence lost feasibility.')
   actual_A=tuple(witness['actual_degrees']);actual_B=tuple(witness['actual_orders'])
   key=(actual_A,actual_B)
   if key in seen:continue
   seen.add(key);record=evaluate(p,n,P,actual_A,actual_B,grid)
   infinity=[P*(1-2*n)*k-2*degree for k,degree in enumerate(actual_A)]
   replay=filter_pressure(p,n,P,record['J'],record['t'],[2*order for order in actual_B],infinity)
   assert abs(replay-record['pressure'])<1e-9
   record.update({'pressure':replay,'degree_bounds':A,'priority':priority,'actual_degrees':actual_A,'actual_orders':actual_B,'witness':witness})
   records.append(record)
   if best is None or record['pressure']<best['pressure']-1e-9:
    best=record;print(json.dumps({'degree_vectors':index+1,**{k:v for k,v in record.items() if k!='witness'}}),flush=True)
  if (index+1)%250==0:print(json.dumps({'degree_vectors':index+1,'distinct_profiles':len(seen),'feasibility_decisions':calls,'seconds':time.monotonic()-start}),flush=True)
 return {'status':'Finite degree-grid and priority search; exact polynomial witnesses, numerical pressure, no obstruction claimed.','prime':p,'n':n,'ramification':P,'common_clearing':L,'degree_offsets':levels,'degree_vectors':len(degrees),'priorities':priorities,'charge_grid':charges,'scale_grid':scales.tolist(),'distinct_profiles':len(seen),'feasibility_decisions':calls,'witness_constructions':constructions,'best':best,'top':sorted(records,key=lambda row:row['pressure'])[:30],'seconds':time.monotonic()-start}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--prime',type=int,default=5);parser.add_argument('--n',type=int,default=32)
 parser.add_argument('--ramification',type=int,default=1);parser.add_argument('--limit',type=int)
 parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
 beta_polynomials(args.prime,args.n,args.ramification,(0,)*args.prime)
 if args.limit is not None and args.limit<1:raise ValueError('Use a positive optional degree-vector limit.')
 result=run_search(args.prime,args.n,args.ramification,args.limit)
 args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({key:value for key,value in result.items() if key in ('prime','n','ramification','degree_vectors','distinct_profiles','feasibility_decisions','seconds')}))


if __name__=='__main__':main()
