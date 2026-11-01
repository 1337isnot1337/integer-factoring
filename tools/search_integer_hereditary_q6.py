#!/usr/bin/env python3
"""Find an integer six-letter point of the hereditary profile relaxation.

The MILP is discovery only; exact cap and outward-pressure replay verifies
the retained point. No tensor or factoring improvement is constructed.
"""
import argparse,json,time
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
import numpy as np
from scipy.optimize import milp,LinearConstraint,Bounds
from scipy.special import logsumexp
from hereditary_adjugate_caps import constraints,inequalities
from check_cell_minor_caps import all_caps,pressure,Dyadic


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--max-order',type=int,default=80)
    parser.add_argument('--seconds',type=float,default=120)
    args=parser.parse_args()
    if args.max_order<1 or args.seconds<=0:parser.error('Use positive order and time budgets.')
    q=6;n=q*q;d=6;E=[0,0,12,16,18,24];K=args.max_order
    rows=[];rhs=[]
    def prefix(k):return np.r_[np.ones(k),np.zeros(n-k)]
    def put(row,b):rows.append(np.r_[row,np.zeros(n)]);rhs.append(float(b))
    for i in range(n-1):
     row=np.zeros(n);row[i]=1;row[i+1]=-1;put(row,0)
    for k,(cap,_) in all_caps(F(d),E,E).items():
     if k:put(prefix(k),cap)
    for cap in constraints(d,E,E):
     row=np.zeros(n)
     for k,c in cap['first_prefixes']+cap['second_prefixes']:row+=c*prefix(k)
     put(row,cap['rhs'])
    for i in range(n):
     for k in range(K):
      first=np.exp(-k/6);slope=np.exp(-(k+1)/6)-first
      row=np.zeros(2*n);row[i]=slope;row[n+i]=-1
      rows.append(row);rhs.append(slope*k-first)
    lower=np.zeros(2*n);upper=np.r_[np.full(n,K),np.ones(n)];upper[0]=0;lower[n]=1
    began=time.monotonic()
    result=milp(np.r_[np.zeros(n),np.ones(n)],integrality=np.r_[np.ones(n),np.zeros(n)],
     bounds=Bounds(lower,upper),constraints=LinearConstraint(np.asarray(rows),-np.inf,np.asarray(rhs)),
     options=dict(time_limit=args.seconds,mip_rel_gap=1e-8))
    print('solver',result.message,'seconds',time.monotonic()-began,flush=True)
    if result.x is not None:
     target=[int(round(v)) for v in result.x[:n]];assert target==sorted(target) and target[0]==0
     assert all(sum(target[:k])<=cap for k,(cap,_) in all_caps(F(d),E,E).items())
     caps=inequalities(d,target,target,E,E);assert all(c['slack']>=0 for c in caps)
     arithmetic=Dyadic(100);exact=pressure(arithmetic,F(d,6),[F(v,6) for v in target],[F(v,6) for v in E],[F(v,6) for v in E])
     score=d/6+2*logsumexp(-np.asarray(target)/6)-2*np.sqrt(2)*logsumexp(-np.asarray(E)/(6*np.sqrt(2)))
     report=dict(scope='Integer point of all compiled necessary inequalities; no realizing tensor.',deficit=d,energies=E,target=target,scale=[1,6],pressure_interval=list(exact),pressure_precision=100,numerical_pressure=float(score),minimum_cap_slack=str(min(c['slack'] for c in caps)),solver_message=str(result.message),seconds=time.monotonic()-began)
     report['max_order']=K
     report['source_sha256']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in (Path(__file__).name,'hereditary_adjugate_caps.py','check_cell_minor_caps.py','charged_gaussian_floor.py')}
     args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__ == "__main__":
    main()
