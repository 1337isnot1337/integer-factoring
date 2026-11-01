#!/usr/bin/env python3
"""Search exact valuation initial forms of higher-spin crossing paths."""
from collections import Counter
from fractions import Fraction
import argparse,json,random,time
from math import factorial
from pathlib import Path
from fused_rmatrix_probe import fused_polynomial,vp,punit,ptrim,profile_initial


def main():
    a=argparse.ArgumentParser();a.add_argument('--q',type=int,default=6);a.add_argument('--p',type=int,default=2);a.add_argument('--count',type=int,default=2000);a.add_argument('--seed',type=int,default=18931);a.add_argument('--output',default='arithmetic-gauges.json');args=a.parse_args()
    q,p=args.q,args.p;k=q-1
    raw,denominator=fused_polynomial(k);f=factorial(k)
    entries=[]
    for (x,y,u,v),polynomial in raw.items():
        for d,c in enumerate(polynomial):
            if c:entries.append((x,y,u,v,d,vp(c/f,p),punit(c/f,p)))
    rng=random.Random(args.seed);seen=set();best=None;records=[];degrees=Counter();slopes=Counter();start=time.time()
    for iteration in range(args.count):
        if iteration==0:hx=hy=[0]*q
        elif iteration<2*q:
            hx=[int(i==iteration//2) for i in range(q)];hy=hx[:] if iteration%2 else [0]*q
        else:
            extent=1+(iteration%4==0)
            hx=[rng.randrange(-extent,extent+1) for _ in range(q)]
            hy=hx[:] if iteration%3==0 else [rng.randrange(-extent,extent+1) for _ in range(q)]
        slope=max(Fraction(-valuation-hx[x]-hy[y]+hy[u]+hx[v],d) for x,y,u,v,d,valuation,unit in entries if d)
        chosen=[]
        for x,y,u,v,d,valuation,unit in entries:
            vv=valuation+hx[x]+hy[y]-hy[u]-hx[v]+d*slope
            assert vv>=0
            if vv==0:chosen.append((x,y,u,v,d,unit))
        key=tuple(chosen)
        if key in seen:continue
        seen.add(key)
        values={};degree=max(row[4] for row in chosen)
        for x,y,u,v,d,unit in chosen:
            entry=values.setdefault((x,y,u,v),[0]*(degree+1));entry[d]=unit
        values={key:ptrim(entry) for key,entry in values.items()}
        degrees[degree]+=1;slopes[str(slope)]+=1
        result=profile_initial(values,degree,q,p)
        result.update({'slope':str(slope),'hx':hx,'hy':hy,'iteration':iteration})
        if best is None or result['B_float']<best['B_float']-1e-8:
            best=result;records.append(result);print(json.dumps(result),flush=True)
        if len(seen)%100==0:
            print('progress',iteration,'unique',len(seen),'elapsed',time.time()-start,flush=True)
        report={'q':q,'p':p,'iterations':iteration+1,'seed':args.seed,'unique_initial_forms':len(seen),'degree_counts':dict(degrees),'slope_counts':dict(slopes),'improvements':records,'best':best,'scope':'finite exploratory screen; floating objective used only for ranking; no improved exponent claimed'}
        Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print('done',args.count,'unique',len(seen),'elapsed',time.time()-start,flush=True)


if __name__=='__main__':main()
