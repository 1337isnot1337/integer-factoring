#!/usr/bin/env python3
"""Exact rational quantum-spin crossing paths at an integer square parameter."""
from fractions import Fraction as F
import argparse,json,time,math
from pathlib import Path
from fused_rmatrix_probe import mmul,poly_linear,arithmetic_initial_form,profile_initial,vp_int,vp

def quantum_source(k,t):
    if not isinstance(k,int) or k<1 or not isinstance(t,int) or t<2:
        raise ValueError('k must be a positive integer and t an integer at least two')
    Q=F(t*t);q=k+1
    def qnum(a):return (Q**a-Q**(-a))/(Q-Q**(-1))
    values={};norm=math.prod(Q**(2*l)-1 for l in range(1,k+1))
    for w in range(2*k+1):
        basis=[(i,w-i) for i in range(q) if 0<=w-i<q];n=len(basis);idx={v:i for i,v in enumerate(basis)}
        C=[[F(0) for _ in basis] for _ in basis]
        for col,(i,j) in enumerate(basis):
            C[col][col]+=qnum(k-w)*qnum(k-w+1)
            raises=[((i-1,j),qnum(i)*F(t)**(-k+2*j)),((i,j-1),qnum(j)*F(t)**(k-2*i))]
            for (ii,jj),e in raises:
                if not e:continue
                lowers=[((ii+1,jj),qnum(k-ii)*F(t)**(-k+2*jj)),((ii,jj+1),qnum(k-jj)*F(t)**(k-2*ii))]
                for dest,f in lowers:
                    if f and dest in idx:C[idx[dest]][col]+=e*f
        coeff=[[[F(0) for _ in basis] for _ in basis] for _ in range(k+1)]
        projector_sum=[[F(0) for _ in basis] for _ in basis]
        for spin in range(abs(k-w),k+1):
            eig=qnum(spin)*qnum(spin+1);P=[[F(i==j) for j in range(n)] for i in range(n)];den=F(1)
            for other in range(abs(k-w),k+1):
                if other==spin:continue
                lam=qnum(other)*qnum(other+1)
                P=mmul(P,[[C[i][j]-(lam if i==j else 0) for j in range(n)] for i in range(n)]);den*=eig-lam
            P=[[x/den for x in row] for row in P]
            assert mmul(C,P)==[[eig*x for x in row] for row in P]
            poly=[F(1)]
            for l in range(1,k+1):poly=poly_linear(poly,Q**(2*l)-1,-1 if l<=spin else Q**(2*l))
            for i in range(n):
                for j in range(n):
                    projector_sum[i][j]+=P[i][j]
                    for d,c in enumerate(poly):coeff[d][i][j]+=c*P[i][j]/norm
        assert projector_sum==[[F(i==j) for j in range(n)] for i in range(n)]
        for i,(x,y) in enumerate(basis):
            for j,(u,v) in enumerate(basis):
                entry=tuple(coeff[d][i][j] for d in range(k+1))
                if any(entry):values[y,x,u,v]=entry # uncheck by output swap
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):assert values.get((x,y,u,v),(0,))[0]==int(x==v and y==u)
    return values

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--q',default='3,4,5,6,7,8,9');ap.add_argument('--t',type=int,default=2);ap.add_argument('--primes',default='3,5,7,11,13');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    records=[];args.output.parent.mkdir(parents=True,exist_ok=True)
    for q in map(int,args.q.split(',')):
        start=time.monotonic();raw=quantum_source(q-1,args.t)
        print(json.dumps({'built_q':q,'t':args.t,'seconds':time.monotonic()-start}),flush=True)
        for p in map(int,args.primes.split(',')):
            fact=[0]
            for i in range(1,q):fact.append(fact[-1]+vp_int(i,p))
            Q=F(args.t*args.t);qfacts=[0]
            for i in range(1,q):qfacts.append(qfacts[-1]+vp((Q**i-Q**(-i))/(Q-Q**(-1)),p))
            qbinoms=[qfacts[-1]-qfacts[i]-qfacts[q-1-i] for i in range(q)]
            heights=[('zero',[0]*q,[0]*q)]
            for kind,h in [('factorial',fact),('quantum_factorial',qfacts),('quantum_binomial',qbinoms)]:
                for sign in [1,-1]:
                    for relative in [1,-1]:heights.append((f'{kind}:{sign}:{relative}',[sign*x for x in h],[relative*sign*x for x in h]))
            seen=set()
            for name,h,g in heights:
                if (tuple(h),tuple(g)) in seen:continue
                seen.add((tuple(h),tuple(g)))
                val,degree,slope=arithmetic_initial_form(raw,q,p,(h,g));rec=profile_initial(val,degree,q,p)
                rec.update(t=args.t,gauge=name,h=h,g=g,critical_slope=str(slope));records.append(rec)
                print(json.dumps({k:rec[k] for k in ['q','prime','gauge','degree','critical_slope','rho_float','certificate']}),flush=True)
                args.output.write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':main()
