#!/usr/bin/env python3
"""Exact fused sl2 rational R-polynomials and finite-field Smith profiles.

Exploratory local script.  No rank or factoring claim is inferred from
floating-point minimization.  Every coefficient and Smith order is exact.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import factorial, log, exp, isqrt
from pathlib import Path


def mmul(a,b):
    n=len(a)
    out=[[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if a[i][k]:
                for j in range(n):
                    if b[k][j]: out[i][j]+=a[i][k]*b[k][j]
    return out


def poly_linear(poly,a,b):
    out=[0]*(len(poly)+1)
    for i,c in enumerate(poly):
        out[i]+=c*a
        out[i+1]+=c*b
    return out


def fused_polynomial(k):
    """R(u)=sum_j prod_(l<=j)(u+l)prod_(l>j)(u-l) P_j.

    P_j are the exact rational eigenspace projectors for total Casimir
    eigenvalue j(j+1) on Sym^k Q^2 tensor Sym^k Q^2.  Work in each total
    weight block to keep all rational matrices of size at most k+1.
    """
    q=k+1
    values={}
    denom=1
    for w in range(2*k+1):
        basis=[(i,w-i) for i in range(q) if 0<=w-i<q]
        n=len(basis); index={v:i for i,v in enumerate(basis)}
        C=[[Fraction(0) for _ in basis] for _ in basis]
        # C=(H_total/2)(H_total/2+1)+F_total E_total.
        for col,(i,j) in enumerate(basis):
            C[col][col]+=(k-w)*(k-w+1)
            for (ii,jj),e in (((i-1,j),i),((i,j-1),j)):
                if not e: continue
                for result,f in (((ii+1,jj),k-ii),((ii,jj+1),k-jj)):
                    if f and result in index: C[index[result]][col]+=e*f
        coeff=[[[Fraction(0) for _ in basis] for _ in basis] for _ in range(k+1)]
        for spin in range(abs(k-w),k+1):
            eigen=spin*(spin+1)
            P=[[Fraction(i==j) for j in range(n)] for i in range(n)]
            denominator=1
            for other in range(abs(k-w),k+1):
                if other==spin: continue
                lam=other*(other+1)
                shifted=[[C[i][j]-(lam if i==j else 0) for j in range(n)] for i in range(n)]
                P=mmul(P,shifted); denominator*=eigen-lam
            polynomial=[1]
            for l in range(1,k+1): polynomial=poly_linear(polynomial,l if l<=spin else -l,1)
            for degree,c in enumerate(polynomial):
                for i in range(n):
                    for j in range(n): coeff[degree][i][j]+=Fraction(c,denominator)*P[i][j]
        for i,(x,y) in enumerate(basis):
            for j,(u,v) in enumerate(basis):
                entry=tuple(coeff[d][i][j] for d in range(k+1))
                if any(entry): values[x,y,u,v]=entry
                for c in entry: denom=__import__('math').lcm(denom,c.denominator)
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    entry=values.get((x,y,u,v),(0,)*(k+1))
                    assert entry[0]==factorial(k)*int(x==v and y==u)
                    assert entry[-1]==int(x==u and y==v)
    return values,denom


def ptrim(a):
    while a and a[-1]==0: a.pop()
    return a


def psub(a,b,p):
    out=a+[0]*max(0,len(b)-len(a))
    for i,c in enumerate(b): out[i]=(out[i]-c)%p
    return ptrim(out)


def pmul(a,b,p,P):
    if not a or not b:return []
    out=[0]*min(P,len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:len(out)-i]): out[i+j]=(out[i+j]+x*y)%p
    return ptrim(out)


def pinv(a,p,P):
    assert a and a[0]
    inv=pow(a[0],-1,p); out=[inv]
    for i in range(1,P):out.append(-inv*sum(a[j]*out[i-j] for j in range(1,min(i+1,len(a))))%p)
    return ptrim(out)


def smith(matrix,p,precision):
    a=[[ptrim(entry[:precision]) for entry in row] for row in matrix]
    n=len(a); orders=[]
    for k in range(n):
        best=None
        for i in range(k,n):
            for j in range(k,n):
                if a[i][j]:
                    v=next(h for h,x in enumerate(a[i][j]) if x)
                    if best is None or v<best[0]:best=(v,i,j)
        if best is None: raise ValueError('precision exhausted or singular')
        order,i,j=best;orders.append(order)
        a[k],a[i]=a[i],a[k]
        for row in a:row[k],row[j]=row[j],row[k]
        inverse=pinv(a[k][k][order:],p,precision)
        for j in range(k,n):a[k][j]=pmul(a[k][j],inverse,p,precision)
        assert a[k][k]==[0]*order+[1]
        for i in range(k+1,n):
            factor=a[i][k][order:]
            for j in range(k+1,n):a[i][j]=psub(a[i][j],pmul(factor,a[k][j],p,precision),p)
            a[i][k]=[]
    return orders


def translated(entry,center,p):
    out=[0]*len(entry)
    for degree,co in enumerate(entry):
        if co:
            b=1
            for h in range(degree+1):
                out[h]=(out[h]+co*b*pow(center,degree-h,p))%p
                if h<degree:b=b*(degree-h)//(h+1)
    return ptrim(out)


def cut_blocks(values,q,cut):
    blocks=[]
    if cut==0:
        for w in range(2*q-1):
            basis=[(i,w-i) for i in range(q) if 0<=w-i<q]
            blocks.append([[list(values.get((x,y,u,v),[])) for u,v in basis] for x,y in basis])
    else:
        for w in range(1-q,q):
            left=[(x,u) for x in range(q) for u in range(q) if x-u==w]
            right=[(y,v) for y in range(q) for v in range(q) if v-y==w]
            blocks.append([[list(values.get((x,y,u,v),[])) for y,v in right] for x,u in left])
    return blocks


def optimize(degrees,orders):
    def score(t):
        # exponents small in recorded cases
        return sum(exp(t*d) for d in degrees)*sum(exp(-t*d) for d in orders)
    lo=0.;hi=10.
    for _ in range(120):
        a=(2*lo+hi)/3;b=(lo+2*hi)/3
        if score(a)<score(b):hi=b
        else:lo=a
    t=(lo+hi)/2
    return score(t),exp(t)


def objective_lower_certificate(degrees,orders,q):
    """An all-s>=1 integer AM-GM lower bound; no numerical optimization."""
    coefficients=Counter(a-b for a in degrees for b in orders)
    lower=coefficients[0]+sum(c+coefficients[-a] if c>=coefficients[-a]
        else isqrt(4*c*coefficients[-a]) for a,c in coefficients.items() if a>0)
    # sqrt(2)<99/70, so q^(2sqrt(2))<q^(99/35).
    return {'integer_lower_bound':lower,'exceeds_half_target':lower**35>q**99}


def probe(k,primes):
    raw,denom=fused_polynomial(k);q=k+1
    print('built',q,'denominator',denom,flush=True)
    records=[]
    for p in primes:
        if k>=p or denom%p==0: continue
        values={v:ptrim([int(x.numerator)*pow(int(x.denominator),-1,p)%p for x in co]) for v,co in raw.items()}
        # column orders at z=0 for z^k R(a(1-z)/z) are the reverse-R orders;
        # nonzero rescaling and translation in u do not change these orders.
        col_blocks=cut_blocks(values,q,1)
        orders=[]
        for block in col_blocks:
            reverse=[[(entry+[0]*(k+1-len(entry)))[::-1] for entry in row] for row in block]
            orders+=smith(reverse,p,len(block)*k+1)
        for a in range(1,p):
            degrees=[]
            for block in cut_blocks(values,q,0):
                translated_block=[[translated(entry,-a,p) for entry in row] for row in block]
                local_orders=smith(translated_block,p,len(block)*k+1)
                degrees += [k-v for v in local_orders]
            B,s=optimize(degrees,orders);rho=log(B)/(4*log(q))
            record={'q':q,'prime':p,'a':a,'degree_profile':dict(sorted(Counter(degrees).items())),'column_orders':dict(sorted(Counter(orders).items())),'B_float':B,'s_float':s,'rho_float':rho}
            records.append(record)
            print(json.dumps(record),flush=True)
    return records


def probe_mobius(k,primes):
    raw,denom=fused_polynomial(k);q=k+1
    records=[]
    print('mobius built',q,'denominator',denom,flush=True)
    for p in primes:
        if k>=p or denom%p==0: continue
        values={v:ptrim([int(x.numerator)*pow(int(x.denominator),-1,p)%p for x in co]) for v,co in raw.items()}
        profiles={}
        for cut in (0,1):
            blocks=cut_blocks(values,q,cut)
            for at in list(range(1,p))+['infinity']:
                orders=[]
                for block in blocks:
                    if at=='infinity':
                        local=[[(entry+[0]*(k+1-len(entry)))[::-1] for entry in row] for row in block]
                    else:
                        local=[[translated(entry,at,p) for entry in row] for row in block]
                    orders+=smith(local,p,len(block)*k+1)
                profiles[cut,at]=orders
        candidates=[]
        points=list(range(1,p))+['infinity']
        for first in points:
            for second in points:
                if first==second:continue
                for cut in (0,1):
                    degrees=[k-v for v in profiles[cut,first]]
                    orders=profiles[1-cut,second]
                    B,s=optimize(degrees,orders);rho=log(B)/(4*log(q))
                    candidates.append({'q':q,'prime':p,'infinity_point':first,'zero_point':second,'first_cut':cut,'degree_profile':dict(sorted(Counter(degrees).items())),'column_orders':dict(sorted(Counter(orders).items())),'B_float':B,'s_float':s,'rho_float':rho,'certificate':objective_lower_certificate(degrees,orders,q)})
        candidate=min(candidates,key=lambda x:x['B_float'])
        candidate['parameter_pairs_checked']=len(candidates)
        candidate['parameter_pairs_excluded_exactly']=sum(c['certificate']['exceeds_half_target'] for c in candidates)
        records.append(candidate)
        print(json.dumps(candidate),flush=True)
    return records


def vp_int(n,p):
    n=abs(n); v=0
    assert n
    while n%p==0:n//=p;v+=1
    return v


def vp(co,p):
    return vp_int(co.numerator,p)-vp_int(co.denominator,p)


def punit(co,p):
    n,d=co.numerator,co.denominator
    while n%p==0:n//=p
    while d%p==0:d//=p
    return n*pow(d,-1,p)%p


def arithmetic_initial_form(raw,q,p,heights=None):
    """Reduce a rational polynomial with P(0)=crossing at a valuation face.

    This performs rational arithmetic only; no p-adic or ramified-field
    package is needed.  Rational heights and a rational critical slope are
    legal, with their common denominator providing a ramification degree.
    """
    assert p>=2 and all(p%d for d in range(2,isqrt(p)+1))
    raw={key:tuple(Fraction(c) for c in entry) for key,entry in raw.items()}
    assert all(len(key)==4 and all(isinstance(x,int) and 0<=x<q for x in key) for key in raw)
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    entry=raw.get((x,y,u,v),(Fraction(0),))
                    assert entry and entry[0]==int(x==v and y==u)
    if heights is None:heights=([0]*q,[0]*q)
    hx,hy=([Fraction(h) for h in row] for row in heights)
    assert len(hx)==len(hy)==q
    slope=None
    for (x,y,u,v),entry in raw.items():
        weight=hx[x]+hy[y]-hy[u]-hx[v]
        for d,co in enumerate(entry):
            if not co or not d:continue
            value=(-vp(co,p)-weight)/d
            if slope is None or value>slope:slope=value
    if slope is None:slope=Fraction(0)
    result={}
    for (x,y,u,v),entry in raw.items():
        weight=hx[x]+hy[y]-hy[u]-hx[v]
        out=[]
        for d,co in enumerate(entry):
            if not co:out.append(0);continue
            valuation=vp(co,p)+weight+slope*d
            assert valuation>=0
            out.append(punit(co,p) if valuation==0 else 0)
        out=ptrim(out)
        if out:result[x,y,u,v]=out
    degree=max(map(len,result.values()))-1
    for key,entry in result.items():
        x,y,u,v=key
        assert entry[0]==int(x==v and y==u)
    return result,degree,slope


def initial_form(raw,k,p,heights=None):
    return arithmetic_initial_form({key:tuple(c/factorial(k) for c in entry) for key,entry in raw.items()},k+1,p,heights)


def profile_initial(values,degree,q,p):
    profiles={}
    for cut in (0,1):
        blocks=cut_blocks(values,q,cut)
        for at in list(range(1,p))+['infinity']:
            orders=[]
            for block in blocks:
                if at=='infinity':local=[[(entry+[0]*(degree+1-len(entry)))[::-1] for entry in row] for row in block]
                else:local=[[translated(entry,at,p) for entry in row] for row in block]
                orders+=smith(local,p,len(block)*degree+1)
            profiles[cut,at]=orders
    candidates=[]
    points=list(range(1,p))+['infinity']
    for first in points:
        for second in points:
            if first==second:continue
            for cut in (0,1):
                degrees=[degree-v for v in profiles[cut,first]]
                orders=profiles[1-cut,second]
                B,s=optimize(degrees,orders);rho=log(B)/(4*log(q))
                candidates.append({'q':q,'prime':p,'degree':degree,'infinity_point':first,'zero_point':second,'first_cut':cut,'degree_profile':dict(sorted(Counter(degrees).items())),'column_orders':dict(sorted(Counter(orders).items())),'B_float':B,'s_float':s,'rho_float':rho,'certificate':objective_lower_certificate(degrees,orders,q)})
    result=min(candidates,key=lambda x:x['B_float'])
    result['parameter_pairs_checked']=len(candidates)
    result['parameter_pairs_excluded_exactly']=sum(c['certificate']['exceeds_half_target'] for c in candidates)
    return result


def probe_initial(k,primes):
    raw,denom=fused_polynomial(k);q=k+1
    print('initial built',q,'denominator',denom,flush=True)
    records=[]
    for p in primes:
        if p>q:continue
        values,degree,slope=initial_form(raw,k,p)
        record=profile_initial(values,degree,q,p);record['slope']=str(slope)
        records.append(record);print(json.dumps(record),flush=True)
    return records


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-q',type=int,default=8);parser.add_argument('--output',default='fused-rmatrix-screen.json');parser.add_argument('--mobius',action='store_true');parser.add_argument('--initial',action='store_true');args=parser.parse_args()
    records=[]
    for k in range(1,args.max_q):
        records+=(probe_initial if args.initial else probe_mobius if args.mobius else probe)(k,(2,3,5,7,11,13))
        Path(args.output).write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':main()
