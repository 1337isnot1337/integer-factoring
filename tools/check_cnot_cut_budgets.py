#!/usr/bin/env python3
"""Exact polynomial positivity certificates for column-cut budget gates."""
from math import comb
import argparse
import hashlib
import json
from pathlib import Path
from collections import Counter

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:c[i+j]+=x*y
 return c

def row_degrees(t):
 d={0:1}
 for _ in range(t):
  out=Counter()
  for k,v in d.items():
   for e,m in [(3,1),(0,2),(-3,1)]:out[k+e]+=m*v
  d=dict(out)
 return d

def complement(q,m):
 if m==0:return q*q
 if m in (1,2):return q*(q-1)//2
 if m>=2*q-1:return 0
 best=0
 for A in range(2,2*m-1):
  B=2*m-A
  for s in range(1,min(q-1,(2*q-2)//B,(2*q-A)//2)+1):
   u=min((2*q-2*s)//A,(2*q-B*s)//2)
   if u>0:best=max(best,s*u)
 return best

def geometric_witness(q,m,value):
 if m < 3 or not value:return None
 for A in range(2,2*m-1):
  B=2*m-A
  for s in range(1,min(q-1,(2*q-2)//B,(2*q-A)//2)+1):
   u=min((2*q-2*s)//A,(2*q-B*s)//2)
   if u>0 and s*u==value:
    assert 2*s+A*u<=2*q and B*s+2*u<=2*q
    return {'A':A,'B':B,'s':s,'u':u}
 raise AssertionError('missing geometric witness')

def gate_polynomial(t,T):
 q=2**t;n=q*q;L=3*t;d=row_degrees(t)
 g=[max(k,complement(q,min(n,sum(m*max(e-k+1,0)for e,m in d.items()))))for k in range(1,L+1)]
 G=list(reversed(g))
 v=mul([1,-2,1],G);v[1]+=n
 r=[0]*(6*t+1)
 for j in range(2*t+1):r[3*j]=comb(2*t,j)
 p=mul(r,v);p[2*L+1]-=T
 return p,g

def shift(p):
 return [sum(p[j]*comb(j,i)for j in range(i,len(p)))for i in range(len(p))]

def certificate(t,T,limit=10000):
 p,g=gate_polynomial(t,T)
 a=shift(p)
 for m in range(limit+1):
  if min(a)>0:
   q=2**t;L=3*t;n=q*q;ds=row_degrees(t)
   cuts=[]
   for k in range(1,L+1):
    r=min(n,sum(v*max(d-k+1,0)for d,v in ds.items()))
    c=complement(q,r)
    if c:
     cuts.append({'k':k,'row_tail_bound':r,'column_prefix_bound':c,
       'geometric_witness':geometric_witness(q,r,c)})
   assert T**35>q**99 and 99**2>8*35**2
   assert all(v>0 for v in a)
   digest=hashlib.sha256(json.dumps(a,separators=(',',':')).encode()).hexdigest()
   return {'t':t,'q':q,'degree':L,'T':T,'multiplier_power':m,
     'positive_polynomial_degree':len(a)-1,'smallest_coefficient':min(a),
     'positive_coefficients_sha256':digest,'column_cut_bounds':g,
     'complement_constraints':cuts,'target_comparison_verified':True}
  a=[a[0]]+[a[i]+a[i-1]for i in range(1,len(a))]+[a[-1]]
 raise RuntimeError('no certificate')

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path,default=Path('cnot-cut-budgets.json'))
 args=parser.parse_args()
 out=[]
 for t,T in [(4,2600),(5,19000),(6,130000)]:
  r=certificate(t,T);out.append(r);print(json.dumps(r),flush=True)
 result={'scope':'Exact ordinary-proof certificates for the zero-basepoint Smith criterion; no factoring exponent improvement.',
   'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
   'certificates':out}
 args.output.write_text(json.dumps(result,indent=2)+'\n')
