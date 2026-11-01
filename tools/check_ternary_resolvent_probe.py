#!/usr/bin/env python3
"""Independent ordinary-integer replay of selected ternary resolvent profiles."""
import argparse,json,time,hashlib
from pathlib import Path
def rank(A,p):
 A=[r[:] for r in A];m=len(A);n=len(A[0]) if m else 0;r=0
 for c in range(n):
  z=next((i for i in range(r,m) if A[i][c]%p),None)
  if z is None:continue
  A[r],A[z]=A[z],A[r];inv=pow(A[r][c],-1,p)
  A[r]=[(v*inv)%p for v in A[r]]
  for i in range(r+1,m):
   if A[i][c]%p:
    z=A[i][c];A[i]=[(v-z*w)%p for v,w in zip(A[i],A[r])]
  r+=1
  if r==m:break
 return r


def matmul(A,B):
 n=len(A);return [[sum(A[i][k]*B[k][j] for k in range(n))%3 for j in range(n)] for i in range(n)]
def build(U,V):
 q=D=3;n=27;B=[[0]*n for _ in range(n)]
 for x in range(q):
  for y in range(q):
   for a in range(D):
    r=(x*q+y)*D+a
    for u in range(q):
     for v in range(q):
      for b in range(D):
       c=(u*q+v)*D+b
       if y==v:B[r][c]=(int(x==b and a==u)+U[3*x+a]*V[3*b+u])%3
 N=[[0]*n for _ in range(n)]
 for r in range(n):
  for u in range(q):
   for v in range(q):
    for b in range(D):
     c=(u*q+v)*D+b;N[r][c]=(int(r==c)-B[r][(v*q+u)*D+b])%3
 C=B;coeff=[]
 for k in range(28):
  if not any(map(any,C)):break
  assert k<27
  coeff.append(C);C=matmul(N,C)
 L=len(coeff)-1
 for x in range(q):
  for y in range(q):
   for u in range(q):
    for v in range(q):
     for a in range(D):
      for b in range(D):
       r=(x*q+y)*D+a;c=(u*q+v)*D+b
       assert sum(M[r][c] for M in coeff)%3==int(x==v and y==u and a==b)
 return coeff

def flatten(coeff,cut_side,reverse=False):
 q=D=3;n=q*q;m=D*D*n;L=len(coeff)-1;cut,side=divmod(cut_side,2);out=[[[0]*(L+1) for _ in range(n)] for _ in range(m)]
 for x in range(q):
  for y in range(q):
   for u in range(q):
    for v in range(q):
     r,c=(x*q+y,u*q+v) if cut==0 else (x*q+u,y*q+v)
     if side:r,c=c,r
     for a in range(D):
      for b in range(D):
       values=[M[(x*q+y)*D+a][(u*q+v)*D+b] for M in coeff]
       out[D*D*r+D*a+b][c]=values[::-1] if reverse else values
 return out

def toeplitz(A,h):
 m,n=len(A),len(A[0]);L=len(A[0][0])-1;rows=[]
 for degree in range(h):
  for i in range(m):
   row=[0]*(h*n)
   for j in range(n):
    for power in range(min(L,degree)+1):row[(degree-power)*n+j]=A[i][j][power]
   if any(row):rows.append(row)
 return rank(rows,3)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--evidence',type=Path,default=Path(__file__).resolve().parents[1]/'evidence'/'ternary-resolvent-probe.json');ap.add_argument('--output',type=Path);a=ap.parse_args();start=time.monotonic();checks=[]
 evidence=json.loads(a.evidence.read_text());assert evidence['field']==3 and evidence['q']==3 and evidence['auxiliary']==3
 for rec in evidence['records']:
  C=build(rec['U'],rec['V']);assert len(C)-1==rec['degree']
  for profile in rec['profiles']:
   for reverse,orders in ((False,profile['nu']),(True,[rec['degree']-d for d in profile['d']])):
    A=flatten(C,profile['cut_side'],reverse)
    for h in sorted({1,2,3,max(orders)+1}):
     observed=toeplitz(A,h);expected=sum(max(h-x,0) for x in orders);assert observed==expected
     checks.append({'model':rec['model'],'cut_side':profile['cut_side'],'reverse':reverse,'precision':h,'rank':observed})
 out={'status':'passed','checks':checks,'count':len(checks),'seconds':time.monotonic()-start,'evidence_sha256':hashlib.sha256(a.evidence.read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Ordinary mod3 dense matrix products and Toeplitz ranks independently replay both full packed-polynomial Smith probes; no asymptotic improvement claimed.'}
 print(json.dumps(out,indent=2))
 if a.output:a.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
