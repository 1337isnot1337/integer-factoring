#!/usr/bin/env python3
"""Prime-field trace-shear crossing, including asymmetric odd-field cut ranks."""
import argparse,json,time
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

def run(p,q,t,direct=False):
 t%=p;kappa=(1+t*q)%p;assert t and kappa
 ti=pow(t,-1,p);ki=pow(kappa,-1,p);n=q*q;eta=[int(i//q==i%q) for i in range(n)]
 def shear(z,a):
  tr=sum(z[q*i+i] for i in range(q))%p
  return [(x+a*tr*y)%p for x,y in zip(z,eta)]
 def units(x,y,v,u,z):
  out=[0]*n;out[q*x+v]=z[q*y+u];return out
 checks=0
 for x in range(q):
  for y in range(q):
   for u in range(q):
    for v in range(q):
     first=units(x,y,v,u,eta)
     second=shear(units(x,u,v,y,shear(eta,-t*ki)),t)
     total=[(-ti*a+kappa*ti*b)%p for a,b in zip(first,second)]
     assert total==([int(x==v and y==u)*a for a in eta]);checks+=1
 def F(a):return [[(int(i==j)+a*eta[i]*eta[j])%p for j in range(n)] for i in range(n)]
 def G(a):return [[(int(i//q==j%q and i%q==j//q)+a*int(i==j))%p for j in range(n)] for i in range(n)]
 assert rank(F(t),p)==rank(F(-t*ki),p)==n
 c1=n*rank(G(t),p);r2=n*rank(G(-t*ki),p)
 directranks=[]
 if direct:
  for left,right,expected in ((-t*ki,t,c1),(t,-t*ki,r2)):
   rows=[[0]*(q**4) for _ in range(q**4)]
   for x in range(q):
    for y in range(q):
     for u in range(q):
      for v in range(q):
       for i in range(q):
        for j in range(q):
         z=[0]*n;z[q*i+j]=1
         out=shear(units(x,y,v,u,shear(z,right)),left)
         for a in range(q):
          for b in range(q):rows[((x*q+u)*q+a)*q+i][((y*q+v)*q+b)*q+j]=out[q*a+b]
   rr=rank(rows,p);assert rr==expected;directranks.append(rr)
 return {'p':p,'q':q,'t':t,'kappa':kappa,'auxiliary_order':n,'row_profile':[1,r2],'column_profile':[c1,1],'direct_opposite_ranks_c1_r2':directranks,'endpoint_checks':checks}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);a=ap.parse_args();start=time.monotonic();rs=[run(3,3,2,True),run(3,6,2),run(3,5,2),run(3,8,2),run(5,5,4)];out={'status':'passed','records':rs,'seconds':time.monotonic()-start,'scope':'General trace-shear catalytic crossing, exact prime-field endpoint and opposite ranks. Odd characteristic has asymmetric ranks unless special parameter congruences hold.'};print(json.dumps(out,indent=2))
 if a.output:a.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
