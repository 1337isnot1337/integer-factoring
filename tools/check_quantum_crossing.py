#!/usr/bin/env python3
"""Independent highest-weight control for rational quantum-spin paths."""
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path
from quantum_crossing_initial_forms import quantum_source
from fused_rmatrix_probe import arithmetic_initial_form,cut_blocks,smith,translated,profile_initial
from check_arithmetic_crossing import toeplitz_rank

def independent_check(k,t,raw):
 q=k+1;Q=F(t*t)
 def bracket(a):return (Q**a-Q**(-a))/(Q-Q**(-1))
 norm=math.prod(Q**(2*l)-1 for l in range(1,k+1))
 vectors=[];checks=0
 for spin in range(k+1):
  w=k-spin;c=F(1);v={}
  for i in range(w+1):
   v[i,w-i]=c
   if i<w:c=-c*bracket(w-i)/bracket(i+1)*F(t)**(2*(spin+1))
  # Direct E check, independent of the total-Casimir construction.
  e={}
  for (i,j),a in v.items():
   if i:e[i-1,j]=e.get((i-1,j),0)+a*bracket(i)*F(t)**(-k+2*j)
   if j:e[i,j-1]=e.get((i,j-1),0)+a*bracket(j)*F(t)**(k-2*i)
  assert all(x==0 for x in e.values())
  spectral=[F(1)]
  for l in range(1,k+1):
   a=Q**(2*l)-1;b=-1 if l<=spin else Q**(2*l);new=[F(0)]*(len(spectral)+1)
   for d,x in enumerate(spectral):new[d]+=a*x;new[d+1]+=b*x
   spectral=new
  spectral=[x/norm for x in spectral]
  for depth in range(2*spin+1):
   assert v
   vectors.append(v)
   w=k-spin+depth
   for x in range(q):
    y=w-x
    if not 0<=y<q:continue
    for d,target in enumerate(spectral):
     got=sum(raw.get((y,x,u,vv),(F(0),)*(k+1))[d]*a for (u,vv),a in v.items())
     assert got==target*v.get((x,y),F(0)),(q,t,spin,depth,x,y,d)
     checks+=1
   nv={}
   for (i,j),a in v.items():
    if i<k:nv[i+1,j]=nv.get((i+1,j),F(0))+a*bracket(k-i)*F(t)**(-k+2*j)
    if j<k:nv[i,j+1]=nv.get((i,j+1),F(0))+a*bracket(k-j)*F(t)**(k-2*i)
   v={ij:a for ij,a in nv.items() if a}
  assert not v
 assert len(vectors)==q*q
 # Modular independent-basis check. A nonzero modular determinant proves
 # rational independence; no source projectors are used here.
 p=1000003;rows=[]
 for v in vectors:
  rows.append([int(v.get((i,j),0).numerator)*pow(int(v.get((i,j),0).denominator),-1,p)%p for i in range(q) for j in range(q)])
 n=q*q
 for col in range(n):
  pos=next((i for i in range(col,n) if rows[i][col]),None);assert pos is not None
  rows[col],rows[pos]=rows[pos],rows[col];inv=pow(rows[col][col],-1,p)
  for i in range(col+1,n):
   scale=rows[i][col]*inv%p
   if scale:rows[i]=[(a-scale*b)%p for a,b in zip(rows[i],rows[col])]
 return {'q':q,'t':t,'highest_weight_vectors':n,'spectral_coefficient_checks':checks,'basis_prime':p}


def initial_control(q,t,p):
 raw=quantum_source(q-1,t)
 values,degree,slope=arithmetic_initial_form(raw,q,p)
 endpoint_checks=0
 for x in range(q):
  for y in range(q):
   for u in range(q):
    for v in range(q):
     entry=values.get((x,y,u,v),[])
     assert (entry[0] if entry else 0)==int(x==v and y==u)
     endpoint_checks+=1
 count=0;profiles=[]
 for cut in (0,1):
  for point in (1,'infinity'):
   orders=[]
   for block in cut_blocks(values,q,cut):
    if point=='infinity':local=[[(entry+[0]*(degree+1-len(entry)))[::-1] for entry in row] for row in block]
    else:local=[[translated(entry,point,p) for entry in row] for row in block]
    spectrum=smith(local,p,len(local)*degree+1)
    for precision in range(1,max(spectrum)+2):
     assert toeplitz_rank(local,precision,p)==sum(max(0,precision-v) for v in spectrum)
     count+=1
    orders+=spectrum
   profiles.append({'cut':cut,'point':point,'orders':sorted(orders)})
 best=profile_initial(values,degree,q,p)
 assert best['parameter_pairs_checked']==best['parameter_pairs_excluded_exactly']
 return {'q':q,'t':t,'prime':p,'degree':degree,'critical_slope':str(slope),
         'endpoint_checks':endpoint_checks,'independent_toeplitz_checks':count,
         'profiles':profiles,'best_local_profile':best}


def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--output',type=Path)
 args=parser.parse_args()
 records=[]
 for t in (2,3):
  for q in range(2,10):
   raw=quantum_source(q-1,t)
   records.append(independent_check(q-1,t,raw))
 initial=[initial_control(*case) for case in [(3,2,2),(3,2,3),(4,2,5),(5,2,7),(6,3,2)]]
 failures=0
 for k,t in [(0,2),(2,0),(2,1)]:
  try:quantum_source(k,t)
  except ValueError:failures+=1
  else:raise AssertionError('invalid source parameter accepted')
 report={'scope':'Exact rational quantum-spin sources and selected arithmetic initial forms; no improved factoring exponent.',
         'source_controls':records,'initial_form_controls':initial,'invalid_parameter_controls':failures}
 text=json.dumps(report,indent=2)+'\n'
 if args.output:args.output.write_text(text)
 print(text,end='')
if __name__=='__main__':main()
