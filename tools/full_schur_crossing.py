#!/usr/bin/env python3
"""Optional NumPy replay of complete integrally normalized Schur pieces."""
from __future__ import annotations
import argparse,hashlib,json,math,time
from collections import Counter,defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
from check_half_piece_ranks import partitions,opcoef,afam,modrank
from check_schur_construction import schur_construct

@lru_cache(None)
def z(part):
 out=1
 for k,m in Counter(part).items():out*=k**m*math.factorial(m)
 return out
@lru_cache(None)
def characters(h):
 ps=partitions(h);rows=[]
 for lam in ps:
  poly=schur_construct(lam);row=[]
  for mu in ps:
   c=poly.get(mu,0)*z(mu);assert c.denominator==1 if isinstance(c,Fraction) else True;row.append(int(c))
  rows.append(row)
 for i,a in enumerate(rows):
  for j,b in enumerate(rows):assert sum(Fraction(x*y,z(mu)) for x,y,mu in zip(a,b,ps))==int(i==j)
 return rows

def mulmod(a,b,modulus):
 assert a.shape[1]==b.shape[0]
 shared=a.shape[1];bound=(modulus-1)**2
 batch=max(1,(2**62-1)//max(1,bound))
 if shared<=batch:return (a@b)%modulus
 out=np.zeros((a.shape[0],b.shape[1]),dtype=np.int64)
 for start in range(0,shared,batch):out=(out+a[:,start:start+batch]@b[start:start+batch,:])%modulus
 return out

def base_piece(h,s,modulus):
 ps=partitions(h);d=len(ps);channels=[(i,sigma) for i in range(2*h+1) for sigma in partitions(s+i)]
 left=np.zeros((d*d,len(channels)),dtype=np.int64);right=np.zeros((len(channels),d*d),dtype=np.int64)
 for c,(i,sigma) in enumerate(channels):
  for x,px in enumerate(ps):
   terms=[(mu,a*math.comb(sum(mu)+h,i)) for mu,a in afam(px) if sum(mu)+h>=i]
   for y,py in enumerate(ps):left[x*d+y,c]=((-1)**i*sum(a*opcoef(mu,sigma,py) for mu,a in terms))%modulus
  for u,pu in enumerate(ps):
   for v,pv in enumerate(ps):right[c,u*d+v]=opcoef(pv,pu,sigma)%modulus
 return mulmod(left,right,modulus)

def transform(tensor,matrix,axis,modulus):
 axes=[axis]+[i for i in range(4) if i!=axis];reverse=np.argsort(axes)
 shaped=np.transpose(tensor,axes);shape=shaped.shape
 flat=mulmod(matrix,shaped.reshape(shape[0],-1),modulus)
 return np.transpose(flat.reshape(shape),reverse)

def run(h,p,pieces_path=None):
 assert h>=1 and p>=2 and all(p%k for k in range(2,math.isqrt(p)+1))
 start=time.monotonic();ps=partitions(h);d=len(ps);fac=math.factorial(h);v=0;unit=fac
 while unit%p==0:unit//=p;v+=1
 power=p**(3*v);modulus=p*power
 assert modulus<2**31
 chi=characters(h);scale=np.array([[int(Fraction(fac*c,z(mu)))%modulus for c,mu in zip(row,ps)] for row in chi],dtype=np.int64);outchange=np.array(chi,dtype=np.int64)%modulus
 invunit=pow(unit**3%p,-1,p);total=np.zeros((d*d,d*d),dtype=np.int64);records=[];divchecks=0;binary_pieces=[]
 for s in range(-2*h,2*h+1):
  a=base_piece(h,s,modulus).reshape(d,d,d,d)
  for axis in (0,1,3):a=transform(a,scale,axis,modulus)
  a=transform(a,outchange,2,modulus)
  assert not np.any(a%power),f'integrality failed h={h},p={p},s={s}'
  divchecks+=a.size;a=(a//power*invunit)%p
  row=a.reshape(d*d,d*d);col=a.transpose(0,2,1,3).reshape(d*d,d*d)
  rr=modrank(row.tolist(),p);cc=modrank(col.tolist(),p);total=(total+row)%p
  records.append({'s':s,'row':rr,'column':cc})
  if p==2:binary_pieces.append([sum(int(v)<<j for j,v in enumerate(r)) for r in row.tolist()])
  print(json.dumps({'progress_h':h,'p':p,'s':s,'r':rr,'c':cc,'seconds':round(time.monotonic()-start,3)}),flush=True)
 expected=np.zeros((d*d,d*d),dtype=np.int64)
 for x in range(d):
  for y in range(d):expected[x*d+y,y*d+x]=1
 assert np.array_equal(total,expected)
 if pieces_path is not None:
  assert p==2
  pieces_path.write_text(json.dumps({'h':h,'q':d,'field':2,'shift_minimum':-2*h,'pieces':binary_pieces},indent=2)+'\n')
 S=sum(math.sqrt(r['row']*r['column']) for r in records);scale_exact=10**6;bound=sum((math.isqrt(r['row']*r['column']*scale_exact**2)+1) for r in records if r['row']*r['column'])
 upper=Fraction(bound,scale_exact);improve=(upper**10 < d**14) if d>1 else False
 return {'h':h,'prime':p,'dimension':d,'common_scale':1,'internal_cancellation_modulus':modulus,'integral_divisibility_checks':divchecks,'pieces':records,'rank_sum':S,'rho_rank_sum':math.log(S)/(2*math.log(d)) if d>1 else None,'strict_rho_7_10_certificate':improve,'seconds':time.monotonic()-start,'operator_cache':str(opcoef.cache_info())}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--weights',type=int,nargs='+',default=[1,2,3,4]);ap.add_argument('--prime',type=int,default=2);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--pieces',type=Path);a=ap.parse_args();rs=[]
 for h in a.weights:
  r=run(h,a.prime,a.pieces);rs.append(r);a.output.write_text(json.dumps({'records':rs,'scope':'Finite full-state Schur crossing ranks; all integrality and endpoint assertions checked. No improved exponent unless an explicit certificate is separately verified.'},indent=2)+'\n');print(json.dumps({'completed':r}),flush=True)
if __name__=='__main__':main()
