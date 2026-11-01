#!/usr/bin/env python3
"""Optional exact finite homogeneous Schur profiles at binary power weights."""
from __future__ import annotations
import argparse,json,math,time
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
from full_schur_crossing import partitions,opcoef,characters,z,mulmod,transform

def rank2(matrix):
 pivots={}
 for row in matrix:
  n=int.from_bytes(np.packbits(row.astype(np.uint8),bitorder='little').tobytes(),'little')
  while n:
   k=n.bit_length()-1
   if k in pivots:n^=pivots[k]
   else:pivots[k]=n;break
 return len(pivots)

def inverse2(matrix):
 n=len(matrix);rows=[sum(int(x)<<j for j,x in enumerate(row))|(1<<(n+i)) for i,row in enumerate(matrix)]
 for k in range(n):
  pivot=next((i for i in range(k,n) if rows[i]>>k&1),None)
  assert pivot is not None,'singular full homogeneous Gram'
  rows[k],rows[pivot]=rows[pivot],rows[k]
  for i in range(n):
   if i!=k and rows[i]>>k&1:rows[i]^=rows[k]
 return np.array([[(row>>(n+j))&1 for j in range(n)] for row in rows],dtype=np.int64)

def base(h,s,modulus,exchange):
 ps=partitions(h);d=len(ps);indices=[0,2*h];channels=[(i,sigma) for i in indices for sigma in partitions(2*h-s-i if exchange else s+i)]
 left=np.zeros((d*d,len(channels)),dtype=np.int64);right=np.zeros((len(channels),d*d),dtype=np.int64)
 for c,(i,sigma) in enumerate(channels):
  for x,px in enumerate(ps):
   for y,py in enumerate(ps):left[x*d+y,c]=opcoef(px,py,sigma)%modulus if exchange else opcoef(px,sigma,py)%modulus
  for u,pu in enumerate(ps):
   for v,pv in enumerate(ps):right[c,u*d+v]=opcoef(pv,sigma,pu)%modulus if exchange else opcoef(pv,pu,sigma)%modulus
 a=mulmod(left,right,modulus)
 if exchange:a=a.reshape(d,d,d,d).transpose(0,2,1,3).reshape(d*d,d*d)
 return a

def run(h,output_pieces=None):
 assert h>0 and h&(h-1)==0
 start=time.monotonic();ps=partitions(h);d=len(ps);fac=math.factorial(h);v=(fac&-fac).bit_length()-1;power=1<<(3*v);modulus=power*2
 assert modulus<2**31,'This finite NumPy replay needs a smaller weight.'
 chi=characters(h);scaled=np.array([[int(Fraction(fac*c,z(mu)))%modulus for c,mu in zip(row,ps)] for row in chi],dtype=np.int64);outchange=np.array(chi,dtype=np.int64)%modulus
 gpower=1<<(2*v);gmod=2*gpower;gbase=np.array([[opcoef(x,(),y)%gmod for y in ps] for x in ps],dtype=np.int64);gx=scaled%gmod;g=mulmod(mulmod(gx,gbase,gmod),gx.T,gmod)
 assert not np.any(g%gpower);g=(g//gpower)%2;ginv=inverse2(g)
 total=np.zeros((d*d,d*d),dtype=np.int64);records=[];allbits=[];count=0
 for s in range(-2*h,2*h+1):
  a=base(h,s,modulus,exchange=(s>0)).reshape(d,d,d,d)
  for axis in (0,1,3):a=transform(a,scaled,axis,modulus)
  a=transform(a,outchange,2,modulus)
  assert not np.any(a%power),f'homogeneous integrality h{h},s{s}'
  count+=a.size;a=(a//power)%2
  row=a.reshape(d*d,d*d);col=a.transpose(0,2,1,3).reshape(d*d,d*d);r,c=rank2(row),rank2(col)
  total=(total+row)%2;records.append({'s':s,'row':r,'column':c})
  if output_pieces is not None:
   normalized=transform(a,ginv,0,2).reshape(d*d,d*d)
   allbits.append([sum(int(value)<<j for j,value in enumerate(r)) for r in normalized])
  print(json.dumps({'h':h,'s':s,'r':r,'c':c,'seconds':round(time.monotonic()-start,3)}),flush=True)
 expected=np.zeros_like(total)
 for x in range(d):
  for y in range(d):
   for v0 in range(d):expected[x*d+y,y*d+v0]=g[x,v0]
 assert np.array_equal(total,expected)
 normalized=transform(total.reshape(d,d,d,d),ginv,0,2).reshape(d*d,d*d)
 assert all(int(normalized[x*d+y,u*d+v0])==int(x==v0 and y==u) for x in range(d) for y in range(d) for u in range(d) for v0 in range(d))
 if output_pieces is not None:output_pieces.write_text(json.dumps({'h':h,'q':d,'field':2,'shift_minimum':-2*h,'pieces':allbits},indent=2)+'\n')
 S=sum(math.sqrt(r['row']*r['column']) for r in records)
 return {'h':h,'field':2,'dimension':d,'scale':1,'angle_indices':[0,2*h],'gram_rank':rank2(g),'integrality_checks':count,'pieces':records,'rank_sum':S,'rho_rank_sum':math.log(S)/(2*math.log(d)) if d>1 else None,'max_rank_product':max(r['row']*r['column'] for r in records),'seconds':time.monotonic()-start,'operator_cache':str(opcoef.cache_info())}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--weights',type=int,nargs='+',default=[2,4,8]);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--pieces',type=Path);a=ap.parse_args();rs=[]
 for h in a.weights:
  r=run(h,a.pieces);rs.append(r);a.output.write_text(json.dumps({'records':rs,'scope':'Exact finite full-Schur homogeneous crossings in characteristic two; no stronger asymptotic rate claimed.'},indent=2)+'\n');print(json.dumps({'complete':r}),flush=True)
if __name__=='__main__':main()
