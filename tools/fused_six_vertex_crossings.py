#!/usr/bin/env python3
"""Literal finite-field fusion of checked six-vertex polynomial matrices.

Each requested instance is constructed in the full qubit space, checked for
invariance of the selected symmetric subspace, divided by its exact common
polynomial factor, and checked again for a nonzero scalar identity endpoint.
There is no unverified all-spin endpoint or rank assumption.
"""
from f16_identity_operator_profiles import MUL,INV,add,multiply
from itertools import product
import json,time

def power(a,n):
 if n<0:return power(INV[a],-n)
 out=1
 while n:
  if n&1:out=MUL[out][a]
  a=MUL[a][a];n//=2
 return out
def pscale(p,c):return {k:MUL[c][v] for k,v in p.items() if MUL[c][v]}
def canon(m,k):return ((1<<k)-1)<<(m-k)
def sym(m,theta,k):
 result={}
 for word in range(1<<m):
  if word.bit_count()!=k:continue
  inv=sum(((word>>i)&1)*(1-((word>>j)&1)) for i in range(m) for j in range(i+1,m))
  result[word]=power(theta,-inv)
 return result
def gate(vector,position,theta,factor):
 out={};c=theta^INV[theta]
 for word,polynomial in vector.items():
  a=(word>>position)&1;b=(word>>(position+1))&1
  if a==b:terms=[(word,{0:INV[theta],1:MUL[theta][factor]})]
  elif a==0:terms=[(word,{1:MUL[c][factor]}),(word^(3<<position),{0:1,1:factor})]
  else:terms=[(word,{0:c}),(word^(3<<position),{0:1,1:factor})]
  for target,coefficient in terms:
   value=add(out.get(target,{}),multiply(polynomial,coefficient))
   if value:out[target]=value
   elif target in out:del out[target]
 return out
def pdiv(a,b):
 a=dict(a);result={}
 while a and max(a)>=max(b):
  shift=max(a)-max(b);scalar=MUL[a[max(a)]][INV[b[max(b)]]];result[shift]=scalar
  a=add(a,{degree+shift:MUL[scalar][value] for degree,value in b.items()})
 return result,a
def gcd(a,b):
 while b:a,b=b,pdiv(a,b)[1]
 return pscale(a,INV[a[max(a)]]) if a else a
def construct(m,theta):
 size=m+1;embedding=[sym(m,theta,k) for k in range(size)]
 result=[[{} for _ in range(size*size)] for _ in range(size*size)];invariant=True
 for a,b in product(range(size),repeat=2):
  vector={x|(y<<m):{0:MUL[c][d]} for x,c in embedding[a].items() for y,d in embedding[b].items()}
  for i in range(m-1,-1,-1):
   for j in range(m):vector=gate(vector,i+j,theta,power(theta,2*(i-j)))
  reconstructed={}
  for c,d in product(range(size),repeat=2):
   value=vector.get(canon(m,c)|(canon(m,d)<<m),{})
   result[c*size+d][a*size+b]=value
   for x,v in embedding[c].items():
    for y,w in embedding[d].items():
     target=x|(y<<m);term=pscale(value,MUL[v][w])
     if term:reconstructed[target]=term
  if reconstructed!=vector:invariant=False
 common={}
 for row in result:
  for value in row:
   if value:common=gcd(common,value)
 if not common:return None
 for row in result:
  for i,value in enumerate(row):
   row[i],rem=pdiv(value,common);assert not rem
 endpoints=[[0 for _ in row] for row in result]
 for i,row in enumerate(result):
  for j,value in enumerate(row):
   for x in value.values():endpoints[i][j]^=x
 scalar=endpoints[0][0]
 endpoint_ok=bool(scalar) and all(endpoints[i][j]==(scalar if i==j else 0) for i in range(size*size) for j in range(size*size))
 if endpoint_ok:result=[[pscale(value,INV[scalar]) for value in row] for row in result]
 return {'m':m,'theta':theta,'invariant':invariant,'common_degree':max(common),'degree':max(max(value,default=0) for row in result for value in row),'endpoint_identity':endpoint_ok,'matrix':result}

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from native_matrix_smith_f4 import MatrixOracle
from f4_lattice_linear import MUL as MUL4

def transform(matrix,degree,a,b):
 if b is None:num={0:a,1:1^a};den={0:1}
 elif a is None:num={0:1^b,1:b};den={1:1}
 else:
  d=MUL[1^a][INV[1^b]];num={0:a,1:MUL[b][d]};den={0:1,1:d}
 num={k:v for k,v in num.items() if v};den={k:v for k,v in den.items() if v}
 powers_num=[{0:1}];powers_den=[{0:1}]
 for _ in range(degree):powers_num.append(multiply(powers_num[-1],num));powers_den.append(multiply(powers_den[-1],den))
 factors=[multiply(powers_num[k],powers_den[degree-k]) for k in range(degree+1)]
 atone=0
 for v in den.values():atone^=v
 normal=INV[power(atone,degree)]
 out=[]
 for row in matrix:
  new=[]
  for value in row:
   result={}
   for k,c in value.items():result=add(result,factors[k],c)
   new.append(pscale(result,normal))
  out.append(new)
 return out
def expanded(matrix):
 size=len(matrix);out=[[(0,0)]*(2*size) for _ in range(2*size)]
 for i,row in enumerate(matrix):
  for j,p in enumerate(row):
   entries=[[{},{}],[{},{}]]
   for k,v in p.items():
    a,b=v&3,v>>2
    coeff=((a,MUL4[2][b]),(b,a^b))
    for r,c in product(range(2),repeat=2):
     if coeff[r][c]:entries[r][c][k]=coeff[r][c]
   for r,c in product(range(2),repeat=2):
    out[2*i+r][2*j+c]=(sum((v&1)<<k for k,v in entries[r][c].items()),sum((v>>1)<<k for k,v in entries[r][c].items()))
 return out
def oracles(matrix,q,precision=512):
 cuts=[]
 for side in range(2):
  cut=[[{} for _ in range(q*q)] for _ in range(q*q)]
  for x,y,u,v in product(range(q),repeat=4):
   row,col=(q*x+y,q*u+v) if side==0 else(q*x+u,q*y+v)
   cut[row][col]=matrix[q*y+x][q*u+v]
  cuts.append(MatrixOracle(expanded(cut),2*q*q,precision=precision))
 return cuts
def evaluate(matrix,q):
 out=[]
 for side,oracle in enumerate(oracles(matrix,q)):
  p=oracle.profile(infinity=bool(side));orders=p['orders']
  assert all(abs(orders[i]-orders[i+1])<1e-9 for i in range(0,len(orders),2))
  out.append(tuple(round(v) for v in orders[::2]))
 a,b=map(np.array,out)
 f=lambda t:float(logsumexp(-t*a)+logsumexp(-t*b))
 best=minimize_scalar(f,bounds=(1e-8,12),method='bounded',options={'xatol':1e-10})
 return {'rho':best.fun/(4*np.log(q)),'scale':float(best.x),'profiles':out}
