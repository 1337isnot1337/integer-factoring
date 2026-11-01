#!/usr/bin/env python3
"""Exact fixed active core for a deletion-closed physical code and full tail.

The code forbids joint occupation of the two lowest physical modes.
Its alphabet has 3*2**(M-2) letters, but its complete natural profiles
use a rank-144 core of ambient order 256 and Gaussian spectators.
The core is split into occupation blocks of order at most 70.
Every polynomial is constructed completely over F2; numerical pressure
searches are explicitly separate from the exact profile computation.
"""
import argparse,time,json,math
from pathlib import Path
from hashlib import sha256
from check_cancellation_module_automata import pmul
from native_matrix_smith_f4 import MatrixOracle
from math import comb

def addv(a,b):return [x^y for x,y in zip(a,b)]
def unitv(n,i):return [int(j==i) for j in range(n)]
def images(M,k=2):
 assert isinstance(M,int) and M>=2*k and not M&(M-1)
 assert isinstance(k,int) and k>=1 and not k&(k-1)
 n=4*k;x=1<<(2*k);S=sum(1<<(2*k*b) for b in range(M//k-1));T=S<<(2*k)
 Rx=[];Ry=[];Cv=[];Cu=[]
 for a in range(k):
  rcx=unitv(n,a);rax=unitv(n,a);rcy=unitv(n,2*k+a);ray=unitv(n,2*k+a)
  ccv=unitv(n,a);cav=unitv(n,a);ccu=unitv(n,2*k+a);cau=unitv(n,2*k+a)
  for c in range(k):
   fca=3<<(2*(c&~a));fac=3<<(2*(a&~c))
   rcx[2*k+c]^=fca;rcx[3*k+c]^=pmul(fca,T)
   ray[c]^=fac;ray[k+c]^=fac
   cav[2*k+c]^=fca;cav[3*k+c]^=fca<<(2*k)
   ccu[c]^=fac;ccu[k+c]^=pmul(fac,S)
  Rx.append((rcx,rax));Ry.append((rcy,ray));Cv.append((ccv,cav));Cu.append((ccu,cau))
 # Exact canonical anticommutators after each local transport.
 for sides in ((Rx,Ry),(Cv,Cu)):
  ops=sides[0]+sides[1]
  for i,(c,a) in enumerate(ops):
   for j,(d,b) in enumerate(ops):
    value=0
    for u,v in zip(c,b):value^=pmul(u,v)
    assert value==int(i==j)
 return (Rx,Ry,Cv,Cu),[2*i+1 for i in range(2*k)]+[2*M+2*i+1 for i in range(2*k)]

def act(matrix,creation,vector):
 n=len(vector);q=1<<n;out=[[0]*q for _ in range(q)]
 for source,row in enumerate(matrix):
  nonzero=[(i,v) for i,v in enumerate(row) if v]
  if not nonzero:continue
  for mode,coefficient in enumerate(vector):
   if not coefficient or bool(source>>mode&1)==creation:continue
   target=source^(1<<mode);dest=out[target]
   for column,value in nonzero:dest[column]^=pmul(coefficient,value)
 return out
def transpose(matrix):return [list(row) for row in zip(*matrix)]
def projector(matrix,operators,right=False):
 # Pi=I+n0*n1, n_i=c_i*a_i.  On the right transpose reverses the product.
 current=transpose(matrix) if right else matrix
 if right:
  sequence=[(False,operators[0][0]),(True,operators[0][1]),(False,operators[1][0]),(True,operators[1][1])]
 else:
  sequence=[(False,operators[1][1]),(True,operators[1][0]),(False,operators[0][1]),(True,operators[0][0])]
 for creation,vector in sequence:current=act(current,creation,vector)
 if right:current=transpose(current)
 return [[a^b for a,b in zip(row,other)] for row,other in zip(matrix,current)]

def core_matrix(M):
 (rx,ry,cv,cu),orders=images(M)
 q=1<<len(orders)
 matrix=[[0]*q for _ in range(q)]
 for state in range(q):matrix[state][state]=1<<sum(value for i,value in enumerate(orders) if state>>i&1)
 for ops,right in ((rx,False),(ry,False),(cv,True),(cu,True)):
  matrix=projector(matrix,ops,right)
 degree=max(value.bit_length()-1 for row in matrix for value in row)
 assert degree<=16*M+40
 return matrix,degree

def binary_rank(rows):
 pivots={}
 for row in rows:
  while row:
   pivot=row.bit_length()-1
   if pivot in pivots:row^=pivots[pivot]
   else:pivots[pivot]=row;break
 return len(pivots)

def core(M):
 matrix,degree=core_matrix(M)
 q=len(matrix)
 # Rank generating polynomial (1+2z)^2 (1+z)^4.
 ranks=[sum(a*b for i,a in enumerate((1,4,4)) for j,b in enumerate((1,4,6,4,1)) if i+j==charge) for charge in range(9)]
 assert sum(ranks)==144
 result=[];blocks=[]
 for charge,rank in enumerate(ranks):
  indices=[i for i in range(q) if i.bit_count()==charge]
  assert all(matrix[i][j]==0 for i in indices for j in range(q) if j.bit_count()!=charge)
  endpoint=[sum((matrix[i][j].bit_count()%2)<<c for c,j in enumerate(indices)) for i in indices]
  assert binary_rank(endpoint)==rank
  if not rank:
   assert all(not matrix[i][j] for i in indices for j in indices);continue
  block=[[(matrix[i][j],0) for j in indices] for i in indices]
  bits=max(v[0].bit_length() for row in block for v in row)
  precision=128
  while precision<bits:precision*=2
  attempts=[]
  while precision<=8192:
   try:p=MatrixOracle(block,rank,precision=precision).profile()['orders'];attempts.append([precision,True]);break
   except ArithmeticError:
    attempts.append([precision,False]);precision*=2
  else:raise ArithmeticError('Complete active profile exceeded the native precision budget.')
  assert all(x==round(x) and x>=0 for x in p)
  result.extend(map(round,p));blocks.append(dict(charge=charge,size=len(indices),rank=rank,endpoint_rank=rank,
   coefficient_degree=bits-1,precision=precision,precision_attempts=attempts,orders=list(map(round,p))))
 assert len(result)==144
 return sorted(result),blocks,degree

def pressure(M,coreorders,gaussian=False):
 import numpy as np
 from scipy.special import logsumexp
 from scipy.optimize import minimize_scalar,brentq
 orders=np.asarray(coreorders)
 spectator=np.asarray([2*j+1 for j in range(4,M)]+[2*M+2*j+1 for j in range(4,M)])
 physical=np.asarray([2*j+1 for j in range(2,M)])
 def rho(logu):
  u=math.exp(logu)
  z=float(np.logaddexp(0,-u*np.arange(1,4*M,2)).sum()) if gaussian else float(logsumexp(-u*orders)+np.logaddexp(0,-u*spectator).sum())
  def equation(r):
   A=float(np.logaddexp(0,-u*np.arange(1,2*M,2)/r).sum()) if gaussian else float(logsumexp([0,-u/r,-3*u/r])+np.logaddexp(0,-u*physical/r).sum())
   return 2*z-4*r*A
  return brentq(equation,.2,1.5,xtol=1e-13)
 grid=np.linspace(math.log(1e-6),math.log(2),101);values=[rho(x) for x in grid]
 candidates=[(values[0],float(grid[0])),(values[-1],float(grid[-1]))];brackets=[]
 for i in range(1,len(grid)-1):
  if values[i]<=values[i-1] and values[i]<=values[i+1]:
   fit=minimize_scalar(rho,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-10})
   candidates.append((float(fit.fun),float(fit.x)))
   brackets.append(dict(rho=float(fit.fun),u=math.exp(float(fit.x)),success=bool(fit.success)))
 best,logu=min(candidates)
 return dict(rho=best,u=math.exp(logu),u_range=[1e-6,2],grid_points=101,local_minima=brackets,
  scope='Numerical grid and bounded local refinements only; no global or strict rate certificate.')

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--modes',default='4,8,16,32,64,128,256,512')
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 names=(Path(__file__).name,'check_cancellation_module_automata.py','native_matrix_smith_f4.py',
  'native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp')
 report=dict(schema=1,complete=False,factoring_exponent_improved=False,scope=__doc__,records=[],
  source_sha256={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names})
 for M in map(int,args.modes.split(',')):
  start=time.monotonic();orders,blocks,degree=core(M)
  record=dict(M=M,orders=orders,blocks=blocks,physical_alphabet_dimension=3*2**(M-2),
   active_ambient_order=256,active_generic_rank=144,actual_coefficient_degree=degree,
   coefficient_degree_bound=16*M+40,physical_hot_heights=[0,1,3],
   physical_tail_heights=[2*j+1 for j in range(2,M)],
   spectator_orders=[2*j+1 for j in range(4,M)]+[2*M+2*j+1 for j in range(4,M)],
   pressure=pressure(M,orders),gaussian_pressure=pressure(M,orders,True),seconds=time.monotonic()-start)
  report['records'].append(record);args.output.write_text(json.dumps(report,indent=2)+'\n')
  print(dict(M=M,complete_active_orders=len(orders),rho=record['pressure']['rho'],
   gaussian_rho=record['gaussian_pressure']['rho'],seconds=record['seconds']),flush=True)
 report['complete']=True;args.output.write_text(json.dumps(report,indent=2)+'\n')
