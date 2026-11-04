#!/usr/bin/env python3
"""One exact confluent q7 quantum-spin tensor; no factoring improvement.

Regularizes the projector formulas in quantum_crossing_initial_forms.py
by paired Hermite conditions at a seventh root in F8.
"""
from pathlib import Path
import argparse,json
from functools import reduce
from operator import xor
from binary_extension_smith import BinaryExtension,ExtensionMatrixOracle
from check_binary_extension_smith import rank
F=BinaryExtension(11);mul=F.multiply;inv=F.inverse
T=[[mul(a,b) for b in range(8)] for a in range(8)]
assert F.power(2,7)==1 and 2!=1
def fp(a,n):return F.power(a,n) if n>=0 else F.power(inv(a),-n)
def add(a,b):return tuple(x^y for x,y in zip(a,b))
def jm(a,b):return tuple(reduce(xor,(mul(a[j],b[i-j]) for j in range(i+1)),0) for i in range(3))
def ji(a):
 c=inv(a[0]);b=[c,0,0]
 for k in range(1,3):b[k]=mul(c,reduce(xor,(mul(a[j],b[k-j]) for j in range(1,k+1)),0))
 assert jm(a,b)==(1,0,0)
 return tuple(b)
def jp(a,n):
 if n<0:return jp(ji(a),-n)
 z=(1,0,0)
 for _ in range(n):z=jm(z,a)
 return z
Q=(2,1,0)
def qnum(n):
 n=abs(n);out=(0,0,0)
 for j in range(n):out=add(out,jp(Q,n-1-2*j))
 return out
lam=[jm(qnum(j),qnum(j+1)) for j in range(7)]
def fpoly(j):
 out=[(1,0,0)]
 for ell in range(1,7):
  power=jp(Q,2*ell);a=add(power,(1,0,0));b=(1,0,0) if ell<=j else power
  nxt=[(0,0,0)]*(len(out)+1)
  for d,c in enumerate(out):nxt[d]=add(nxt[d],jm(c,a));nxt[d+1]=add(nxt[d+1],jm(c,b))
  out=nxt
 return out
fs=[fpoly(j) for j in range(7)];D=(1,0,0)
for ell in range(1,7):D=jm(D,add(jp(Q,2*ell),(1,0,0)))
assert D[0]==1
conditions=[];V=[];pairs=[]
for j in range(3):
 other=6-j;diff=add(lam[j],lam[other]);assert diff[:2]==(0,0) and diff[2]
 assert all(add(a,b)[:2]==(0,0) for a,b in zip(fs[j],fs[other]))
 value=lam[j][0]
 V.append([fp(value,a) for a in range(7)])
 conditions.append([mul(c[0],inv(D[0])) for c in fs[j]])
 V.append([fp(value,a-1) if a%2 else 0 for a in range(7)])
 conditions.append([mul(a[2]^b[2],inv(mul(D[0],diff[2]))) for a,b in zip(fs[j],fs[other])])
 pairs.append(dict(spins=[j,other],eigenvalue=value,collision_order=2,leading_difference=diff[2]))
V.append([fp(lam[3][0],a) for a in range(7)]);conditions.append([mul(c[0],inv(D[0])) for c in fs[3]])
assert len({l[0] for l in lam})==4
A=[row[:]+[int(i==j) for j in range(7)] for i,row in enumerate(V)]
for j in range(7):
 pivot=next(i for i in range(j,7) if A[i][j]);A[pivot],A[j]=A[j],A[pivot];c=inv(A[j][j]);A[j]=[mul(c,a) for a in A[j]]
 for i in range(7):
  if i!=j:
   c=A[i][j]
   if c:A[i]=[a^mul(c,b) for a,b in zip(A[i],A[j])]
Vinv=[row[7:] for row in A]
rcoef=[[reduce(xor,(mul(Vinv[a][i],conditions[i][d]) for i in range(7)),0) for a in range(7)] for d in range(7)]
assert rcoef[0]==[1,0,0,0,0,0,0]
for d in range(7):
 for i in range(7):
  assert reduce(xor,(mul(V[i][a],rcoef[d][a]) for a in range(7)),0)==conditions[i][d]
def mm(a,b):return [[reduce(xor,(mul(x,y) for x,y in zip(row,col)),0) for col in zip(*b)] for row in a]
def meval(coeff,powers):return [[reduce(xor,(mul(c,powers[a][i][j]) for a,c in enumerate(coeff)),0) for j in range(len(powers[0]))] for i in range(len(powers[0]))]
values={}
for w in range(13):
 basis=[(i,w-i) for i in range(7) if 0<=w-i<7];n=len(basis);index={b:i for i,b in enumerate(basis)};C=[[0]*n for _ in range(n)]
 for col,(i,j) in enumerate(basis):
  C[col][col]^=mul(qnum(6-w)[0],qnum(7-w)[0])
  raises=[((i-1,j),mul(qnum(i)[0],fp(2,-3+j))),((i,j-1),mul(qnum(j)[0],fp(2,3-i)))]
  for (ii,jj),e in raises:
   if not e:continue
   lowers=[((ii+1,jj),mul(qnum(6-ii)[0],fp(2,-3+jj))),((ii,jj+1),mul(qnum(6-jj)[0],fp(2,3-ii)))]
   for dest,f in lowers:
    if f and dest in index:C[index[dest]][col]^=mul(e,f)
 powers=[[[int(i==j) for j in range(n)] for i in range(n)]]
 for a in range(1,7):powers.append(mm(powers[-1],C))
 coeff=[meval(c,powers) for c in rcoef]
 for i,(x,y) in enumerate(basis):
  for j,(u,v) in enumerate(basis):
   entry=[coeff[d][i][j] for d in range(7)]
   while entry and not entry[-1]:entry.pop()
   if entry:values[y,x,u,v]=entry
for x in range(7):
 for y in range(7):
  for u in range(7):
   for v in range(7):assert (values.get((x,y,u,v),[0])[0])==int(x==v and y==u)
assert max(len(p)-1 for p in values.values())==6
from fused_rmatrix_probe import cut_blocks
ranks={}
for point in range(8):
 ranks[point]=[]
 for cut in (0,1):
  total=0
  for block in cut_blocks(values,7,cut):
   evaluated=[[reduce(xor,(mul(c,fp(point,d)) for d,c in enumerate(entry)),0) for entry in row] for row in block]
   total+=rank(evaluated,T)
  ranks[point].append(total)

def changed(p,point):
 if point=='infinity':return (p+[0]*(7-len(p)))[::-1]
 out=[0]*len(p)
 for d,c in enumerate(p):
  for j in range(d+1):
   if j&d==j:out[j]^=mul(c,F.power(point,d-j))
 return out

def toeplitz_nullity(matrix,h):
 n=len(matrix);A=[[0]*(n*h) for _ in range(n*h)]
 for a in range(h):
  for b in range(a+1):
   for i in range(n):
    for j in range(n):
     p=matrix[i][j]
     if a-b<len(p):A[a*n+i][b*n+j]=p[a-b]
 return n*h-rank(A,T)

profiles={};totals={};checks=0
for cut in (0,1):
 profiles[str(cut)]={};totals[str(cut)]={}
 for point in list(range(8))+['infinity']:
  orders=[]
  for block in cut_blocks(values,7,cut):
   matrix=[[changed(p,point) for p in row] for row in block]
   local=ExtensionMatrixOracle(matrix,11,len(block)).profile()['orders']
   for h in range(1,max(local)+2):
    assert toeplitz_nullity(matrix,h)==sum(min(h,v) for v in local)
    checks+=1
   orders+=local
  if point!='infinity':assert sum(v==0 for v in orders)==ranks[point][cut]
  profiles[str(cut)][str(point)]=sorted(orders)
  totals[str(cut)][str(point)]=sum(orders)
assert sum(totals['0'].values())==294
assert sum(totals['1'].values())==190
assert max(totals['0'].values())==84
assert max(totals['1'].values())==66
# Each non-F8 point has at least two distinct Frobenius conjugates.
remaining=294-sum(totals['1'].values())
assert remaining//2==52<66
assert 84+66<294
result=dict(scope='One degree-six confluent q7 tensor over F8. Determinant budgets exclude its bare two-place criterion over every finite coefficient extension; no general crossing obstruction or factoring gain.',field_modulus=11,root=2,degree=6,eigenvalue_jets=lam,collisions=pairs,normalization_jet=D,hermite_coefficients=rcoef,endpoint_entries_checked=2401,point_ranks=ranks,profiles=profiles,determinant_totals=totals,independent_toeplitz_checks=checks,extension_certificate=dict(homogeneous_determinant_degree=294,first_cut_known_total=294,first_cut_maximum=84,second_cut_known_total=190,second_cut_remaining=remaining,second_cut_nonrational_maximum=52,second_cut_maximum=66,two_cut_maximum=150,positive_jensen_exponent_numerator=144,pressure_lower_bound=2401))
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('degree','endpoint_entries_checked','independent_toeplitz_checks','extension_certificate')}))
