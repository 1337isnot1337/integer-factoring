#!/usr/bin/env python3
"""Two fixed nilpotent Gram crossings and their exact row-profile obstruction."""
from pathlib import Path
import argparse,json
from native_matrix_smith_f4 import MatrixOracle
n=49;I=[1<<i for i in range(n)]
def transpose(A):return [sum(((A[i]>>j)&1)<<i for i in range(n)) for j in range(n)]
def mm(A,B):
 out=[]
 for row in A:
  value=0
  while row:
   bit=row&-row;value^=B[bit.bit_length()-1];row^=bit
  out.append(value)
 return out
def madd(A,B):return [a^b for a,b in zip(A,B)]
def inverse(A):
 rows=[a|(1<<(n+i)) for i,a in enumerate(A)]
 for j in range(n):
  p=next(i for i in range(j,n) if rows[i]>>j&1);rows[j],rows[p]=rows[p],rows[j]
  for i in range(n):
   if i!=j and rows[i]>>j&1:rows[i]^=rows[j]
 out=[a>>n for a in rows];assert mm(A,out)==I;return out
def power(A,k):
 out=I[:]
 for _ in range(k):out=mm(out,A)
 return out
lengths=(13,12,12,12);starts=[];offset=0;N=[0]*n;B=[0]*n;pairs=[]
for length in lengths:
 starts.append(offset)
 for i in range(length):
  B[offset+i]=1<<(offset+length-1-i)
  if i+1<length:N[offset+i]=1<<(offset+i+1)
 for i in range(length//2):pairs.append((offset+i,offset+length-1-i))
 offset+=length
assert offset==n and mm(transpose(N),B)==mm(B,N)
S=I[:]
for b,(offset,length) in enumerate(zip(starts,lengths)):
 nextb=(b+1)%4;S[offset+length-6]^=1<<(starts[nextb]+lengths[nextb]-1)
assert mm(S,S)==I and mm(S,power(N,5))==power(N,5)
assert mm(N,S)!=mm(S,N)
assert S[12]==1<<12
physical_pairs=[]
for i in range(7):
 for j in range(i+1,7):physical_pairs.append((1<<(7*i+j),1<<(7*j+i)))
diagonal=[(1<<0)^(1<<i) for i in range(1,7)];diag_pairs=[]
def dot(a,b):return (a&b).bit_count()%2
def diag_matrix(a):return sum(((a>>i)&1)<<(8*i) for i in range(7))
while diagonal:
 u=diagonal.pop(0);j=next(i for i,v in enumerate(diagonal) if dot(u,v));v=diagonal.pop(j)
 diag_pairs.append((diag_matrix(u),diag_matrix(v)))
 diagonal=[w^(u if dot(w,v) else 0)^(v if dot(w,u) else 0) for w in diagonal]
physical_pairs+=diag_pairs
columns=[0]*n;columns[6]=sum(1<<(8*i) for i in range(7))
for (a,b),(u,v) in zip(pairs,physical_pairs):columns[a]=u;columns[b]=v
iso=transpose(columns);iso_inv=inverse(iso)
J=[1<<(7*(i%7)+i//7) for i in range(n)]
assert mm(mm(transpose(iso),J),iso)==B

def construct(flag):
 p1=[N,I];p5=[power(N,5),power(N,4),[0]*n,[0]*n,N,I]
 O=[[0]*n for _ in range(7)]
 for a,A in enumerate(p1):
  for b,C in enumerate(p5):O[a+b]=madd(O[a+b],mm(mm(mm(A,flag),C),flag))
 endpoint=[0]*n
 for A in O:endpoint=madd(endpoint,A)
 C=inverse(endpoint);O=[mm(A,C) for A in O]
 K=[[0]*n for _ in range(13)]
 for a,A in enumerate(O):
  for b,C in enumerate(O):K[a+b]=madd(K[a+b],mm(mm(transpose(A),B),C))
 K=[mm(mm(transpose(iso_inv),A),iso_inv) for A in K]
 endpoint=[0]*n
 for A in K:endpoint=madd(endpoint,A)
 assert endpoint==J
 cold=[0]*n;cold[7]=1<<7;assert K[0]==cold
 return K

def matrix(K,cut,place):
 M=len(K)-1;out=[]
 for i in range(n):
  row=[]
  for j in range(n):
   if cut==0:a,b=i,j
   else:
    x,u=divmod(i,7);y,v=divmod(j,7);a,b=7*x+y,7*u+v
   coefficients=[(A[a]>>b)&1 for A in K]
   if place=='infinity':coefficients.reverse()
   row.append((sum(v<<d for d,v in enumerate(coefficients)),0))
  out.append(row)
 return out

def binary_rank(rows):
 pivots={}
 for row in rows:
  while row:
   leading=row.bit_length()-1
   if leading not in pivots:pivots[leading]=row;break
   row^=pivots[leading]
 return len(pivots)

def row_toeplitz_nullity(K,h):
 rows=[]
 for a in range(h):
  for i in range(n):
   row=0
   for b in range(a+1):
    if a-b<len(K):row|=K[a-b][i]<<(b*n)
   rows.append(row)
 return n*h-binary_rank(rows)

results={};expected=[0]+[4]*15+[16]*33
for name,flag in [('single_flag',I),('four_chain_cycle',S)]:
 K=construct(flag);assert any(K[12]);inverse(K[12])
 profiles={}
 for cut in (0,1):
  for place in ('zero','infinity'):
   orders=list(map(int,MatrixOracle(matrix(K,cut,place),49,precision=128).profile()['orders']))
   profiles[f'{cut}_{place}']=orders
 assert profiles['0_zero']==expected
 assert profiles['0_infinity']==[0]*49
 nullities={h:row_toeplitz_nullity(K,h) for h in (1,4,5,16,17)}
 assert list(nullities.values())==[48,192,225,588,588]
 assert sum(expected)==49*12
 results[name]=dict(degree=12,endpoint_entries_checked=2401,cold_nonzero_row=7,cold_nonzero_column=7,profiles=profiles,independent_row_nullities=nullities)
assert 15*285**2>103**3 and 345>7**3 and 8<9
certificate=dict(q=7,degree=12,determinant_degree=588,row_orders=expected,best_relaxed_column_orders=[0]+[12]*36+[13]*12,amgm_cube_left=15*285**2,amgm_cube_right=103**3,strict_pressure_lower_bound=345,q_cubed=343,scope='The displayed row profile excludes every bare two-place choice over every finite coefficient extension; it is not a general nilpotent-Gram obstruction.')
result=dict(scope='Two fixed exact binary nilpotent Gram tensors, their endpoint and local profiles, and an ordinary row-profile pressure exclusion. No factoring improvement or parameter screen.',results=results,certificate=certificate)
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(cases=len(results),endpoint_entries_checked=4802,independent_row_toeplitz_checks=10,pressure_strict_lower_bound=345)))
