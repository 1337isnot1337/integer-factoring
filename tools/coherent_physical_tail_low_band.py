#!/usr/bin/env python3
"""Exact low-band injectivity witnesses for arbitrary constant F4 hot maps.

All 85 projective functionals on the four lowest occupation states are
paired on the two physical legs.  Both ambient-height weighting choices
are retained independently on each leg.  Every successful rank16
specialization certifies rational-function rank16 of the leading band.
"""
import numpy as np, itertools, time,json
P=0x11d

def mul(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&256:a^=P
 return out
MUL=np.array([[mul(a,b) for b in range(256)]for a in range(256)],dtype=np.uint8)
def power(a,n):
 if n<0:return power(power(a,254),-n)
 out=1
 while n:
  if n&1:out=mul(out,a)
  a=mul(a,a);n>>=1
 return out
OMEGA=power(2,85); F4=(0,1,OMEGA,OMEGA^1)
assert mul(OMEGA,OMEGA)==OMEGA^1
assert all(power(a,255)==1 for a in range(1,256))


def images(t,side):
 from coherent_physical_tail_core import finite_field_images
 groups=finite_field_images(2,t,0,mul,power)
 selected=groups[:2] if side=='row' else groups[2:]
 result=[]
 for group in selected:
  for creation,annihilation in group:
   if side=='column':
    creation,annihilation=(annihilation[1],annihilation[0]),(creation[1],creation[0])
   assert all(not any(vector[12:]) for op in (creation,annihilation) for vector in op)
   result.append(tuple(tuple(np.array(vector[:12],dtype=np.uint8) for vector in op)
                       for op in (creation,annihilation)))
 def dot(a,b):return np.bitwise_xor.reduce(MUL[a,b],initial=np.uint8(0))
 def pairing(A,B):return dot(A[0],B[1])^dot(A[1],B[0])
 for i,(c,a) in enumerate(result):
  assert dot(*c)==dot(*a)==0
  for j,(d,b) in enumerate(result):
   assert pairing(c,d)==pairing(a,b)==0
   assert pairing(c,b)==int(i==j)
 return result

INDICES=[]
for i in range(12):
 lo=np.array([s for s in range(4096)if not s>>i&1]);INDICES.append((lo,lo^(1<<i)))
def act(A,op):
 out=np.zeros_like(A)
 for i,(c,a) in enumerate(zip(*op)):
  lo,hi=INDICES[i]
  if c:out[hi]^=MUL[c,A[lo]]
  if a:out[lo]^=MUL[a,A[hi]]
 return out

def rowbasis(A):
 A=A.copy(); pivots=[]; row=0
 for col in range(A.shape[1]):
  nz=np.flatnonzero(A[row:,col])
  if not len(nz):continue
  ix=row+int(nz[0]);A[[row,ix]]=A[[ix,row]];pivots.append(col)
  A[row]=MUL[power(int(A[row,col]),-1),A[row]]
  rows=np.flatnonzero(A[:,col]);rows=rows[rows!=row]
  A[rows]^=MUL[A[rows,col,None],A[row,None,:]]
  row+=1
  if row==A.shape[0]:break
 return pivots

def decomposed(t,side):
 ops=images(t,side);V=np.zeros((4096,16),dtype=np.uint8)
 for i in range(16):V[i,i]=1
 W=[]
 for state in range(16):
  A=V
  for i,(_,a) in enumerate(ops):
   if state>>i&1:A=act(A,a)
  for c,a in ops:A=act(act(A,c),a)
  W.append(A)
 W=np.stack(W)
 nonzero=np.flatnonzero(np.any(W,axis=(0,2)));W=W[:,nonzero,:]
 # Rows of the stacked coordinate evaluation matrix distinguish every image.
 # Select independent coordinate functionals by column pivots on its transpose.
 allv=W.transpose(1,0,2).reshape(len(nonzero),256)
 keep=rowbasis(allv.T)
 W=W[:,keep,:]
 return W,nonzero[np.array(keep)]

def lines():
 out=[]
 for first in range(4):
  for tail in itertools.product(range(4),repeat=3-first):
   out.append((0,)*first+(1,)+tail)
 return out
LINES=lines();assert len(LINES)==85

def rank(A):
 A=A.copy();r=0
 for col in range(A.shape[1]):
  ix=np.flatnonzero(A[r:,col])
  if not len(ix):continue
  j=r+int(ix[0]);A[[r,j]]=A[[j,r]]
  A[r]=MUL[power(int(A[r,col]),-1),A[r]]
  rows=np.flatnonzero(A[r+1:,col])+r+1
  A[rows]^=MUL[A[rows,col,None],A[r,None,:]]
  r+=1
  if r==A.shape[1]:return r
 return r

def screen(t,side,weights=(True,True),decomposition=None):
 start=time.time();W,coordinates=decomposed(t,side) if decomposition is None else decomposition

 ww=np.array([[F4[c]for c in row]for row in LINES],dtype=np.uint8)
 height=np.array([power(t,-h)for h in (0,1,3,4)],dtype=np.uint8)
 first=MUL[ww,height] if weights[0] else ww
 second=MUL[ww,height] if weights[1] else ww
 X=np.zeros((85,4,W.shape[1],16),dtype=np.uint8)
 for i in range(4):
  for j in range(4):X[:,j]^=MUL[first[:,i,None,None],W[i+4*j]]
 bad=[];hist={}
 for p in range(85):
  for q in range(85):
   A=np.bitwise_xor.reduce(MUL[second[q,:,None,None],X[p]],axis=0)
   r=rank(A);hist[r]=hist.get(r,0)+1
   if r<16:bad.append([p,q,r])

 return dict(t=t,side=side,weights=weights,shape=list(W.shape),coordinate_rows=coordinates.tolist(),rank_histogram=hist,failures=bad,seconds=time.time()-start)

if __name__=='__main__':
 import argparse
 from pathlib import Path
 from hashlib import sha256
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 report=dict(schema=1,complete=False,factoring_exponent_improved=False,
  field_polynomial=P,F4_embedding=list(F4),specialization_t=2,
  physical_hot_heights=[0,1,3,4],projective_functionals=[list(x) for x in LINES],records=[],
  source_sha256={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
   for name in (Path(__file__).name,'coherent_physical_tail_core.py','check_cancellation_module_automata.py')},
  scope='Exact finite-field witnesses for the leading T-band of constant F4 physical maps. No all-scale floor or complete high-band profile is claimed.')
 for side in ('row','column'):
  decomposition=decomposed(2,side)
  for weights in ((False,False),(False,True),(True,False),(True,True)):
   result=screen(2,side,weights,decomposition)
   report['records'].append(result)
   args.output.write_text(json.dumps(report,indent=2)+'\n')
   print({k:result[k] for k in ('side','weights','rank_histogram','seconds')},flush=True)
 report['complete']=True
 report['exact_pair_checks']=sum(sum(r['rank_histogram'].values()) for r in report['records'])
 args.output.write_text(json.dumps(report,indent=2)+'\n')
