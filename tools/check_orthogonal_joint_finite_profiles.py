#!/usr/bin/env python3
"""Complete F4/F8/F16 controls; not an all-field profile theorem."""
import json,time,itertools
from pathlib import Path
import numpy as np
from joint_energy_conjugations import model,kernel,binary_product
from binary_extension_smith import BinaryExtension,ExtensionMatrixOracle
D=model();L=D['leading'].basis;squares=[binary_product(a,a) for a in L];K=kernel([sum(((x>>i)&1)<<j for j,x in enumerate(squares)) for i in range(625)],10);rad=[]
for v in K:
 A=0
 for i,a in enumerate(L):
  if v>>i&1:A^=a
 rad.append(A)
assert K==[6,24,64,384]
P,Q=[np.array([[int(a>>(25*i+j)&1) for j in range(25)] for i in range(25)],dtype=np.uint8) for a in rad[:2]]
coeff=np.zeros((19,25,25),dtype=np.uint8)
for k,a in D['weighted'].items():
 for i in range(25):
  for j in range(25):coeff[k+8,i,j]=(a>>(25*i+j))&1
ref=json.load(open('evidence/joint-energy-conjugations.json'))['baseline_infinity']
for mod in (7,11,19):
 f=BinaryExtension(mod);M=np.array([[f.multiply(a,b) for b in range(f.size)] for a in range(f.size)],dtype=np.uint8);patterns={};records=[];start=time.monotonic()
 for a,b in itertools.product(range(f.size),repeat=2):
  U=np.eye(25,dtype=np.uint8)^M[a,P]^M[b,Q]
  left=np.bitwise_xor.reduce(M[U[None,:,:,None],coeff[:,None,:,:]],axis=2)
  out=np.bitwise_xor.reduce(M[left[:,:,:,None],U[None,None,:,:]],axis=2)
  C=np.zeros_like(out)
  for x,y,u,v in itertools.product(range(5),repeat=4):C[:,5*x+u,5*y+v]=out[:,x+5*y,v+5*u]
  A=[[C[:,i,j].tolist() for j in range(25)] for i in range(25)]
  orders=tuple(x-8 for x in ExtensionMatrixOracle(A,mod,25).profile()['orders'])
  patterns.setdefault(orders,[]).append([a,b]);cum1=cum2=0;dom=True
  for x,y in zip(orders,ref):cum1+=x;cum2+=y;dom&=cum1<=cum2
  assert dom,(mod,a,b,orders)
  records.append({'parameters':[a,b],'profile_index':list(patterns).index(orders)})
 result={'completed':True,'field_polynomial':mod,'cases':f.size**2,'profile_classes':[{'orders':list(k),'parameters':v} for k,v in patterns.items()],'all_prefix_dominated_by_reference':True,'seconds':time.monotonic()-start,'scope':'Complete finite-field two-parameter orthogonal orbit; not an all-field pressure theorem.'}
 Path('evidence/orthogonal-two-parameter-'+str(f.size)+'.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'field':f.size,'cases':f.size**2,'profile_classes':len(patterns),'seconds':result['seconds']}),flush=True)
