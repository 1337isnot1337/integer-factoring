#!/usr/bin/env python3
"""Exact controls for a division-free, truncated-ring Schur constructor."""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
from check_integral_schur_modes import schur as independent_schur

@lru_cache(None)
def partitions(n,maximum=None):
 if n<0:return ()
 if not n:return ((),)
 if maximum is None:maximum=n
 return tuple((k,)+tail for k in range(min(n,maximum),0,-1) for tail in partitions(n-k,k))
@lru_cache(None)
def centralizer(mu):
 v=1
 for k,m in Counter(mu).items():v*=k**m*factorial(m)
 return v

class Ring:
 def __init__(self,weight):self.weight=weight;self.multiplies=0;self.coefficient_products=0;self.max_bits=1
 def inspect(self,p):
  for c in p.values():self.max_bits=max(self.max_bits,abs(c.numerator).bit_length(),c.denominator.bit_length())
  return p
 def add(self,*signed):
  out=defaultdict(Fraction)
  for sign,p in signed:
   for mu,c in p.items():out[mu]+=sign*c
  return self.inspect({mu:c for mu,c in out.items() if c})
 def mul(self,a,b):
  self.multiplies+=1;out=defaultdict(Fraction)
  for mu,c in a.items():
   for nu,d in b.items():
    if sum(mu)+sum(nu)>self.weight:continue
    self.coefficient_products+=1;out[tuple(sorted(mu+nu,reverse=True))]+=c*d
  return self.inspect({mu:c for mu,c in out.items() if c})
 def h(self,n):
  if n<0 or n>self.weight:return {}
  return {mu:Fraction(1,centralizer(mu)) for mu in partitions(n)}

def characteristic(matrix,ring):
 """Coefficients of det(I-tA), using Schur-complement moments only."""
 n=len(matrix)
 if n==0:return [{():Fraction(1)}]
 b=[row[1:] for row in matrix[1:]];previous=characteristic(b,ring)
 moments=[];r=matrix[0][1:];vector=[row[0] for row in matrix[1:]]
 for j in range(n-1):
  moment={}
  for x,y in zip(r,vector):moment=ring.add((1,moment),(1,ring.mul(x,y)))
  moments.append(moment)
  if j<n-2:
   following=[]
   for row in b:
    value={}
    for x,y in zip(row,vector):value=ring.add((1,value),(1,ring.mul(x,y)))
    following.append(value)
   vector=following
 output=[]
 for k in range(n+1):
  value=previous[k] if k<n else {}
  if k>=1:value=ring.add((1,value),(-1,ring.mul(matrix[0][0],previous[k-1])))
  for j in range(k-1):value=ring.add((1,value),(-1,ring.mul(previous[k-2-j],moments[j])))
  output.append(value)
 return output

def schur_construct(lam,ring=None):
 ring=Ring(sum(lam)) if ring is None else ring
 matrix=[[ring.h(lam[i]-i+j) for j in range(len(lam))] for i in range(len(lam))]
 coefficient=characteristic(matrix,ring)[-1]
 return {mu:(-1)**len(lam)*c for mu,c in coefficient.items()}

def verify(maximum):
 counts=Counter();records=[]
 for h in range(maximum+1):
  ps=partitions(h);table=[];max_bits=1;operations=products=0
  for lam in ps:
   ring=Ring(h);actual=schur_construct(lam,ring)
   expected={tuple(sorted(mu,reverse=True)):Fraction(c) for mu,c in independent_schur(tuple(sorted(lam)))}
   assert actual==expected,(h,lam,actual,expected)
   assert all(sum(mu)==h for mu in actual)
   row=[]
   for mu in ps:
    value=actual.get(mu,Fraction())*centralizer(mu);assert value.denominator==1;row.append(int(value));counts['character_integrality']+=1
   table.append(row);counts['schur_polynomials']+=1;max_bits=max(max_bits,ring.max_bits);operations+=ring.multiplies;products+=ring.coefficient_products
  for i,a in enumerate(table):
   for j,b in enumerate(table):assert sum(Fraction(x*y,centralizer(mu)) for x,y,mu in zip(a,b,ps))==int(i==j);counts['hall_pairings']+=1
  records.append({'weight':h,'partition_dimension':len(ps),'ring_multiplications':operations,'coefficient_products':products,'maximum_rational_bits':max_bits})
 # Full ring truncation is needed: a degree>h entry cannot contribute to
 # a determinant of weight h when all weights are nonnegative.
 ring=Ring(2);matrix=[[ring.h(1),ring.h(3)],[ring.h(0),ring.h(1)]]
 assert characteristic(matrix,ring)[2]==ring.mul(ring.h(1),ring.h(1))
 counts['weighted_truncation_control']+=1
 return {'status':'passed','checks':dict(counts),'weights':records,'scope':'Finite independent determinant and lattice controls. Uniform operation and bit bounds are proved in the accompanying construction note.','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-weight',type=int,default=6);ap.add_argument('--output',type=Path);a=ap.parse_args();out=verify(a.max_weight);text=json.dumps(out,indent=2)+'\n';print(text)
 if a.output:a.output.write_text(text)
if __name__=='__main__':main()
