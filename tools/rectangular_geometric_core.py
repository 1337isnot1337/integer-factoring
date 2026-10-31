"""Closed geometric concentration, specialized at z=c with Z formal.

Every sequence has a finite sum of z^(e*j)(1+z)^(f*j) terms and
finitely many endpoint deltas. M is a power of two, so the full-range
sum is computed using r(z)^M=r(Z), never by an M-dependent loop.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache

def mul2(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length()-1)
        b ^= bit
    return result


def mul4(a, b):
    first, second = mul2(a[0], b[0]), mul2(a[1], b[1])
    cross = mul2(a[0]^a[1], b[0]^b[1])
    return first^second, cross^first

MUL=((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))
INV=(0,1,3,2)


def degree(a):return max(a).bit_length()-1
def scalar(a,c):
 if c==0:return(0,0)
 if c==1:return a
 if c==2:return(a[1],a[0]^a[1])
 return(a[0]^a[1],a[0])
def lead(a):
 k=degree(a)
 return((a[0]>>k)&1)+2*((a[1]>>k)&1)if k>=0 else 0
def pdiv(a,b):
 if not any(b):raise ZeroDivisionError
 q=(0,0);a=tuple(a);db=degree(b);ib=INV[lead(b)]
 while any(a)and degree(a)>=db:
  shift=degree(a)-db;c=MUL[lead(a)][ib];s=scalar(b,c)
  a=(a[0]^(s[0]<<shift),a[1]^(s[1]<<shift))
  q=(q[0]^((c&1)<<shift),q[1]^(((c>>1)&1)<<shift))
 return q,a
def gcd(a,b):
 while any(b):a,b=b,pdiv(a,b)[1]
 return scalar(a,INV[lead(a)])if any(a)else(0,0)


@dataclass(frozen=True)
class Rat:
 n:tuple
 d:tuple=(1,0)
 def __bool__(self):return any(self.n)
 def __add__(self,b):return radd(self,rat(b))
 __radd__=__add__
 def __mul__(self,b):return rmul(self,rat(b))
 __rmul__=__mul__
 def __truediv__(self,b):return rdiv(self,rat(b))
 def __rtruediv__(self,b):return rdiv(rat(b),self)
 def __pow__(self,k):
  if k<0:return normalize(self.d,self.n)**(-k)
  result=ONE;base=self
  while k:
   if k&1:result=result*base
   base=base*base;k>>=1
  return result
 def order(self):
  if not self:return None
  return ((self.n[0]|self.n[1])&-(self.n[0]|self.n[1])).bit_length()-((self.d[0]|self.d[1])&-(self.d[0]|self.d[1])).bit_length()


def rat(a):
 if isinstance(a,Rat):return a
 if isinstance(a,int)and 0<=a<4:return Rat((a&1,a>>1))
 raise TypeError(a)
@lru_cache(maxsize=100000)
def normalize(n,d):
 if not any(d):raise ZeroDivisionError
 if not any(n):return Rat((0,0))
 common=gcd(n,d);n=pdiv(n,common)[0];d=pdiv(d,common)[0];unit=INV[lead(d)]
 return Rat(scalar(n,unit),scalar(d,unit))
@lru_cache(maxsize=200000)
def radd(a,b):
 if not a:return b
 if not b:return a
 x=mul4(a.n,b.d);y=mul4(b.n,a.d)
 return normalize((x[0]^y[0],x[1]^y[1]),mul4(a.d,b.d))
@lru_cache(maxsize=200000)
def rmul(a,b):
 if not a or not b:return ZERO
 return normalize(mul4(a.n,b.n),mul4(a.d,b.d))
@lru_cache(maxsize=100000)
def rdiv(a,b):
 if not b:raise ZeroDivisionError
 return normalize(mul4(a.n,b.d),mul4(a.d,b.n))


ZERO,ONE=rat(0),rat(1)
Z=Rat((2,0));UZ=ONE+Z
Z0=2
GEOMETRIC_RATIOS=set()
PIVOT_CHECKS=[]


def fpow(a,k):
 if not a and k<0:raise ZeroDivisionError
 if k<0:a=INV[a];k=-k
 result=1
 while k:
  if k&1:result=MUL[result][a]
  a=MUL[a][a];k>>=1
 return result
def base_value(base):return MUL[fpow(Z0,base[0])][fpow(1^Z0,base[1])]
def base_macro(base):return Z**base[0]*UZ**base[1]
def mpow(base,index):return base_macro(base)**index[0]*rat(fpow(base_value(base),index[1]))
def mono(a=0,b=0,c=0,d=0):return Z**a*UZ**c*rat(MUL[fpow(Z0,b)][fpow(1^Z0,d)])


class Sequence:
 def __init__(self,terms=None,deltas=None):
  self.terms={k:rat(v)for k,v in (terms or {}).items()if rat(v)}
  self.deltas={k:rat(v)for k,v in (deltas or {}).items()if rat(v)}
 def __add__(self,other):
  terms=dict(self.terms);deltas=dict(self.deltas)
  for k,v in other.terms.items():terms[k]=terms.get(k,ZERO)+v
  for k,v in other.deltas.items():deltas[k]=deltas.get(k,ZERO)+v
  return Sequence(terms,deltas)
 def scale(self,value):return Sequence({k:v*value for k,v in self.terms.items()},{k:v*value for k,v in self.deltas.items()})
 def base_at(self,index):return sum((v*mpow(k,index)for k,v in self.terms.items()),ZERO)
 def at(self,index):return self.base_at(index)+self.deltas.get(index,ZERO)
 def __mul__(self,other):
  terms={}
  for a,x in self.terms.items():
   for b,y in other.terms.items():
    k=(a[0]+b[0],a[1]+b[1]);terms[k]=terms.get(k,ZERO)+x*y
  deltas={i:self.at(i)*other.at(i)+self.base_at(i)*other.base_at(i)for i in self.deltas.keys()|other.deltas.keys()}
  return Sequence(terms,deltas)
 def total(self):
  answer=sum(self.deltas.values(),ZERO)
  for base,value in self.terms.items():
   if base==(0,0):continue # M is even.
   constant=base_value(base)
   GEOMETRIC_RATIOS.add((base,constant))
   if constant==1:raise ArithmeticError(('bad specialization for geometric ratio',base,Z0))
   answer=answer+value*(ONE+base_macro(base))/rat(1^constant)
  return answer
 def outside_sum(self,excluded):return self.total()+sum((self.at(i)for i in excluded),ZERO)
 def delta_add(self,index,value):return self+Sequence(deltas={index:value})
 def outside_add(self,value,excluded):
  answer=self+value
  for index in excluded:answer=answer.delta_add(index,value.at(index))
  return answer


def seq(base=(0,0),value=ONE):return Sequence({base:value})
def coord(index,value=ONE):return Sequence(deltas={index:value})
LOW=(0,0);ONEINDEX=(0,1);TWO=(0,2);TOP=(1,-1);HIGH=(1,-2);HIGHB=(1,-3)


def many_rows(images,block,pivot,coefficient,excluded):
 for name,value in images.items():
  creation,annihilation=value[block],value[block+2]
  cp=creation.at(pivot);change=(annihilation*coefficient).outside_sum(excluded)
  value[block]=creation.outside_add(coefficient.scale(cp),excluded)
  value[block+2]=annihilation.delta_add(pivot,change)


def many_columns(images,block,pivot,coefficient,excluded):
 for name,value in images.items():
  creation,annihilation=value[block],value[block+2]
  ap=annihilation.at(pivot);change=(creation*coefficient).outside_sum(excluded)
  value[block]=creation.delta_add(pivot,change)
  value[block+2]=annihilation.outside_add(coefficient.scale(ap),excluded)


def concentrate_four(images,block,s,kind):
 kept={LOW,ONEINDEX,HIGH,TOP};first_excluded={LOW,ONEINDEX,TOP}
 f1=seq(s,mpow(s,(0,-1)))
 many_rows(images,block,ONEINDEX,f1,first_excluded)
 f2=(seq()+f1).scale(ONE/(ONE+mpow(s,(1,-3))))
 many_rows(images,block,HIGH,f2,kept)
 if kind=='constant':
  pivot=ONE+f1.outside_sum(first_excluded)
  expected=(ONE+mpow(s,(1,-2)))/(ONE+mpow(s,(0,1)))
  assert pivot==expected
  PIVOT_CHECKS.append({'kind':'constant','ratio':s,'order':pivot.order()})
  f3=seq(value=ONE/pivot);target=ONEINDEX
 else:
  inverse=seq((-s[0],-s[1]))
  pivot=inverse.at(HIGH)+(f2*inverse).outside_sum(kept)
  expected=mpow(s,(-1,2))*(ONE+mpow(s,(1,-2)))/((ONE+mpow(s,(0,1)))*(ONE+mpow(s,(1,-3))))
  assert pivot==expected
  PIVOT_CHECKS.append({'kind':'inverse','ratio':s,'order':pivot.order()})
  f3=inverse.scale(ONE/pivot);target=HIGH
 many_columns(images,block,target,f3,kept)
 return sorted(kept)


def concentrate_six(images,block,s):
 kept={LOW,ONEINDEX,TWO,HIGHB,HIGH,TOP}
 f1=seq(s,mpow(s,(0,-1)))
 many_rows(images,block,ONEINDEX,f1,kept-{HIGH})
 f2=(seq()+f1).scale(ONE/(ONE+mpow(s,(1,-3))))
 many_rows(images,block,HIGH,f2,kept)
 f3=seq((-s[0],-s[1]),mpow(s,HIGHB))
 many_columns(images,block,HIGHB,f3,kept-{TWO})
 f4=(seq()+f3).scale(ONE/(ONE+mpow(s,(1,-5))))
 many_columns(images,block,TWO,f4,kept)
 return sorted(kept)


def initial_images(place):
 images={}
 def put(name,c1=None,c2=None,a1=None,a2=None):images[name]=[v if v is not None else Sequence()for v in(c1,c2,a1,a2)]
 for index,label in ((LOW,'0'),(TOP,'M-1')):
  high=index==TOP
  if place=='infinity':
   put('out_v'+label+'_create',c2=coord(index,mono(d=-1)))
   put('out_y'+label+'_create',c1=coord(index),c2=seq((0,2))if high else coord(LOW))
   put('out_v'+label+'_annihilate',a1=coord(TOP,mono(c=2,d=-1))if high else seq(value=mono(d=1)),a2=coord(index,mono(d=1)))
   put('out_y'+label+'_annihilate',a1=coord(index))
   sigma=mono(a=-2+2*index[0],b=1+2*index[1])
   put('in_u'+label+'_create',c1=coord(index,sigma))
   put('in_u'+label+'_annihilate',a1=coord(index,ONE/sigma),a2=coord(TOP,mono(a=-2,b=1,d=2))if high else seq((0,-2),mono(b=-1,c=2)))
   put('in_x'+label+'_create',c1=seq(value=mono(b=-1,d=1))if high else seq((2,0),mono(a=-2,b=1,d=1)),c2=coord(TOP,mono(a=2,b=-1,d=-1))if high else seq((2,0),mono(b=1,d=-1)))
   put('in_x'+label+'_annihilate',a2=seq((-2,0),mono(b=-1,d=1))if high else coord(LOW,mono(b=-1,d=1)))
  else:
   sigma=mono(a=2*index[0],b=1+2*index[1]);qindex=mono(c=2*index[0],d=2*index[1])
   put('left_y'+label+'_create',c2=coord(index))
   put('left_y'+label+'_annihilate',a2=coord(index),c1=seq((2,0),mono(d=1))if high else seq(value=mono(d=1)))
   put('right_u'+label+'_create',c2=coord(index,sigma),a1=seq(value=mono(b=-1,d=1))if high else seq((-2,0),mono(b=-1,d=1)))
   put('right_u'+label+'_annihilate',a2=coord(index,ONE/sigma))
   put('right_v'+label+'_create',c1=seq((2,0),mono(b=1))if high else coord(LOW,mono(b=1)),a2=coord(TOP,mono(d=1)*qindex/sigma)if high else seq((-2,0),mono(b=-1,d=1)))
   put('right_v'+label+'_annihilate',a1=coord(TOP,ONE/sigma)if high else seq((-2,0),mono(b=-1)))
   put('left_x'+label+'_create',c1=coord(index))
   put('left_x'+label+'_annihilate',a1=coord(index),c2=seq((2,0),mono(d=1))if high else seq(value=mono(d=1)))
 return images


def pairing(a,b):return sum(((a[k]*b[k+2]+a[k+2]*b[k]).total()for k in(0,1)),ZERO)


def concentrated(place,specialization=2):
 global Z0
 Z0=specialization;GEOMETRIC_RATIOS.clear();PIVOT_CHECKS.clear();images=initial_images(place)
 before={(name,other):pairing(a,b)for name,a in images.items()for other,b in images.items()}
 if place=='infinity':
  for value in images.values():
   value[1]=value[1]*seq((0,-2));value[3]=value[3]*seq((0,2))
  active=[concentrate_four(images,0,(2,0),'constant'),concentrate_six(images,1,(2,-2))]
 else:active=[concentrate_six(images,0,(2,0)),concentrate_four(images,1,(2,0),'inverse')]
 for (name,other),expected in before.items():assert pairing(images[name],images[other])==expected,(place,name,other)
 for value in images.values():
  for block in(0,1):
   for kind in(block,block+2):
    assert not value[kind].terms,(place,kind,value[kind].terms)
    assert set(value[kind].deltas)<=set(active[block])
 orders=[]
 for block in(0,1):
  for a,b in active[block]:orders.append(((-2+2*block+2*a,1+2*b)if place=='infinity'else(2*a,1+2*b)))
 projected={name:([value[k].at(i)for k in(0,1)for i in active[k]],[value[k+2].at(i)for k in(0,1)for i in active[k]])for name,value in images.items()}
 return {'place':place,'active':active,'orders':orders,'images':projected,'specialization':specialization,
         'geometric_ratios':[{'exponents':list(base),'value':value}for base,value in sorted(GEOMETRIC_RATIOS)],
         'pivot_checks':list(PIVOT_CHECKS)}


if __name__=='__main__':
 for place in('infinity','zero'):
  data=concentrated(place)
  print(place,'operators',len(data['images']),'active',data['active'],'orders',data['orders'])
