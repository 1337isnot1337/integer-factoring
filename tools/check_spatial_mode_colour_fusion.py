#!/usr/bin/env python3
"""Exact spatial-mode colour fusion and its stopped fixed-q8 pressure.

Matrix/rank and literal-control scaffolding adapted with credit from
Sol worker 2's check_fixed_curve_quartic_floor (PR111,705b3e3).
No complete profile, alphabet search or stronger factoring bound.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb
from pathlib import Path
import check_fixed_curve_quartic_floor as s
import check_colour_alternation as exact

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/spatial-mode-colour-fusion.json'

@lru_cache(None)
def alphas(word):
 colors=[c for _,c in word]
 if len(set(colors))!=4 or s.charge(colors):return (0,0)
 out=[]
 for axis in (0,1):
  aa=sorted((c,v) for v,c in word if not (v>>axis&1))
  bb=sorted((c,v) for v,c in word if v>>axis&1)
  if len(aa)!=2:out.append(0);continue
  exp=(6,10,0) if axis==0 else (10,6,0)
  out.append(int(all(exp in exact.KERNEL[a][b] for (_,a),(_,b) in zip(aa,bb))))
 return tuple(out)

@lru_cache(None)
def coefficient(word,axis):
 ax,ay=alphas(word)
 # Exact degree32 coefficient before specialization/truncation.
 terms=set()
 for enabled,exp in ((ax,(12,20)),(ax^ay,(16,16)),(ay,(20,12))):
  if enabled:terms.add(exp)
 assert all(x+y==32 and x>=12 and y>=12 for x,y in terms)
 assert len(terms)%2==0 # endpoint111 remains unchanged
 value=0
 for exp in terms:value^=1<<exp[axis]
 return value&s.MASK

def entry(aa,bb,left,right,pair,axis):
 total=len(aa)+len(bb)
 if not total:return 1
 if not aa or not bb or total>6:return 0
 out=0
 if len(aa)==len(bb):
  for perm in permutations(bb):
   if any(c!=d for (_,c),(_,d) in zip(aa,perm)):continue
   p=1
   for (i,_),(j,_) in zip(aa,perm):p=s.pmul(p,pair[i][j])
   out^=p
 if total==4:
  word=tuple(sorted([(left[i],c) for i,c in aa]+[(right[j],c) for j,c in bb]))
  out^=coefficient(word,axis)
 if total==6:
  for i,c in aa:
   for j,d in bb:
    if c!=d:continue
    word=tuple(sorted([(left[k],v) for k,v in aa if k!=i]+[(right[k],v) for k,v in bb if k!=j]))
    out^=s.pmul(coefficient(word,axis),pair[i][j])
 return out

def matrix(axis,target):
 left,right=s.cube.axis_parts(axis)
 def evaluate(poly):
  out=0
  for a,b,c in poly:out^=1<<(2*(a if axis==0 else b)+c)
  return out&s.MASK
 pair=[[evaluate(s.kernel[i][j]) for j in right] for i in left]
 words=[w for w in product(range(8),repeat=4) if s.charge(w)==target]
 active=[tuple((i,c) for i,c in enumerate(w) if c) for w in words]
 blocks=[[0]*512 for _ in range(19)]
 for i,aa in enumerate(active):
  for j,bb in enumerate(active):
   p=entry(aa,bb,left,right,pair,axis)
   while p:
    bit=p&-p;blocks[bit.bit_length()-1][i]|=1<<j;p^=bit
 return blocks

def literal_controls(axis):
 """Actual support expansion versus the fast cross-matching formula."""
 left,right=s.cube.axis_parts(axis)
 def evaluate(poly):
  value=0
  for a,b,c in poly:value^=1<<(2*(a if axis==0 else b)+c)
  return value&s.MASK
 pair=[[evaluate(s.kernel[i][j]) for j in right] for i in left]
 edge={(a,b):pair[i][j] for i,a in enumerate(left) for j,b in enumerate(right)}
 edge.update({(b,a):v for (a,b),v in list(edge.items())})
 @lru_cache(None)
 def pfaffian(active):
  if not active:return 1
  (a,color),rest=active[0],active[1:];out=0
  for i,(b,other) in enumerate(rest):
   if color==other and (a,b) in edge:
    out^=s.pmul(edge[a,b],pfaffian(rest[:i]+rest[i+1:]))
  return out
 def literal(active):
  out=pfaffian(active)
  for support in combinations(range(len(active)),4):
   quartet=tuple(active[i] for i in support)
   if not (any(v in left for v,_ in quartet) and any(v in right for v,_ in quartet)):
    continue
   rest=tuple(v for i,v in enumerate(active) if i not in support)
   out^=s.pmul(coefficient(quartet,axis),pfaffian(rest))
  return out
 panel=set(product(range(2),repeat=8))
 for quartet in combinations(range(1,8),4):
  if s.charge(quartet):continue
  for support in combinations(range(8),4):
   if not (set(support)&set(left) and set(support)&set(right)):continue
   for colors in permutations(quartet):
    word=[0]*8
    for v,c in zip(support,colors):word[v]=c
    panel.add(tuple(word))
    for a in set(left)-set(support):
     for b in set(right)-set(support):
      for c in range(1,8):
       extended=word.copy();extended[a]=extended[b]=c;panel.add(tuple(extended))
 for word in panel:
  aa=tuple((i,word[v]) for i,v in enumerate(left) if word[v])
  bb=tuple((i,word[v]) for i,v in enumerate(right) if word[v])
  active=tuple((v,c) for v,c in enumerate(word) if c)
  assert entry(aa,bb,left,right,pair,axis)==literal(active),(axis,word)
 return len(panel)


def pressure_certificate(combined):
 assert combined==Counter({0:1,6:14,12:35,15:84,16:70,17:14,18:208})
 dimension=4096
 tail=dimension-sum(combined.values())
 known_sum=sum(k*v for k,v in combined.items())
 mean=Fraction(56*dimension-known_sum,tail)
 assert tail==3670 and known_sum==6866 and mean==Fraction(22251,367)<61
 # t=u^3, subtract (1+7*t^(32/3))^3 and factor u^18.
 coefficients={3*k-18:v for k,v in combined.items() if k}
 coefficients[165]=tail
 for k,v in {14:-21,46:-147,78:-343}.items():
  coefficients[k]=coefficients.get(k,0)+v
 degree=300
 bernstein=[sum(Fraction(a*comb(k,j),comb(degree,j))
                for j,a in coefficients.items() if j<=k) for k in range(degree+1)]
 assert min(bernstein)>0
 for j in range(degree+1):
  assert comb(degree,j)*sum((-1)**(j-k)*comb(j,k)*bernstein[k]
                          for k in range(j+1))==coefficients.get(j,0)
 encoded=[[b.numerator,b.denominator] for b in bernstein]
 return dict(dimension=dimension,determinant_degree=56*dimension,
             known_orders=sum(combined.values()),known_order_sum=known_sum,
             tail_count=tail,tail_mean=str(mean),tail_cap=61,
             polynomial=coefficients,bernstein_degree=degree,
             positive_coefficients=len(bernstein),minimum=str(min(bernstein)),
             coefficient_sha256=sha256(json.dumps(encoded).encode()).hexdigest(),
             paired_type_pressure_floor='3/4',target='1/sqrt(2)')


def compute():
 record=dict(scope='Uniform spatial-mode quartet construction, ordinary all-q leading prefix and fixed-q8 STOP certificate; no full profile or factoring improvement.',
             precision=19,cases=[],axes=[])
 for axis in (0,1):
  combined=Counter()
  for target in range(8):
   blocks=matrix(axis,target)
   result=s.ranks(blocks)
   predicted={0:1,12:7} if target==0 else {6:2,12:4}
   assert {k:v for k,v in result['smith_prefix'].items() if k<=12}==predicted
   check=s.column_rank(blocks,19)
   assert check==result['toeplitz_ranks'][19]
   combined.update(result['smith_prefix'])
   record['cases'].append(dict(axis=axis,charge=target,dimension=512,
                               **result,independent_column_rank_19=check))
  record['axes'].append(dict(axis=axis,whole_prefix=dict(combined),
                             literal_entry_controls=literal_controls(axis),
                             pressure=pressure_certificate(combined)))
 names=['tools/check_spatial_mode_colour_fusion.py',
        'tools/check_fixed_curve_quartic_floor.py','tools/check_colour_alternation.py',
        'tools/check_quartic_charge_fusion.py','tools/cube_homogeneous_gaussian.py']
 record['source_sha256']={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
 return json.loads(json.dumps(record))


if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path)
 args=parser.parse_args()
 result=compute()
 if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
 else:assert result==json.loads(RECEIPT.read_text()),'receipt differs; inspect before regeneration'
 print('Spatial-mode colour fusion: whole leading layer35, direct entries and q8 STOP floor PASS')
