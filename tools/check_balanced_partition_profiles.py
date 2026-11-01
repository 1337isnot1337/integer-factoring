#!/usr/bin/env python3
"""Exact profiles for even-uniform and weighted multipartite crossings."""

import argparse
from functools import lru_cache
from itertools import permutations, product
import json
from pathlib import Path
import random
import time

from check_balanced_partition_crossings import F4
from check_permutation_homogeneous_labels import compose, rank_two
from even_partition_crossing import partition_labels, coefficient_terms, coefficient, univariate
from fermionic_piece_polynomial import polynomial_product as pmul

def even_two_block_profile(vertices, block_size):
    start=time.monotonic();a=vertices;k=block_size;q=a//k
    labels=partition_labels(a,k);d=len(labels);index={label:i for i,label in enumerate(labels)};n=d*d

    def transform(mask,action):return sum(1<<action[i]for i in range(a)if mask>>i&1)
    def label_action(action):return tuple(index[tuple(sorted(transform(mask,action)for mask in label))]for label in labels)
    small=tuple(permutations(range(k)));actions=[]
    for p,r,swap in product(small,small,range(2)):
     action=tuple(p[i]+(k if swap else 0)for i in range(k))+tuple(r[i]+(0 if swap else k)for i in range(k))
     actions.append(label_action(action))
    normalizations=[]
    for label in labels:
     physical=[i for mask in label for i in range(a)if mask>>i&1]
     inverse=[0]*a
     for i,j in enumerate(physical):inverse[j]=i
     normalizations.append(label_action(inverse))
     assert normalizations[-1][index[label]]==0
    orbit_of=[-1]*(d**3);representatives=[]
    for aa,bb,uu in product(range(d),repeat=3):
     code=(aa*d+bb)*d+uu
     if orbit_of[code]>=0:continue
     oi=len(representatives);representatives.append((aa,bb,uu))
     for action in actions:orbit_of[(action[aa]*d+action[bb])*d+action[uu]]=oi
    assert min(orbit_of)>=0

    profiles=[(l,t)for l in range(4*q+1)for t in range(5)]
    term_data=[];coeff_data=[];uniform=[]
    for aa,bb,uu in representatives:
     terms=coefficient_terms(labels[aa],labels[bb],labels[uu],labels[0],a,k)
     term_data.append(terms)
     coeff_data.append(tuple(coefficient(terms,l,(t,)*a)for l,t in profiles))
     uniform.append(univariate(terms,(1,)*a))


    for (aa,bb,uu),terms in zip(representatives,term_data):
     expected=int(aa==bb and uu==0)^int(aa==uu and bb==0)^int(aa==0 and bb==uu)
     assert coefficient(terms,2*q,(2,)*a)==expected
    matrices=[[0]*n for _ in profiles];cols=[[0]*n for _ in profiles]
    univdegree=4*(q+a);umats=[[0]*n for _ in range(univdegree+1)];ucols=[[0]*n for _ in range(univdegree+1)]
    for yy,norm in enumerate(normalizations):
     for aa,bb,uu in product(range(d),repeat=3):
      oi=orbit_of[(norm[aa]*d+norm[bb])*d+norm[uu]]
      r=aa*d+yy;c=bb*d+uu;rr=aa*d+uu;cc=bb*d+yy
      rb=1<<c;cb=1<<cc
      for i,value in enumerate(coeff_data[oi]):
       if value:matrices[i][r]|=rb;cols[i][rr]|=cb
      poly=uniform[oi]
      assert poly.bit_count()%2==int(aa==bb and uu==yy)
      while poly:
       bit=poly&-poly;poly^=bit;degree=bit.bit_length()-1
       umats[degree][r]|=rb;ucols[degree][rr]|=cb

    result={'vertices':a,'block_size':k,'labels':d,'quartet_orbits':len(representatives),'uniform_physical_profiles':[], 'uniform_curve_profiles':[]}
    for (l,t),rows,columns in zip(profiles,matrices,cols):
     r,c=rank_two(rows),rank_two(columns)
     if r:result['uniform_physical_profiles'].append([l,t,r,c])
    for i,(rows,columns)in enumerate(zip(umats,ucols)):
     r,c=rank_two(rows),rank_two(columns)
     if r:result['uniform_curve_profiles'].append([i,r,c])
    result['full_crossing_entries_checked']=d**4
    return result

def scale(lo,hi,c):
 if c==0:return 0,0
 if c==1:return lo,hi
 if c==2:return hi,lo^hi
 return lo^hi,lo

def mul(x,y):
 a,b=x;c,d=y;bd=pmul(b,d)
 return pmul(a,c)^bd,pmul(a,d)^pmul(b,c)^bd

def rank_four(rows):
 basis={};r=0
 for lo,hi in rows:
  while lo|hi:
   bit=(lo|hi).bit_length()-1;c=((lo>>bit)&1)+2*((hi>>bit)&1)
   if bit in basis:
    u,v=scale(*basis[bit],c);lo^=u;hi^=v
   else:basis[bit]=scale(lo,hi,(0,1,3,2)[c]);r+=1;break
 return r


def multipartite_profile(part_size):
    start=time.monotonic();b=part_size;k=3;w=(1,2,3);perms=tuple(permutations(range(b)));pairs=tuple(product(perms,repeat=2));d=len(pairs);n=d*d
    index={p:i for i,p in enumerate(pairs)}
    multiplication=[[index[tuple(compose(p[i],q[i])for i in range(2))]for q in pairs]for p in pairs]
    labels=[tuple((i,sigma[i],tau[i])for i in range(b))for sigma,tau in pairs]
    place=[5**i for i in range(b)]

    def kernel(e,f):
     c=int(e==f)
     for j in range(k):
      if e[j]==f[j]:c^=w[j]
     return c

    def polynomial(aa,bb,uu,yy):
     # Every U particle must have an identical creation partner elsewhere.
     if any(e not in aa and e not in bb and e not in yy for e in uu):return(0,0)
     particles=[(g,e)for g,label in (('U',uu),('B',bb),('A',aa),('Y',yy))for e in label]
     size=len(particles);entries=[[(0,0)]*size for _ in range(size)]
     for i,(g,e)in enumerate(particles):
      for j in range(i+1,size):
       h,f=particles[j];pair=g,h;value=0;c=1
       if pair==('U','B')and e==f:value=1^(1<<place[f[0]])
       elif pair==('U','A')and e==f:value=(1<<place[e[0]])^(1<<(place[e[0]]+place[f[0]]))
       elif pair==('U','Y')and e==f:value=1<<place[f[0]]
       elif pair==('B','A'):c=kernel(e,f);value=1<<place[f[0]]
       elif pair==('B','Y'):c=kernel(e,f);value=(1<<place[f[0]])^(1<<(place[f[0]]+place[e[0]]))
       elif pair==('A','Y'):c=kernel(e,f);value=1^(1<<place[e[0]])
       entries[i][j]=scale(value,0,c)
     @lru_cache(None)
     def rec(rem):
      if not rem:return 1,0
      bit=rem&-rem;i=bit.bit_length()-1;rest=rem^bit;todo=rest;lo=hi=0
      while todo:
       bit=todo&-todo;todo^=bit;j=bit.bit_length()-1
       if entries[i][j]!=(0,0):
        x,y=mul(entries[i][j],rec(rest^bit));lo^=x;hi^=y
      return lo,hi
     return rec((1<<size)-1)

    # Fixing Y to identity uses independent relabeling of the two non-row parts.
    records=[]
    for aa,bb,uu in product(range(d),repeat=3):
     lo,hi=polynomial(labels[aa],labels[bb],labels[uu],labels[0])
     assert(lo.bit_count()%2)+2*(hi.bit_count()%2)==int(aa==bb and uu==0)
     if lo|hi:records.append((aa,bb,uu,lo,hi))

    mats={};cols={}
    for aa,bb,uu,lo,hi in records:
     both=lo|hi;codes=[]
     while both:
      bit=both&-both;both^=bit;code=bit.bit_length()-1;c=int(bool(lo&bit))+2*int(bool(hi&bit));codes.append((code,c))
      if code not in mats:mats[code]=[[0,0]for _ in range(n)];cols[code]=[[0,0]for _ in range(n)]
     for yy,action in enumerate(multiplication):
      a1,b1,u1=action[aa],action[bb],action[uu]
      r=a1*d+yy;c=b1*d+u1;rr=a1*d+u1;cc=b1*d+yy
      for code,value in codes:
       if value&1:mats[code][r][0]^=1<<c;cols[code][rr][0]^=1<<cc
       if value&2:mats[code][r][1]^=1<<c;cols[code][rr][1]^=1<<cc

    profile=[]
    for code in sorted(mats):
     r,c=rank_four(mats[code]),rank_four(cols[code]);profile.append([code,r,c])
    maximum=max(r*c for _,r,c in profile)
    res={'labels':d,'pieces':len(profile),'maximum':maximum,'worst':[x for x in profile if x[1]*x[2]==maximum],'profile':profile}
    res['part_size']=part_size
    res['full_crossing_entries_checked']=d**4
    res['normalized_nonzero_quartets']=len(records)
    return res


def reference_rank_four(matrix):
    rows=[row[:] for row in matrix]
    rank=0
    width=len(rows[0]) if rows else 0
    for column in range(width):
        pivot=next((j for j in range(rank,len(rows)) if rows[j][column]),None)
        if pivot is None:
            continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inverse=(0,1,3,2)[rows[rank][column]]
        rows[rank]=[F4[inverse][value] for value in rows[rank]]
        for j in range(rank+1,len(rows)):
            factor=rows[j][column]
            if factor:
                rows[j]=[x^F4[factor][y] for x,y in zip(rows[j],rows[rank])]
        rank+=1
    return rank


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    options=parser.parse_args()
    rng=random.Random(489377)
    for size in range(1,13):
        for _ in range(16):
            matrix=[[rng.randrange(4) for _ in range(size)] for _ in range(size)]
            packed=[(sum((value&1)<<j for j,value in enumerate(row)),
                     sum((value>>1)<<j for j,value in enumerate(row))) for row in matrix]
            assert rank_four(packed)==reference_rank_four(matrix)
    even=[even_two_block_profile(4,2),even_two_block_profile(8,4)]
    weighted=[multipartite_profile(2),multipartite_profile(3)]
    assert [case['maximum'] for case in weighted]==[49,675684]
    assert next(row[2:] for row in even[1]['uniform_physical_profiles']
                if row[:2]==[4,2])==[596,596]
    result={'coefficient_fields':['F2','F4'],
            'independent_F4_rank_checks':192,
            'even_uniform_two_block_cases':even,
            'weighted_three_partite_cases':weighted,
            'scope':'complete specified finite profiles; no asymptotic extrapolation'}
    encoded=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if options.output:
        options.output.write_text(encoded)
    print(json.dumps({'independent_F4_rank_checks':192,
                      'even_uniform_cases':[{k:v for k,v in c.items()
                          if k not in ('uniform_physical_profiles','uniform_curve_profiles')}
                          for c in even],
                      'weighted_cases':[{k:v for k,v in c.items()
                          if k not in ('profile','worst')} for c in weighted]},indent=2))


if __name__=='__main__':
    main()
