#!/usr/bin/env python3
"""Independent original-B determinant controls for the ten-point F8 code."""
import sys,json,random,ctypes
from functools import lru_cache
from pathlib import Path
from hashlib import sha256
from coherent_f8_maps import tensor_vector,mul,height
from coherent_hot_compounds import one_particle
from coherent_hot_compounds_f8 import library
from check_coherent_hot_f4_hyperoval_entries import pmul
@lru_cache(None)
def determinant(a,b):
 rows=[i for i in range(24)if a>>i&1];cols=[i for i in range(24)if b>>i&1];B=one_particle(12);v={0:1}
 for row in rows:
  nxt={}
  for mask,c in v.items():
   for j,col in enumerate(cols):
    if not mask>>j&1 and B[row][col]:nxt[mask|1<<j]=nxt.get(mask|1<<j,0)^pmul(c,B[row][col])
  v=nxt
 return v.get((1<<len(cols))-1,0)
def pack(vs):
 off=[0];masks=[];scalars=[]
 for v in vs:
  for w,c in sorted(v.items()):masks.append(w);scalars.append(c)
  off.append(len(masks))
 return [(ctypes.c_uint32*len(x))(*x)for x in [off,masks,scalars]]
d=json.load(open('evidence/f8-hyperoval-q96-full-profile.json'));p=[{int(k):v for k,v in x.items()}for x in d['embedding']];l=[{int(k):v for k,v in x.items()}for x in d['retraction']];charges=[next(iter(v)).bit_count() for v in p];rng=random.Random(80961);entries=nonbinary=top_entries=0;blocks=[]
for ch in sorted({a+b for a in charges for b in charges}):
 available=[(x,y)for x in range(96)for y in range(96)if charges[x]+charges[y]==ch];special=[v for v in available if max(v)>=93];chosen=[rng.choice(special or available)];remaining=[v for v in available if v!=chosen[0]]
 if remaining:chosen.append(rng.choice(remaining))
 for cut in ('zero','infinity'):
  rows=[tensor_vector(l[x],l[y]if cut=='zero' else p[y],12)for x,y in chosen];cols=[tensor_vector(p[y]if cut=='zero' else l[y],p[x],12)for x,y in chosen];n=len(rows);words=2;orders=(ctypes.c_int64*n)();out=(ctypes.c_uint64*(3*n*n*words))()
  status=library(words,127).coherent_hot_compound_f8(12,n,*pack(rows),*pack(cols),orders,out);assert status==0,status
  for i,row in enumerate(rows):
   for j,col in enumerate(cols):
    planes=[0,0,0]
    for a,ca in row.items():
     for b,cb in col.items():
      term=determinant(a,b);c=mul(ca,cb)
      for k in range(3):
       if c>>k&1:planes[k]^=term
    native=[sum(int(out[(3*(i*n+j)+k)*words+w])<<(64*w)for w in range(words))for k in range(3)]
    assert planes==native,(ch,cut,i,j)
    assert sum((v.bit_count()%2)<<k for k,v in enumerate(planes))==int(i==j)
    entries+=1;nonbinary+=bool(planes[1]or planes[2]);top_entries+=bool(planes[2])
  blocks.append({'charge':ch,'cut':cut,'pairs':chosen})
 print(json.dumps({'charge':ch,'entries':entries,'nonbinary':nonbinary,'top_plane_entries':top_entries}),flush=True)
names=['coherent_f8_maps.py','coherent_hot_compounds_f8.py','coherent_hot_compounds_f8.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','coherent_hot_compounds.py']
r={'completed':True,'complete_entries':entries,'nonbinary_entries':nonbinary,'top_plane_entries':top_entries,'independent_binary_determinants':determinant.cache_info().currsize,'blocks':blocks,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names},'scope':'Selected full polynomial entries of the q96 F8 hyperoval; not an independent full profile replay.'}
Path('evidence/f8-hyperoval-q96-entry-controls.json').write_text(json.dumps(r,indent=2)+'\n')
