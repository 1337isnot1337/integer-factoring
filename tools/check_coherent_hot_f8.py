#!/usr/bin/env python3
"""Independent full-entry and partial-profile F8 controls."""
import ctypes,json
from hashlib import sha256
from pathlib import Path
from functools import lru_cache
from coherent_f8_maps import height,mul,validate_maps,tensor_vector
from coherent_hot_compounds import one_particle
from coherent_hot_compounds_f8 import partial_block,library
from binary_extension_smith import ExtensionMatrixOracle,BinaryExtension
from check_coherent_hot_f4_hyperoval_entries import pmul
@lru_cache(None)
def determinant(a,b):
 rows=[i for i in range(8) if a>>i&1];cols=[i for i in range(8) if b>>i&1];B=one_particle(4);v={0:1}
 for i in rows:
  nxt={}
  for mask,c in v.items():
   for j,k in enumerate(cols):
    if not mask>>j&1 and B[i][k]:nxt[mask|1<<j]=nxt.get(mask|1<<j,0)^pmul(c,B[i][k])
  v=nxt
 return v.get((1<<len(cols))-1,0)
def pack(vs):
 off=[0];masks=[];scalars=[]
 for v in vs:
  for w,c in sorted(v.items()):masks.append(w);scalars.append(c)
  off.append(len(masks))
 return [(ctypes.c_uint32*len(x))(*x)for x in [off,masks,scalars]]
entries=profiles=0
spectra={}
field=BinaryExtension(11)
assert all(mul(a,b)==field.multiply(a,b) for a in range(8) for b in range(8))
for scalar in range(1,8):
 lower=[w for w in range(16) if height(w)<8];p=[{w:1} for w in lower]+[{6:1,9:scalar}];l=[{w:1}for w in lower]+[{6:1}];_,charges=validate_maps(p,l,4);q=len(p)
 for cut in ('zero','infinity'):
  for ch in sorted({a+b for a in charges for b in charges}):
   pairs=[(x,y)for x in range(q)for y in range(q) if charges[x]+charges[y]==ch]
   rows=[tensor_vector(l[x],l[y]if cut=='zero' else p[y],4)for x,y in pairs];cols=[tensor_vector(p[y]if cut=='zero' else l[y],p[x],4)for x,y in pairs];n=len(rows);words=2;cutoff=64
   orders=(ctypes.c_int64*n)();out=(ctypes.c_uint64*(3*n*n*words))();status=library(words,cutoff).coherent_hot_compound_f8(4,n,*pack(rows),*pack(cols),orders,out);assert status==0,status
   matrix=[]
   for i,row in enumerate(rows):
    line=[]
    for j,col in enumerate(cols):
     planes=[0,0,0]
     for a,ca in row.items():
      for b,cb in col.items():
       d=determinant(a,b);c=mul(ca,cb)
       for k in range(3):
        if c>>k&1:planes[k]^=d
     native=[sum(int(out[(3*(i*n+j)+k)*words+w])<<(64*w)for w in range(words))for k in range(3)]
     assert planes==native;(entries:=entries+1)
     line.append([sum(((p>>t)&1)<<k for k,p in enumerate(planes))for t in range(max(p.bit_length()for p in planes))])
    matrix.append(line)
   expected=ExtensionMatrixOracle(matrix,11,n).profile()['orders'];assert sorted(orders)==expected
   spectra[scalar,cut,ch]=expected
   for K in [8,16,24]:
    r=partial_block(4,rows,cols,K);assert r['retained_orders_t']==[x for x in expected if x<=K];assert r['tail_count']==sum(x>K for x in expected)
   profiles+=1
 print(json.dumps({'scalar':scalar,'entries':entries,'profiles':profiles}),flush=True)
for (scalar,cut,ch),orders in spectra.items():assert spectra[mul(scalar,scalar),cut,ch]==orders
invalid=0
for p,l in [([{0:1}],[{0:2}]),([{0:1,1:1}],[{0:1}]),([{1:1}],[{1:1}]),([{0:8}],[{0:1}]),([{0:1.0}],[{0:1}])]:
 try:validate_maps(p,l,4)
 except ValueError:invalid+=1
 else:raise AssertionError('Invalid map accepted')
for args in [(4,[{0:1}],[{0:1}],-1),(0,[{0:1}],[{0:1}],8),(4,[{0:8}],[{0:1}],8),(4,[],[],8)]:
 try:partial_block(*args)
 except ValueError:invalid+=1
 else:raise AssertionError('Invalid block accepted')
assert invalid==9
names=['check_coherent_hot_f8.py','coherent_f8_maps.py','coherent_hot_compounds_f8.py','coherent_hot_compounds_f8.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','binary_extension_smith.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','coherent_hot_compounds.py','check_coherent_hot_f4_hyperoval_entries.py']
out={'completed':True,'full_entries':entries,'profile_blocks':profiles,'independent_binary_determinants':determinant.cache_info().currsize,'partial_cutoffs':[8,16,24],'field_products':64,'frobenius_profile_comparisons':len(spectra),'invalid_inputs_rejected':invalid,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names}}
Path('evidence/coherent-hot-f8-controls.json').write_text(json.dumps(out,indent=2)+'\n')

