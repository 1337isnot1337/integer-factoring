#!/usr/bin/env python3
"""Independent full-polynomial entries for a selected q96 F8 retraction."""
import argparse,ctypes,json,random
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from coherent_f8_maps import tensor_vector,mul,validate_maps
from coherent_hot_compounds import one_particle
from coherent_hot_compounds_f8 import library
from probe_f8_hyperoval_retraction import maps

def pmul(a,b):
    out=0
    while b:
        bit=b&-b;out^=a<<(bit.bit_length()-1);b^=bit
    return out

@lru_cache(None)
def determinant(a,b):
    rows=[i for i in range(24)if a>>i&1];cols=[i for i in range(24)if b>>i&1]
    assert len(rows)==len(cols);B=one_particle(12);values={0:1}
    for row in rows:
        following={}
        for mask,value in values.items():
            for j,col in enumerate(cols):
                if not mask>>j&1 and B[row][col]:following[mask|1<<j]=following.get(mask|1<<j,0)^pmul(value,B[row][col])
        values=following
    return values.get((1<<len(cols))-1,0)

def pack(vectors):
    offsets=[0];masks=[];scalars=[]
    for vector in vectors:
        for word,c in sorted(vector.items()):masks.append(word);scalars.append(c)
        offsets.append(len(masks))
    return [(ctypes.c_uint32*len(v))(*v)for v in [offsets,masks,scalars]]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('profile',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    data=json.loads(args.profile.read_text());assert data['completed']
    for n,h in data['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
    p,l,_,_=maps(data['pivot_columns']);assert json.loads(json.dumps(p))==data['embedding']and json.loads(json.dumps(l))==data['retraction']
    hs,ch=validate_maps(p,l,12);assert hs==data['profile']['heights']and ch==data['profile']['charges']
    rng=random.Random(80963);entries=nonbinary=third=0;blocks=[]
    for charge in sorted({a+b for a in ch for b in ch}):
        pairs=[(x,y)for x in range(96)for y in range(96)if ch[x]+ch[y]==charge];special=[v for v in pairs if max(v)>=93]
        chosen=[rng.choice(special or pairs)];remaining=[v for v in pairs if v not in chosen]
        chosen+=rng.sample(remaining,min(2,len(remaining)))
        for cut in ('zero','infinity'):
            rows=[tensor_vector(l[x],l[y]if cut=='zero'else p[y],12)for x,y in chosen]
            cols=[tensor_vector(p[y]if cut=='zero'else l[y],p[x],12)for x,y in chosen]
            n=len(rows);words=2;orders=(ctypes.c_int64*n)();coeff=(ctypes.c_uint64*(3*n*n*words))()
            status=library(words,127).coherent_hot_compound_f8(12,n,*pack(rows),*pack(cols),orders,coeff);assert status==0
            for i,row in enumerate(rows):
                for j,col in enumerate(cols):
                    planes=[0,0,0]
                    for a,ca in row.items():
                        for b,cb in col.items():
                            term=determinant(a,b);c=mul(ca,cb)
                            for k in range(3):
                                if c>>k&1:planes[k]^=term
                    native=[sum(int(coeff[(3*(i*n+j)+k)*words+w])<<(64*w)for w in range(words))for k in range(3)]
                    assert planes==native
                    assert sum((v.bit_count()%2)<<k for k,v in enumerate(planes))==int(i==j)
                    entries+=1;nonbinary+=bool(planes[1]or planes[2]);third+=bool(planes[2])
            blocks.append({'charge':charge,'cut':cut,'pairs':chosen})
    assert entries==166 and third>0
    names=[Path(__file__).name,'probe_f8_hyperoval_retraction.py','coherent_f8_maps.py','coherent_hot_compounds_f8.py','coherent_hot_compounds_f8.cpp','coherent_hot_compounds.py','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp']
    out={'completed':True,'pivot_columns':data['pivot_columns'],'full_polynomial_entries':entries,'nonbinary_entries':nonbinary,'third_plane_entries':third,'independent_binary_determinants':determinant.cache_info().currsize,'blocks':blocks,'profile_sha256':sha256(args.profile.read_bytes()).hexdigest(),'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names},'scope':'Independent selected complete entries from original-B subset determinants, not a second full matrix/profile reconstruction.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ('blocks','source_sha256','scope')}))
if __name__=='__main__':main()
