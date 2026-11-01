#!/usr/bin/env python3
"""Independent exact selected-entry controls for the F4 q96 hyperoval code.

Subset dynamic programming computes determinants of the original B(t), not
of the native height-factored K(t). This is not another full Smith replay.
"""
import json,random
from pathlib import Path
from functools import lru_cache
from hashlib import sha256
from coherent_hot_compounds import one_particle
from coherent_hot_compounds_f4 import tensor_vector,compound_block,mul


def pmul(a,b):
    out=0
    while b:
        bit=b&-b;out^=a<<(bit.bit_length()-1);b^=bit
    return out


@lru_cache(None)
def determinant(rowmask,colmask):
    rows=[i for i in range(24) if rowmask>>i&1];cols=[i for i in range(24) if colmask>>i&1]
    assert len(rows)==len(cols);n=len(rows);B=one_particle(12)
    values={0:1}
    for row in rows:
        next_values={}
        for mask,value in values.items():
            for j,col in enumerate(cols):
                if not mask>>j&1 and B[row][col]:
                    target=mask|(1<<j);next_values[target]=next_values.get(target,0)^pmul(value,B[row][col])
        values=next_values
    return values.get((1<<n)-1,0)


def main():
    path=Path('evidence/coherent-hot-f4-hyperoval-q96.json');data=json.loads(path.read_text());assert data['completed']
    p=[{int(w):c for w,c in v.items()} for v in data['embedding']];l=[{int(w):c for w,c in v.items()} for v in data['retraction']]
    charges=data['profile']['charges'];rng=random.Random(73296);entries=nonbinary=0;blocks=[]
    for charge in sorted({a+b for a in charges for b in charges}):
        available=[(x,y) for x in range(96) for y in range(96) if charges[x]+charges[y]==charge]
        special=[v for v in available if max(v)>=93]
        chosen=[rng.choice(special or available)]
        remaining=[v for v in available if v!=chosen[0]]
        if remaining:chosen.append(rng.choice(remaining))
        for cut in ('zero','infinity'):
            rows=[tensor_vector(l[x],l[y] if cut=='zero' else p[y],12) for x,y in chosen]
            cols=[tensor_vector(p[y] if cut=='zero' else l[y],p[x],12) for x,y in chosen]
            native=compound_block(12,rows,cols,coefficients=True)
            for i,row in enumerate(rows):
                for j,col in enumerate(cols):
                    result=[0,0]
                    for a,ca in row.items():
                        for b,cb in col.items():
                            f=mul(ca,cb);term=determinant(a,b)
                            if f&1:result[0]^=term
                            if f&2:result[1]^=term
                    assert tuple(result)==native['matrix'][i][j],(charge,cut,i,j)
                    assert (result[0].bit_count()%2 | (result[1].bit_count()%2)<<1)==int(i==j)
                    assert max(x.bit_length()-1 for x in result)<=native['degree_bound_t']
                    entries+=1;nonbinary+=bool(result[1])
            blocks.append({'charge':charge,'cut':cut,'pairs':chosen,'degree_bound_t':native['degree_bound_t']})
    assert nonbinary>0
    names=[Path(__file__).name,'coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','coherent_hot_compounds.py']
    report={'status':'passed','scope':__doc__,'profile_sha256':sha256(path.read_bytes()).hexdigest(),'complete_F4_entries':entries,'nonbinary_entries':nonbinary,'independent_determinants':determinant.cache_info().currsize,'blocks':blocks,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names}}
    Path('evidence/coherent-hot-f4-hyperoval-entry-controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('blocks','source_sha256','scope')}))

if __name__=='__main__':main()
