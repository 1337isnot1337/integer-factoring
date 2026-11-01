#!/usr/bin/env python3
"""Binary-limit and independent dense F4 controls for differential kernels."""
import json
from pathlib import Path
from hashlib import sha256
from coherent_differential_kernels import differential_kernel as binary
from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds import support,height

TABLE=((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))


def dense_rank(rows):
    if not rows:return 0
    a=[r[:] for r in rows];rank=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank];v=a[rank][j];inv=TABLE[v][v];a[rank]=[TABLE[inv][x] for x in a[rank]]
        for i in range(rank+1,len(a)):
            c=a[i][j]
            if c:a[i]=[x^TABLE[c][y] for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank


def main():
    cases=[];binary_controls=0;dense_controls=0;annihilation=0
    for cap in (8,12,16):
        for coefficient in (1,2,3):
            ops=[{6:1,9:coefficient}];code=differential_kernel(ops,4,cap)
            words=[w for w in range(16) if height(w)<=cap];images=[]
            for oi,op in enumerate(ops):
                for target in range(16):
                    images.append([next((c for removed,c in op.items() if w&removed==removed and w^removed==target),0) for w in words])
            assert len(words)-dense_rank(images)==code['q'];dense_controls+=1
            for vec in code['embedding']:
                v=[vec.get(w,0) for w in words]
                for row in images:
                    total=0
                    for x,y in zip(row,v):total^=TABLE[x][y]
                    assert total==0;annihilation+=1
            if coefficient==1:
                old=binary([sum(1<<w for w in ops[0])],4,cap)
                assert code['embedding']==[{w:1 for w in support(v)} for v in old['embedding']]
                assert code['pivot_words']==old['pivot_words'];binary_controls+=1
            cases.append({'cap':cap,'coefficient':coefficient,'q':code['q']})
    names=[Path(__file__).name,'coherent_differential_kernels_f4.py','coherent_hot_compounds_f4.py','coherent_differential_kernels.py']
    r={'status':'passed','cases':cases,'binary_basis_comparisons':binary_controls,'independent_dense_nullities':dense_controls,'literal_F4_annihilation_equations':annihilation,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names}}
    Path('evidence/coherent-differential-kernels-f4-controls.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))

if __name__=='__main__':main()
