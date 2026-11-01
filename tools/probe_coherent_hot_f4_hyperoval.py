#!/usr/bin/env python3
"""One explicit F4 hyperoval boundary embedding, with complete paired profiles.

A positive search target, not a claimed improvement: distribute three logical
boundary directions across six physical words using the [6,3,4] F4 hyperoval.
"""
import argparse,json,time
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from coherent_hot_compounds_f4 import height,mul,basis,profiles,validate_maps
from search_coherent_hot_retractions import pressure


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--indices',default='0,1,2,3,4,5');args=ap.parse_args()
    words=sorted(range(1<<12),key=lambda w:(height(w),w));lower=[w for w in words if height(w)<25]
    boundary=[w for w in words if height(w)==25 and w.bit_count()==3];indices=list(map(int,args.indices.split(',')))
    assert len(lower)==92 and len(boundary)==10
    assert len(indices)==6 and len(set(indices))==6 and all(0<=i<10 for i in indices)
    columns=[(1,t,mul(t,t)) for t in range(4)]+[(0,1,0),(0,0,1)]
    for selected in combinations(columns,3):assert len(basis([{i:v for i,v in enumerate(col) if v} for col in selected]))==3
    ps=[{boundary[indices[j]]:columns[j][i] for j in range(6) if columns[j][i]} for i in range(3)]
    ls=[{boundary[indices[j]]:1} for j in (0,4,5)]
    p=[{w:1} for w in lower+[31]]+ps;l=[{w:1} for w in lower+[31]]+ls
    validate_maps(p,l,12)
    names=[Path(__file__).name,'coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','search_coherent_hot_retractions.py']
    r={'completed':False,'indices':indices,'boundary_words':boundary,'hyperoval_columns':columns,'embedding':p,'retraction':l,'independent_triple_checks':20,'source_sha256':{n:sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},'scope':__doc__}
    args.output.write_text(json.dumps(r,indent=2)+'\n');start=time.monotonic()
    result=profiles(p,l,12);fit=pressure(result)
    r.update(completed=True,profile=result,pressure=fit,seconds=time.monotonic()-start)
    args.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'pressure':fit,'seconds':r['seconds']}),flush=True)

if __name__=='__main__':main()
