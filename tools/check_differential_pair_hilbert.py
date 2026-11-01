#!/usr/bin/env python3
"""Exact weighted-Hilbert and binary/F4 dimension comparisons."""
import json
from pathlib import Path
from hashlib import sha256
from itertools import product
from collections import Counter
from differential_pair_hilbert import hilbert
from coherent_differential_kernels_f4 import differential_kernel
from check_coherent_differential_kernels_f4 import dense_rank


def operators(pairs,C):
    return [{(1<<a)|(1<<b):row[i] for i,(a,b) in enumerate(pairs) if row[i]} for row in C]


def main():
    pairs=[(3,4),(2,5),(1,6),(0,7)];families={'mds_two':[[1,1,1,1],[0,1,2,3]],'weighted_one':[[1,1,2,3]]};checks=[]
    for name,C in families.items():
        for modes,cap in [(8,64),(12,16),(12,20),(12,22),(12,24)]:
            prediction,ranks=hilbert(pairs,C,modes);prediction=Counter({k:v for k,v in prediction.items() if k[0]<=cap})
            actual=differential_kernel(operators(pairs,C),modes,cap);observed=Counter(zip(actual['heights'],actual['charges']))
            assert prediction==observed
            checks.append({'family':name,'modes':modes,'cap':cap,'dimension':actual['q'],'graded_cells':len(observed)})
    dims=[]
    for columns in product([(1,0),(0,1),(1,1)],repeat=4):
        C=[list(row) for row in zip(*columns)]
        if dense_rank(C)!=2:continue
        P,_=hilbert(pairs,C,8);dims.append(sum(P.values()))
    assert len(dims)==78 and min(dims)==96
    md,_=hilbert(pairs,families['mds_two'],8);assert sum(md.values())==92
    C=[[1,0,1,1],[0,1,1,0]];binary=differential_kernel(operators(pairs,C),8,64);assert binary['q']==96
    names=[Path(__file__).name,'differential_pair_hilbert.py','coherent_differential_kernels_f4.py','coherent_hot_compounds_f4.py']
    report={'status':'passed','weighted_cell_comparisons':checks,'binary_nonzero_column_rank_two_matrices':len(dims),'minimum_binary_dimension':96,'F4_MDS_dimension':92,'independent_binary_example_dimension':binary['q'],'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Exact physical-code dimension and weighted Hilbert counts; a smaller dimension alone is not a better crossing pressure.'}
    Path('evidence/differential-pair-hilbert-controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
