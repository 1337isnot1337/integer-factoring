#!/usr/bin/env python3
"""Compare the subset Gram criterion with full differential-kernel pairings."""
from itertools import product
from pathlib import Path
from hashlib import sha256
import json
from coherent_hot_compounds_f4 import mul,dot,profiles
from coherent_differential_kernels_f4 import differential_kernel
from orthogonal_differential_pairs_f4 import orthogonal_retraction,subset_gram_certificate


def main():
    records=[];accepted=[];endpoint_entries=0
    for scalars in product((1,2,3),repeat=4):
        C=[list(scalars),[mul(s,t) for s,t in zip(scalars,(0,1,2,3))]]
        ops=[{(1<<a)|(1<<(7-a)):row[a] for a in range(4) if row[a]} for row in C]
        code=differential_kernel(ops,8,64);assert code['q']==92
        cert=subset_gram_certificate(C);L=orthogonal_retraction(code)
        assert (L is not None)==cert['all_subsets_nondegenerate']
        if L is not None:
            accepted.append(scalars)
            assert all(dot(l,p)==int(i==j) for i,l in enumerate(L) for j,p in enumerate(code['embedding']))
            endpoint_entries+=92**2
            # The projection's ambient matrix is symmetric, not necessarily
            # a representation whose individual basis vectors have norm one.
            P=code['embedding']
            for a in range(256):
                for b in range(a,256):
                    ab=ba=0
                    for p,l in zip(P,L):ab^=mul(p.get(a,0),l.get(b,0));ba^=mul(p.get(b,0),l.get(a,0))
                    assert ab==ba
        records.append({'column_scalars':scalars,'criterion':cert,'full_kernel_nondegenerate':L is not None})
    assert len(accepted)==18
    small_profiles=[]
    for c in (2,3):
        code=differential_kernel([{9:1,6:c}],4,16);L=orthogonal_retraction(code);assert L is not None
        result=profiles(code['embedding'],L,4)
        assert result['zero_orders_t']==result['infinity_orders_t']
        small_profiles.append({'coefficient':c,'q':10,'orders':result['zero_orders_t']})
    names=[Path(__file__).name,'orthogonal_differential_pairs_f4.py','coherent_differential_kernels_f4.py','coherent_hot_compounds_f4.py']
    report={'status':'passed','two_cut_controls':small_profiles,'coefficient_matrices':81,'subset_gram_checks':81*16,'full_kernel_dimension':92,'accepted':accepted,'endpoint_pairing_entries':endpoint_entries,'symmetric_projection_entries':18*256*257//2,'records':records,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Exact construction and Gram criterion controls; no crossing-pressure or factoring improvement asserted.'}
    Path('evidence/orthogonal-differential-pairs-f4.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256','accepted','two_cut_controls')}))

if __name__=='__main__':main()
