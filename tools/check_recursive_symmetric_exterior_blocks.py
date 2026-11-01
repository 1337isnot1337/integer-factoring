#!/usr/bin/env python3
"""Independent exact controls for recursive symmetric exterior blocks."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random
from cube_symmetric_exterior_blocks import multiplication_matrix,series_smith_orders,symmetric_profile
from recursive_symmetric_exterior_blocks import block_types,coefficient_indices,recursive_blocks
ROOT=Path(__file__).resolve().parents[1]


def profile(values,precision):
    answer=Counter()
    for multiplicity,local,kinds in recursive_blocks(values):
        for order in series_smith_orders(multiplication_matrix(local),precision):answer[order]+=multiplicity
    return answer


def native_profile(values,precision):
    from native_matrix_smith_f4 import MatrixOracle
    answer=Counter();blocks=[]
    for multiplicity,local,kinds in recursive_blocks(values):
        n=len(local)
        oracle=MatrixOracle([[(entry,0) for entry in row] for row in multiplication_matrix(local)],n,precision=precision)
        result=oracle.profile();orders=list(map(int,result['orders']));assert result['orders']==orders
        for order in orders:answer[order]+=multiplicity
        blocks.append({'type':list(kinds),'multiplicity':multiplicity,'dimension':n,'orders':dict(sorted(Counter(orders).items()))})
    return answer,blocks


def main():
    pa=argparse.ArgumentParser();pa.add_argument('--large',action='store_true');pa.add_argument('--output',type=Path,default=Path('/tmp/recursive-symmetric-exterior-blocks.json'));args=pa.parse_args()
    rng=random.Random(1967);dense=one_level=0
    for r in range(1,11):
        for trial in range(4):
            values=[1<<rng.randrange(1,12)]+[rng.randrange(1<<13) for _ in range(r)]
            expected=symmetric_profile(values,32);actual=profile(values,32);assert actual==expected;one_level+=1
            if r<=6:
                original=[values[mask.bit_count()] for mask in range(1<<r)]
                assert Counter(series_smith_orders(multiplication_matrix(original),32))==actual;dense+=1
    samples=[]
    for r in range(1,257):
        types=block_types(r);assert len(types)==(r+1)//2
        assert {sum(kind*(1<<j) for j,kind in enumerate(kinds)) for kinds in types}==set(range(1 if r&1 else 2,r+1,2))
        assert sum(mult*(1<<sum(kinds)) for kinds,mult in types.items())==1<<r
        for kinds,mult in types.items():
            minimal=sum(kind*(1<<j) for j,kind in enumerate(kinds));dimension=1<<sum(kinds)
            assert 4*dimension<=(minimal+2)**2<=(r+2)**2
            if r<=32:
                indices=coefficient_indices(kinds);assert len(indices)==dimension
                assert max(i for i in indices if i is not None)==minimal
        if r in [15,31,63,127,255]:samples.append({'generators':r,'full_dimension':1<<r,'distinct_blocks':len(types),'largest_block':max(1<<sum(k) for k in types)})
    large=[]
    if args.large:
        for axes in [16,32,64]:
            initial_degree=128;exponents=[(initial_degree+j)//(j+1) for j in range(axes)];degree_cap=max((j+1)*e for j,e in enumerate(exponents));values=[1<<e for e in exponents]
            orders,blocks=native_profile(values,512);assert sum(orders.values())==1<<(axes-1)
            assert sum(order*mult for order,mult in orders.items())==(1<<(axes-1))*initial_degree
            assert orders==Counter({2*initial_degree-order:mult for order,mult in orders.items()})
            large.append({'axes':axes,'degree_cap':degree_cap,'coefficient_exponents':exponents,'profile':dict(sorted(orders.items())),'blocks':blocks,'scope':'Exact one-variable specialization profile only; no generic multivariable pressure improvement is asserted.'})
    sources=['tools/recursive_symmetric_exterior_blocks.py','tools/check_recursive_symmetric_exterior_blocks.py','tools/cube_symmetric_exterior_blocks.py']
    if args.large:sources+=['tools/native_matrix_smith_f4.py','tools/native_matrix_smith_f4.cpp','tools/native_rectangular_smith_f4.cpp','tools/native_weighted_smith_f4.py']
    report={'scope':'Simultaneous block conjugacy and exact finite profile controls; no improved factoring exponent.',
       'independent_dense_smith_profiles':dense,'independent_one_level_profiles':one_level,
       'all_type_counts_dimension_sums_and_size_bounds_through':256,'dimension_examples':samples,
       'large_profiles_requested':args.large,'large_specialization_profiles':large,
       'source_sha256':{path:sha256((ROOT/path).read_bytes()).hexdigest() for path in sources}}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if isinstance(v,(str,int,bool))}),flush=True)
if __name__=='__main__':main()
