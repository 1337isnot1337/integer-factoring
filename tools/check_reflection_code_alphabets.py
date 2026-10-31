#!/usr/bin/env python3
"""Independent exact controls for physical code-alphabet restrictions."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path

from check_block_smith_crossings import reverse_polynomial,smith_profile,toeplitz_rank
from search_reflection_code_alphabets import Search


def check_case(search,alphabet):
    q=search.q;degree=search.seed['degree'];pairs=list(product(alphabet,repeat=2));records=[]
    for cut in(0,1):
        matrix=[]
        for a,b in pairs:
            row=[]
            for c,d in pairs:
                x,y,u,v=(a,b,c,d)if cut==0 else(a,c,b,d)
                value=search.seed['values'][x+q*y+q*q*u+q*q*q*v][0]
                assert value.bit_count()%2==int(x==v and y==u)
                if cut==0:
                    value<<=(search.energy[a]+search.energy[b])//2
                else:
                    value=reverse_polynomial(value,degree)<<((search.energy[b]+search.energy[d])//2)
                row.append(value)
            matrix.append(row)
        bound=max(value.bit_length()-1 for row in matrix for value in row)
        orders,precision=smith_profile(matrix,bound)
        restored=[value-(degree if cut else 0)for value in orders]
        assert restored==search.profile(alphabet,cut)
        checks=0
        for length in sorted({1}|{value+1 for value in orders}):
            assert toeplitz_rank(matrix,length)==sum(max(0,length-value)for value in orders)
            checks+=1
        records.append({'cut':cut,'orders':dict(Counter(restored)),'independent_precision':precision,'toeplitz_checks':checks})
    return {'interval':search.M,'alphabet':list(alphabet),'profiles':records,'endpoint_entries':len(alphabet)**4}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    results=[];small=Search(4)
    for size in range(2,5):
        for alphabet in combinations(range(4),size):results.append(check_case(small,alphabet))
    larger=Search(8)
    for alphabet in((0,4,8),(0,4,8,12),(0,3,12,15)):
        results.append(check_case(larger,alphabet))
    sources=[Path(__file__),Path(__file__).with_name('search_reflection_code_alphabets.py'),
             Path(__file__).with_name('reflection_dual_rail_gaussian.py'),Path(__file__).with_name('check_block_smith_crossings.py')]
    report={'status':'pass','cases':results,'independent_profiles':2*len(results),
            'toeplitz_checks':sum(p['toeplitz_checks']for r in results for p in r['profiles']),
            'endpoint_entries':sum(r['endpoint_entries']for r in results),
            'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest()for p in sources},
            'scope':'exact coordinate-restriction, endpoint, native/Python Smith and Toeplitz controls; no improved pressure certificate'}
    rendered=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':main()
