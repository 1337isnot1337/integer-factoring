#!/usr/bin/env python3
"""Independent group, full-Toeplitz and crossing controls for rank-two seeds."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import permutations,product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import check_power,local_profiles,rank
from search_multipoint_rank_one_seeds import KEYS,SWAP,seed_from_values
from search_multipoint_rank_two_seeds import low_rank_representatives,low_rank_tensors,run

ROOT = Path(__file__).resolve().parents[1]


def inverse(g):
    a,b,c,d = g
    scalar = pow((a*d-b*c) % 3,-1,3)
    return tuple(x*scalar % 3 for x in (d,-b,-c,a))


def mm(a,b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) % 3
                 for i,j in product(range(2),repeat=2))


def orbit_coverage():
    group = [g for g in product(range(3),repeat=4) if (g[0]*g[3]-g[1]*g[2]) % 3]
    symmetries = [p for p in permutations(range(4))
                  if all(SWAP[KEYS.index(tuple(key[i] for i in p))] == e
                         for key,e in zip(KEYS,SWAP))]
    assert len(symmetries) == 8
    maps = [[KEYS.index(tuple(key[i] for i in p)) for key in KEYS] for p in symmetries]
    covered = set()
    for representative in low_rank_representatives():
        # Every representative is e00(x,y) beta(u,v).
        beta = tuple(representative[KEYS.index((0,0,u,v))]
                     for u,v in product(range(2),repeat=2))
        orbit = set()
        for g,h in product(group,repeat=2):
            hi = inverse(h)
            hit = (hi[0],hi[2],hi[1],hi[3])
            transformed = mm(mm(hit,beta),inverse(g))
            tensor = tuple(g[2*x]*h[2*y]*transformed[2*u+v] % 3
                           for x,y,u,v in KEYS)
            orbit.add(tensor)
        for tensor in orbit:
            covered.update(tuple(tensor[i] for i in mapping) for mapping in maps)
    expected = set(low_rank_tensors())
    assert covered == expected
    # Direct dense ranks check every finite value, including the cut-swap additions.
    for tensor in expected:
        first = [[tensor[KEYS.index((x,y,u,v))] for u,v in product(range(2),repeat=2)]
                 for x,y in product(range(2),repeat=2)]
        second = [[tensor[KEYS.index((x,y,u,v))] for y,v in product(range(2),repeat=2)]
                  for x,u in product(range(2),repeat=2)]
        assert rank(first,3)*rank(second,3) <= 2
    return len(covered)


def check(saved):
    catalogue = run()
    assert saved == json.loads(json.dumps(catalogue))
    assert catalogue['profiled_seeds'] == 32265
    assert catalogue['prefix_certificate_failures'] == 0
    count = 0
    for record in catalogue['catalogue']:
        seed = seed_from_values(*record['values'])
        for point,pair,certificate in zip((0,2,None),record['profiles'],record['prefix_certificates']):
            assert tuple(map(tuple,local_profiles(seed,2,2,point,3))) == pair
            bound,k,l = certificate
            assert k*l >= 8
            assert max(Fraction(sum(pair[0][:k]),k),Fraction(sum(pair[1][:l]),l)) <= Fraction(bound)
        assert sum(Fraction(x[0]) for x in record['prefix_certificates']) <= 2
        check_power(seed,2,2,2,[0,2,None],[2,2,1],3)
        count += 1
    return {'status': 'pass', 'finite_values_covered_by_stabilizer': orbit_coverage(),
            'seeds_replayed': catalogue['profiled_seeds'],
            'independent_toeplitz_profiles': 6*count,
            'reconstructed_two_copy_crossings': count,
            'complete_two_copy_entries': count*2**8,
            'exact_rho_floor': '3/4', 'scope': catalogue['scope'],
            'source_sha256': {name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('search_multipoint_rank_one_seeds.py',
                                           'search_multipoint_rank_two_seeds.py',
                                           'check_multipoint_rank_two_seeds.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue',type=Path,default=ROOT/'evidence/multipoint-rank-two-seeds.json')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check(json.loads(args.catalogue.read_text()))
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
