#!/usr/bin/env python3
"""Replay the finite catalogue, with independent Toeplitz and group controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import check_power, local_profiles
from search_multipoint_rank_one_seeds import (
    KEYS, SWAP, pure_tensors, representatives, run, seed_from_values,
)

ROOT = Path(__file__).resolve().parents[1]


def orbit_control():
    group = [(a,b,c,d) for a,b,c,d in product(range(3), repeat=4)
             if (a*d-b*c) % 3]
    assert len(group) == 48
    pair_orbits = []
    for partner in (0, 1):
        orbit = set()
        for a,b,c,d in group:
            inverse = pow((a*d-b*c) % 3, -1, 3)
            left = (a,c)
            # Inverse transpose acts on the paired factor.
            right = ((d,-b) if partner == 0 else (-c,a))
            right = tuple(x*inverse % 3 for x in right)
            orbit.add(tuple(left[i]*right[j] % 3 for i,j in product(range(2),repeat=2)))
        pair_orbits.append(orbit)
    covered = {(0,)*16}
    sizes = [1]
    for u,v in product(range(2), repeat=2):
        for scalar in ((1,2) if (u,v) == (0,0) else (1,)):
            orbit = {tuple(scalar*p[2*x+v1]*q[2*y+u1] % 3
                           for x,y,u1,v1 in KEYS)
                     for p in pair_orbits[v] for q in pair_orbits[u]}
            assert not covered.intersection(orbit)
            covered.update(orbit)
            sizes.append(len(orbit))
    assert covered == set(pure_tensors())
    assert len(sizes) == len(representatives())
    return sizes


def check(saved):
    catalogue = run()
    assert saved == json.loads(json.dumps(catalogue))
    assert catalogue['profiled_seeds'] == 3078
    assert catalogue['prefix_certificate_failures'] == 0
    counts = 0
    for record in catalogue['catalogue']:
        seed = seed_from_values(*record['values'])
        assert tuple(sum(seed[key]) % 3 for key in KEYS) == SWAP
        for point, pair, certificate in zip((0,2,None), record['profiles'],
                                            record['prefix_certificates']):
            assert tuple(map(tuple,local_profiles(seed,2,2,point,3))) == pair
            threshold,k,l = certificate
            threshold = Fraction(threshold)
            assert k*l >= 8
            assert Fraction(sum(pair[0][:k]), k) <= threshold
            assert Fraction(sum(pair[1][:l]), l) <= threshold
        assert sum(Fraction(x[0]) for x in record['prefix_certificates']) <= 2
        check_power(seed,2,2,2,[0,2,None],[2,2,1],3)
        counts += 1
    return {'status': 'pass', 'orbit_sizes': orbit_control(),
            'seeds_replayed': catalogue['profiled_seeds'],
            'independent_toeplitz_profiles': 6*counts,
            'reconstructed_two_copy_crossings': counts,
            'complete_two_copy_entries': counts*2**8,
            'exact_rho_floor': '3/4',
            'scope': catalogue['scope'],
            'source_sha256': {name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('search_multipoint_rank_one_seeds.py',
                                           'check_multipoint_rank_one_seeds.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path,
                        default=ROOT/'evidence/multipoint-rank-one-seeds.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check(json.loads(args.catalogue.read_text()))
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
