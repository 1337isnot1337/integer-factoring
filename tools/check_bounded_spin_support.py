#!/usr/bin/env python3
"""Exact bounded-word support, channel unranking and both extreme cut factors."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from bounded_spin_crossing import paths, recursive_vector, path_norm
from bounded_spin_support import (supported_word, orbit_count, orbit_word, vertex_factor,
                                  unrank_subset,corner_channel_count,corner_words)
from check_bounded_spin_six_vertex import block_tensor
from check_bounded_spin_rank_algebra import rank

ROOT = Path(__file__).resolve().parents[1]


def check():
    from itertools import combinations
    for size in range(8):
        for weight in range(size+1):
            expected = list(combinations(range(size),weight))
            assert [tuple(unrank_subset(size,weight,j)) for j in range(len(expected))] == expected
    support_vectors,complement_entries,word_controls = 0,0,0
    dimensions = []
    for n in range(2,13,2):
        ps = paths(n)
        support = {w for w in range(1 << n) if supported_word(n,w)}
        representatives = {w for w in support if not w & 1}
        assert len(support) == 2*orbit_count(n)
        assert {orbit_word(n,j) for j in range(orbit_count(n))} == representatives
        word_controls += 1 << n
        for p in ps:
            v = recursive_vector(p)
            assert set(v) <= support
            for w,c in v.items():
                assert v.get(w^((1 << n)-1),0) == (-1)**(n//2)*c
                complement_entries += 1
            support_vectors += 1
        dimensions.append(orbit_count(n))
    assert dimensions == [1,3,10,34,116,396]
    for n in range(6,41,2):
        assert orbit_count(n) == 4*orbit_count(n-2)-2*orbit_count(n-4)
    cases,corner_controls = [],0
    for n in (2,4,6,8):
        ps = paths(n)
        monomials,t = block_tensor(n)
        slots = {m:i for i,m in enumerate(monomials)}
        extreme = t[slots[n//2,0,0]]
        assert np.array_equal(extreme,t[slots[0,n//2,0]])
        assert np.array_equal(extreme,t[slots[0,0,n//2]])
        norms = np.array([path_norm(p) for p in ps],dtype=np.int64)
        normalized = extreme.copy()
        for r,s in product(range(len(ps)),repeat=2):
            normalized[:,:,r,s] = normalized[:,:,r,s]*pow(int(norms[r]*norms[s]),-1,5) % 5
        channels = orbit_count(n)
        row_left = np.array([[vertex_factor(p,q,j,scale=2) for j in range(channels)]
                             for p,q in product(ps,repeat=2)],dtype=np.int64)
        row_right = np.array([[vertex_factor(r,s,j,True,True) for r,s in product(ps,repeat=2)]
                              for j in range(channels)],dtype=np.int64)
        column_left = np.array([[vertex_factor(p,r,j,dual_second=True,scale=2)
                                 for j in range(channels)] for p,r in product(ps,repeat=2)],dtype=np.int64)
        column_right = np.array([[vertex_factor(q,s,j,dual_second=True)
                                  for q,s in product(ps,repeat=2)] for j in range(channels)],dtype=np.int64)
        row_matrix = normalized.reshape(len(ps)**2,-1)
        column_matrix = normalized.transpose(0,2,1,3).reshape(len(ps)**2,-1)
        assert np.array_equal(row_left @ row_right % 5,row_matrix)
        assert np.array_equal(column_left @ column_right % 5,column_matrix)
        actual = [rank(row_matrix),rank(column_matrix)]
        assert max(actual) <= channels
        assert vertex_factor(ps[0],ps[0],channels) == 0
        cases.append({'length':n,'labels':len(ps),'channels':channels,'rank_pair':actual})
        if n >= 4:
            vectors = [recursive_vector(p) for p in ps]
            for b,c in ((1,0),(0,1)):
                channel_words = [corner_words(n,b,c,j) for j in range(corner_channel_count(n,b,c))]
                left,right = [],[]
                for x,y,u,v in channel_words:
                    left.append([2*p.get(x,0)*q.get(y,0) % 5 for p,q in product(vectors,repeat=2)])
                    right.append([p.get(u,0)*q.get(v,0)*pow(int(r*s),-1,5) % 5
                                  for (p,q),(r,s) in zip(product(vectors,repeat=2),product(norms,repeat=2))])
                actual = np.array(left,dtype=np.int64).T @ np.array(right,dtype=np.int64) % 5
                expected = t[slots[n//2-b-c,b,c]].copy()
                for r,s in product(range(len(ps)),repeat=2):
                    expected[:,:,r,s] = expected[:,:,r,s]*pow(int(norms[r]*norms[s]),-1,5) % 5
                assert np.array_equal(actual,expected.reshape(len(ps)**2,-1))
                corner_controls += 1
    return {'status':'pass','prime':5,'support_path_vectors':support_vectors,
            'nonzero_complement_entry_controls':complement_entries,
            'exhaustive_physical_words':word_controls,'orbit_dimensions':dimensions,
            'rank_factorization_cases':cases,'recurrence_checked_through_length':40,
            'complete_first_corner_factorizations':corner_controls,
            'scope':'Uniform explicit support and corner rank bounds; corner bands alone do not form a full crossing decomposition.',
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ('tools/bounded_spin_crossing.py','tools/bounded_spin_support.py',
                              'tools/check_bounded_spin_support.py','tools/check_bounded_spin_six_vertex.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
