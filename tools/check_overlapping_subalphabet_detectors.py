#!/usr/bin/env python3
"""Exact F2 functional certificates for all two/three three-label overlaps."""
from __future__ import annotations

import argparse
from itertools import combinations, permutations, product
import json
from pathlib import Path
import time


# The order of each functional's coordinates is (E,R,C,H_1,...,H_m).
CASES = [
    ("adjacent complement edges", ((0,1,2),(0,1,3)),
     [(3,(1,0,0,0,0)), (2,(1,0,0,0,1)), (2,(1,0,0,1,0)),
      (4,(1,0,1,0,0)), (4,(1,1,0,0,0))]),
    ("disjoint complement edges", ((0,1,2),(2,3,4)),
     [(2,(1,0,0,0,0)), (3,(1,0,0,0,1)), (3,(1,0,0,1,0)),
      (4,(1,0,1,0,0)), (4,(1,1,0,0,0))]),
    ("complement triangle", ((0,1,2),(0,1,3),(0,1,4)),
     [(2,(1,0,0,0,0,0)), (1,(1,0,0,0,0,1)), (2,(1,0,0,0,1,0)),
      (2,(1,0,0,1,0,0)), (4,(1,0,1,0,0,0)), (3,(1,1,0,0,0,0)),
      (1,(1,1,0,0,0,1))]),
    ("complement star", ((0,1,2),(0,1,3),(0,2,3)),
     [(3,(1,0,0,0,0,0)), (2,(1,0,0,0,1,0)), (4,(1,0,1,0,0,0)),
      (3,(1,1,0,0,0,0)), (1,(1,1,0,0,0,1))]),
    ("complement four-vertex path", ((0,1,2),(0,1,3),(0,2,4)),
     [(2,(1,0,0,0,0,0)), (1,(1,0,0,0,0,1)), (3,(1,0,0,0,1,0)),
      (4,(1,0,1,0,0,0)), (3,(1,1,0,0,0,0)), (1,(1,1,0,0,0,1))]),
    ("complement three-vertex path plus edge", ((0,1,2),(2,3,4),(0,1,3)),
     [(1,(1,0,0,0,0,0)), (1,(1,0,0,0,0,1)), (3,(1,0,0,0,1,0)),
      (2,(1,0,0,1,0,0)), (4,(1,0,1,0,0,0)), (2,(1,1,0,0,0,0)),
      (1,(1,1,0,0,0,1))]),
]


def gf2_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def canonical_family(family, label_permutations):
    return min(tuple(sorted(tuple(sorted(perm[label] for label in subset))
                            for subset in family)) for perm in label_permutations)


def make_bases(subsets):
    d = 5
    dimension = 3+len(subsets)
    full = [0]*dimension
    row = [[0]*(d*d) for _ in range(dimension)]
    column = [[0]*(d*d) for _ in range(dimension)]
    for index, (x,y,u,v) in enumerate(product(range(d), repeat=4)):
        e = int(x == v and y == u)
        r = int(x == y and u == v)
        c = int(x == u and y == v)
        values = [e,r,c]+[(e^r^c) if {x,y,u,v} <= set(subset) else 0
                          for subset in subsets]
        for j, value in enumerate(values):
            if value:
                full[j] |= 1 << index
                row[j][d*x+y] |= 1 << (d*u+v)
                column[j][d*x+u] |= 1 << (d*y+v)
    assert gf2_rank(full) == dimension
    return row, column


def linear_combination(bases, coefficients):
    rows = [0]*len(bases[0])
    for basis, coefficient in zip(bases, coefficients):
        if coefficient:
            rows = [a^b for a,b in zip(rows,basis)]
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    label_permutations = list(permutations(range(5)))
    all_subsets = list(combinations(range(5), 3))
    expected = {canonical_family(subsets, label_permutations): name
                for name, subsets, cert in CASES}
    assert len(expected) == 6
    orbit_sizes = {key:0 for key in expected}
    family_counts = {}
    for size in (2,3):
        count = 0
        observed = set()
        for family in combinations(all_subsets,size):
            key = canonical_family(family, label_permutations)
            assert key in expected
            orbit_sizes[key] += 1
            observed.add(key)
            count += 1
        assert len(observed) == (2 if size == 2 else 4)
        family_counts[size] = count
    assert family_counts == {2:45, 3:120}
    results = []
    for name, subsets, certificate in CASES:
        row_bases, column_bases = make_bases(subsets)
        dimension = len(row_bases)
        profiles = []
        for coefficients in product((0,1), repeat=dimension):
            if not any(coefficients):
                continue
            rr = gf2_rank(linear_combination(row_bases, coefficients))
            cr = gf2_rank(linear_combination(column_bases, coefficients))
            detector_rank = sum(weight for weight, functional in certificate
                                if sum(a*b for a,b in zip(coefficients,functional)) % 2)
            assert detector_rank**2 <= rr*cr
            profiles.append({'coefficients':coefficients, 'cut_ranks':[rr,cr],
                             'detector_rank':detector_rank})
        crossing_bound = sum(weight for weight, functional in certificate if functional[0])
        assert crossing_bound >= 13
        assert crossing_bound**2 > 5**3
        results.append({'orbit':name, 'subalphabets':subsets,
                        'labelled_families':orbit_sizes[canonical_family(subsets,label_permutations)],
                        'basis_rank':dimension, 'certificate':certificate,
                        'verified_lower_base':crossing_bound,
                        'nonzero_tensor_cases':len(profiles), 'profiles':profiles})
    assert sum(len(result['profiles']) for result in results) == 314
    result = {'field_size':2, 'alphabet':5, 'subalphabet_size':3,
              'family_counts':family_counts, 'orbit_count':len(results),
              'nonzero_tensor_cases':314, 'uniform_spectral_lower_base':13,
              'scope':'Only F2 observable individual-word rank spectral transfers in the specified entry spaces; arbitrary finite auxiliary size.',
              'orbits':results, 'seconds':time.monotonic()-started}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    summary = {**result, 'orbits':[{key:value for key,value in item.items() if key != 'profiles'}
                                  for item in results]}
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
