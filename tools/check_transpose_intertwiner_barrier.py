#!/usr/bin/env python3
"""Constant partial-transpose intertwiners, exact counts and literal jet ranks."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import flatten,local_polynomial,multiply
from check_multipoint_rank_one_seeds import orbit_control
from search_multipoint_rank_one_seeds import KEYS,SWAP,pure_tensors,representatives

ROOT = Path(__file__).resolve().parents[1]


def mm(a,b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) % 3
                 for i,j in product(range(2),repeat=2))


def forms():
    return [p for p in product(range(3),repeat=4)
            if (p[0]*p[3]-p[1]*p[2]) % 3 and next(x for x in p if x) == 1]


def masks(tensor,intertwiners):
    result = []
    for pair in ('yu','xv'):
        mask = 0
        for bit,p in enumerate(intertwiners):
            for a,b in product(range(2),repeat=2):
                f = tuple(tensor[KEYS.index((a,y,u,b) if pair == 'yu' else (y,a,b,u))]
                          for y,u in product(range(2),repeat=2))
                transpose = f[0],f[2],f[1],f[3]
                if mm(transpose,p) != mm(p,f):
                    break
            else:
                mask |= 1 << bit
        result.append(mask)
    return tuple(result)


def column_rank(matrix):
    """Independent column-span reduction, rather than the search's row oracle."""
    basis = {}
    for column in zip(*matrix):
        column = list(column)
        for pivot in sorted(basis):
            scalar = column[pivot]
            if scalar:
                column = [(a-scalar*b) % 3 for a,b in zip(column,basis[pivot])]
        pivot = next((i for i,x in enumerate(column) if x),None)
        if pivot is not None:
            inverse = pow(column[pivot],-1,3)
            basis[pivot] = [x*inverse % 3 for x in column]
    return len(basis)


def independent_jet_ranks(seed,copies,point,cut):
    # Shift first, then form polynomial Kronecker matrices. The search
    # tensors physical entry tables first, then shifts the complete power.
    local = {key:local_polynomial(value,point,3,3) for key,value in seed.items()}
    coefficients = [flatten(local,2,cut,k) for k in range(4)]
    matrix = [[tuple(coefficients[k][i][j] for k in range(4))
               for j in range(4)] for i in range(4)]
    work = [[[1]]]
    for _ in range(copies):
        n = len(work)
        work = [[multiply(work[i//4][j//4],matrix[i%4][j%4],3)
                 for j in range(4*n)] for i in range(4*n)]
    return [column_rank([[value[k] for value in row] for row in work])
            for k in range(3*copies+1)]


def check():
    intertwiners = forms()
    assert len(intertwiners) == 24
    tensors,reps = pure_tensors(),representatives()
    tensor_masks = [masks(t,intertwiners) for t in tensors]
    histogram = Counter(tensor_masks)
    row_counts = []
    for a in reps:
        am = masks(a,intertwiners)
        direct = sum(bool(am[0]&bm[0]&cm[0] or am[1]&bm[1]&cm[1])
                     for bm in tensor_masks for cm in tensor_masks)
        grouped = sum(bcount*ccount*bool(am[0]&bm[0]&cm[0] or am[1]&bm[1]&cm[1])
                      for bm,bcount in histogram.items() for cm,ccount in histogram.items())
        assert direct == grouped
        row_counts.append(direct)
    assert row_counts == [180225,54833,54833,71409,71409,86625]
    orbit_sizes = orbit_control()
    weighted = sum(a*b for a,b in zip(row_counts,orbit_sizes))
    assert weighted == 32454657
    saved = json.loads((ROOT/'evidence/multipoint-actual-jets.json').read_text())
    assert len(saved['records']) == 200 and saved['copies'] == [1,2,3]
    selected = {0,15,16,17,47,70,72,127,128,165,197,199}
    certified = []
    independent = 0
    minima = {1:10000,2:10000,3:10000}
    for record in saved['records']:
        ai,bi,ci = record['indices']
        a,b,c = reps[ai],tensors[bi],tensors[ci]
        am,bm,cm = masks(a,intertwiners),tensor_masks[bi],tensor_masks[ci]
        common = (am[0]&bm[0]&cm[0],am[1]&bm[1]&cm[1])
        seed = {key:(x,(y-e-z)%3,(-x-y-e)%3,z)
                for key,x,y,z,e in zip(KEYS,a,b,c,SWAP)}
        if any(common):
            place = 0 if common[0] else 1
            bit = (common[place]&-common[place]).bit_length()-1
            certified.append({'signature_index':record['signature_index'],
                              'pair':('yu','xv')[place],'intertwiner':intertwiners[bit]})
        for result in record['results']:
            copies = result['copies']
            ranks = result['complete_jet_cut_ranks']
            assert len(ranks) == 3 and all(len(x) == 3*copies+1 for x in ranks)
            if any(common):
                assert all(a == b for row in ranks for a,b in row)
            if record['signature_index'] in selected:
                for j,point in enumerate((0,2,None)):
                    for cut in (0,1):
                        actual = independent_jet_ranks(seed,copies,point,cut)
                        assert actual == [pair[cut] for pair in ranks[j]]
                        independent += len(actual)
            minima[copies] = min(minima[copies],result['best_rank_product'])
    assert len(certified) == 66
    assert minima == {1:4,2:36,3:729}
    sources = ('check_transpose_intertwiner_barrier.py','search_multipoint_actual_jets.py')
    return {'status':'pass','projective_invertible_intertwiners':24,
            'excluded_cubic_representatives_by_first_value':row_counts,
            'excluded_cubic_representatives':sum(row_counts),
            'total_cubic_representatives':6*513**2,
            'excluded_actual_coefficient_seeds':weighted,
            'total_actual_coefficient_seeds':513**3,
            'actual_jet_witnesses':200,'copies_checked':[1,2,3],
            'certified_rank_equal_witnesses':certified,
            'independent_kronecker_column_rank_values':independent,
            'finite_best_rank_products':minima,
            'scope':'A uniform rank-equivalence obstruction for constant partial-transpose intertwiners; the 200 actual-jet witnesses do not exhaust coefficient seeds.',
            'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in sources},
            'actual_jet_evidence_sha256':sha256((ROOT/'evidence/multipoint-actual-jets.json').read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'certified_rank_equal_witnesses'},indent=2))


if __name__ == '__main__':
    main()
