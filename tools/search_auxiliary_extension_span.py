#!/usr/bin/env python3
"""Full F4 recombinations of small binary auxiliary coefficient spaces."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
from random import Random

import check_cancellation_module_automata as module
from check_cancellation_module_crossing import local_crossing
from check_cubic_full_image_span import rank_dense
from search_auxiliary_coefficient_span import coefficient_matrices,minimum_threshold,rank_sum_bounds

ROOT = Path(__file__).resolve().parents[1]


def combine(matrices,mask):
    # Field basis (omega,omega^2): digit 3 represents 1, so the endpoint
    # mask is still all bits one. Addition is XOR of the packed digits.
    first,second = [0]*len(matrices[0]),[0]*len(matrices[0])
    k = 0
    while mask:
        value = mask & 3
        if value & 1:
            first = [a^b for a,b in zip(first,matrices[k])]
        if value & 2:
            second = [a^b for a,b in zip(second,matrices[k])]
        mask >>= 2
        k += 1
    return list(zip(first,second))


def extension_rank(rows,reverse=False):
    basis = {}
    for a,b in reversed(rows) if reverse else rows:
        while a|b:
            pivot = (((a|b)&-(a|b)) if reverse else (a|b)).bit_length()-1
            value = ((a >> pivot)&1)+2*((b >> pivot)&1)
            if pivot not in basis:
                if value == 1:  # multiply by omega^2
                    a,b = a^b,a
                elif value == 2:  # multiply by omega
                    a,b = b,a^b
                basis[pivot] = a,b
                break
            c,d = basis[pivot]
            if value == 1:
                a,b = a^d,b^c^d
            elif value == 2:
                a,b = a^c^d,b^c
            else:
                a,b = a^c,b^d
    return len(basis)


def dense_rank(rows,columns):
    # Independently implemented field arithmetic uses the basis (1,omega).
    matrix = [[((b >> j)&1)+2*(((a^b) >> j)&1) for j in range(columns)]
              for a,b in rows]
    return rank_dense(matrix,4,reverse=True)


def field_controls():
    random = Random(271828)
    for n,m in product(range(1,9),repeat=2):
        rows = [(random.randrange(1 << m),random.randrange(1 << m)) for _ in range(n)]
        assert extension_rank(rows) == dense_rank(rows,m)
    return 64


def analyze(name,local,degree,length,precision):
    matrices,endpoints = coefficient_matrices(local,degree,length)
    count = 1 << (2*len(matrices[0]))
    pairs = [[extension_rank(combine(cut,mask)) for cut in matrices] for mask in range(count)]
    threshold = minimum_threshold(pairs)
    sums = rank_sum_bounds(pairs,precision)
    controls = set(threshold['decomposition_masks']+sums['upper_decomposition_masks']+[count-1])
    for mask in controls:
        assert [dense_rank(combine(cut,mask),len(cut[0])) for cut in matrices] == pairs[mask]
    for binary_mask in range(1 << len(matrices[0])):
        encoded = sum(3 << (2*k) for k in range(len(matrices[0])) if binary_mask >> k & 1)
        for cut in range(2):
            rows = combine(matrices[cut],encoded)
            assert all(a == b for a,b in rows)
            assert pairs[encoded][cut] == len(module.binary_basis([a for a,b in rows]))
    return {'family':name,'copies':length,'degree_bound':degree*length,
            'coefficient_functionals_enumerated':count,'endpoint_entries':endpoints,
            'rank_pairs_by_functional_mask':pairs,**threshold,**sums,
            'independent_dense_rank_controls':len(controls)}


def verify_record(record,local,degree):
    matrices,endpoints = coefficient_matrices(local,degree,record['copies'])
    pairs = record['rank_pairs_by_functional_mask']
    assert len(pairs) == 1 << (2*len(matrices[0]))
    assert endpoints == record['endpoint_entries']
    for mask,pair in enumerate(pairs):
        assert [extension_rank(combine(cut,mask),reverse=True) for cut in matrices] == pair
    target,threshold = len(pairs)-1,record['minimum_maximum_rank_product']
    dual = record['strict_lower_threshold_separating_covector']
    assert (target & dual).bit_count() % 2 == 1
    assert all((w & dual).bit_count() % 2 == 0 for w,(a,b) in enumerate(pairs) if a*b < threshold)
    value = 0
    for mask in record['decomposition_masks']:
        assert pairs[mask][0]*pairs[mask][1] <= threshold
        value ^= mask
    assert value == target
    scale = record['denominator']
    squares = [a*b*scale**2 for a,b in pairs]
    floors = [isqrt(q) for q in squares]
    potential = record['rank_sum_lower_potential']
    assert len(potential) == len(pairs) and potential[0] == 0
    assert potential[target] == record['optimal_rank_sum_lower_numerator']
    assert all(potential[v^w] <= potential[v]+weight
               for v in range(len(pairs)) for w,weight in enumerate(floors))
    value,upper = 0,0
    for mask in record['upper_decomposition_masks']:
        value ^= mask
        upper += floors[mask]+int(floors[mask]**2 != squares[mask])
    assert value == target and upper == record['optimal_rank_sum_upper_numerator']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    args = parser.parse_args()
    assert args.output or args.verify
    constructors = {'elementary':lambda:module.elementary(1),'module':lambda:(local_crossing(),2),
                    'pairing':module.pairing,'mixed-loop':module.mixed_loop}
    field_checks = field_controls()
    if args.verify:
        saved = json.loads(args.verify.read_text())
        assert saved['field'] == 4 and saved['status'] == 'pass'
        for name,digest in saved['source_sha256'].items():
            assert digest == sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
        for record in saved['records']:
            verify_record(record,*constructors[record['family']]())
        print(json.dumps({'status':'pass','records':len(saved['records']),
                          'field_arithmetic_controls':field_checks,
                          'complete_functionals':sum(r['coefficient_functionals_enumerated'] for r in saved['records'])}))
        return
    records = []
    for name,lengths in (('elementary',(1,2,3,4)),('module',(1,2)),('pairing',(1,2)),('mixed-loop',(1,))):
        local,degree = constructors[name]()
        for length in lengths:
            record = analyze(name,local,degree,length,20)
            records.append(record)
            print(json.dumps({k:v for k,v in record.items() if k not in
                              ('rank_pairs_by_functional_mask','rank_sum_lower_potential')}),flush=True)
    result = {'status':'pass','field':4,'field_polynomial':'omega^2+omega+1',
              'coefficient_digit_basis':['omega','omega^2'],'records':records,
              'field_arithmetic_controls':field_checks,
              'scope':'Complete F4 coefficient recombinations only at the stated finite lengths; no algebraic-closure or asymptotic optimum.',
              'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                                for name in ('search_auxiliary_extension_span.py','search_auxiliary_coefficient_span.py',
                                             'check_cancellation_module_automata.py','check_cancellation_module_crossing.py',
                                             'check_cubic_full_image_span.py')}}
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n')


if __name__ == '__main__':
    main()
