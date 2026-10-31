#!/usr/bin/env python3
"""Actual F2 coefficient-span optima on no-adjacent-ones physical codes."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
from random import Random

from downset_crossing import (coefficient_matrices,column_blocks,is_downset,
                              literal_entry,normalized_column_rank)
from search_auxiliary_coefficient_span import combine,reverse_rank,minimum_threshold,rank_sum_bounds

ROOT = Path(__file__).resolve().parents[1]


def codewords(length,weight=None,prefix=None):
    code = [x for x in range(1 << length) if not x&(x >> 1)
            and (weight is None or x.bit_count() == weight)]
    return code if prefix is None else code[:prefix]


def literal_controls(code,cuts,length):
    if len(code) <= 8:
        quadruples = product(range(len(code)),repeat=4)
    else:
        random = Random(314159)
        quadruples = [tuple(random.randrange(len(code)) for _ in range(4)) for _ in range(128)]
    checked = 0
    for i,j,a,b in quadruples:
        actual = sum(((matrix[i*len(code)+j] >> (a*len(code)+b))&1) << k
                     for k,matrix in enumerate(cuts[0]))
        assert actual == literal_entry(code[i],code[j],code[a],code[b],length)
        checked += 1
    return checked


def build_record(length,weight,prefix):
    code = codewords(length,weight,prefix)
    cuts,metadata = coefficient_matrices(code)
    controls = literal_controls(code,cuts,length)
    pairs = [[reverse_rank(combine(cut,mask)) for cut in cuts]
             for mask in range(1 << (metadata['degree']+1))]
    blocks = column_blocks(code) if is_downset(code) else None
    if blocks is not None:
        assert all(pair[1] == normalized_column_rank(blocks,mask) for mask,pair in enumerate(pairs))
    threshold = minimum_threshold(pairs)
    sums = rank_sum_bounds(pairs,20)
    lower,upper,scale = (sums['optimal_rank_sum_lower_numerator'],
                        sums['optimal_rank_sum_upper_numerator'],sums['denominator'])
    return {'length':length,'weight':weight,'prefix':prefix,'codewords':code,**metadata,
            'literal_kernel_controls':controls,'constant_normal_form_controls':len(pairs) if blocks else 0,
            'rank_pairs_by_functional_mask':pairs,**threshold,**sums,
            'rank_sum_lower_exceeds_rho_three_quarters':lower**2 > scale**2*len(code)**3,
            'rank_sum_upper_below_rho_seven_tenths':upper**5 < scale**5*len(code)**7}


def verify_record(record):
    code = codewords(record['length'],record['weight'],record['prefix'])
    assert record['codewords'] == code
    cuts,metadata = coefficient_matrices(code)
    assert all(record[k] == value for k,value in metadata.items())
    assert literal_controls(code,cuts,record['length']) == record['literal_kernel_controls']
    pairs = record['rank_pairs_by_functional_mask']
    assert len(pairs) == 1 << (metadata['degree']+1)
    blocks = column_blocks(code) if is_downset(code) else None
    for mask,pair in enumerate(pairs):
        assert pair == [reverse_rank(combine(cut,mask)) for cut in cuts]
        if blocks is not None:
            assert pair[1] == normalized_column_rank(blocks,mask)
    target,threshold = len(pairs)-1,record['minimum_maximum_rank_product']
    dual = record['strict_lower_threshold_separating_covector']
    assert (target & dual).bit_count() % 2 == 1
    assert all((w & dual).bit_count() % 2 == 0 for w,(a,b) in enumerate(pairs) if a*b < threshold)
    value = 0
    for mask in record['decomposition_masks']:
        value ^= mask
        assert pairs[mask][0]*pairs[mask][1] <= threshold
    assert value == target
    scale = record['denominator']
    squares = [a*b*scale**2 for a,b in pairs]
    floors = [isqrt(q) for q in squares]
    potential = record['rank_sum_lower_potential']
    assert len(potential) == len(pairs) and potential[0] == 0
    lower = record['optimal_rank_sum_lower_numerator']
    assert potential[target] == lower
    assert all(potential[v^w] <= potential[v]+weight
               for v in range(len(pairs)) for w,weight in enumerate(floors))
    value,upper = 0,0
    for mask in record['upper_decomposition_masks']:
        value ^= mask
        upper += floors[mask]+int(floors[mask]**2 != squares[mask])
    assert value == target and upper == record['optimal_rank_sum_upper_numerator']
    assert record['rank_sum_lower_exceeds_rho_three_quarters'] == (lower**2 > scale**2*len(code)**3)
    assert record['rank_sum_upper_below_rho_seven_tenths'] == (upper**5 < scale**5*len(code)**7)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    args = parser.parse_args()
    assert args.output or args.verify
    if args.verify:
        saved = json.loads(args.verify.read_text())
        assert saved['status'] == 'pass' and saved['prime'] == 2
        for name,digest in saved['source_sha256'].items():
            assert digest == sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
        for record in saved['records']:
            verify_record(record)
        print(json.dumps({'status':'pass','cases':len(saved['records']),
                          'complete_functionals':sum(len(r['rank_pairs_by_functional_mask']) for r in saved['records'])}))
        return
    cases = [(m,None,None) for m in range(1,8)]
    cases += [(8,2,None),(8,3,None),(10,2,None),(10,3,None),
              (3,None,4),(6,None,16),(7,None,32)]
    records = []
    for length,weight,prefix in cases:
        record = build_record(length,weight,prefix)
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k not in
                          ('rank_pairs_by_functional_mask','rank_sum_lower_potential','codewords')}),flush=True)
    result = {'status':'pass','prime':2,'records':records,
              'scope':'Complete finite F2 recombinations of the elementary crossing on the stated nonlinear codes. No uniform upper bound or stronger factoring exponent.',
              'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                                for name in ('search_regular_code_crossings.py','downset_crossing.py',
                                             'search_auxiliary_coefficient_span.py')}}
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n')


if __name__ == '__main__':
    main()
