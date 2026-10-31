#!/usr/bin/env python3
"""All four-bit downsets: literal kernels, constant ranks and Smith factors."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from downset_crossing import (is_downset,coefficient_matrices,column_blocks,
                              normalized_column_rank,verify_local_normal_form,literal_entry)
from search_auxiliary_coefficient_span import combine,reverse_rank

ROOT = Path(__file__).resolve().parents[1]


def all_downsets(bits=4):
    for mask in range(1,1 << (1 << bits)):
        code = tuple(x for x in range(1 << bits) if mask >> x & 1)
        if is_downset(code):
            yield code


def check():
    cases,identities,entries,ranks,smith = 0,0,0,0,0
    local_checked = set()
    for code in all_downsets():
        cuts,metadata = coefficient_matrices(code)
        blocks = column_blocks(code)
        for _,_,outside in blocks:
            if outside not in local_checked:
                identities += verify_local_normal_form(outside)
                local_checked.add(outside)
        orders = sorted(diff.bit_count()+2*s.bit_count() for diff,_,outside in blocks for s in outside)
        expected = sorted(x.bit_count()+u.bit_count() for x,u in product(code,repeat=2))
        assert orders == expected
        smith += len(orders)
        for i,x in enumerate(code):
            for j,y in enumerate(code):
                for a,u in enumerate(code):
                    for b,v in enumerate(code):
                        polynomial = sum(((matrix[i*len(code)+j] >> (a*len(code)+b))&1) << k
                                         for k,matrix in enumerate(cuts[0]))
                        assert polynomial == literal_entry(x,y,u,v,4)
                        entries += 1
        for mask in range(1 << (metadata['degree']+1)):
            assert reverse_rank(combine(cuts[1],mask)) == normalized_column_rank(blocks,mask)
            ranks += 1
        cases += 1
    assert cases == 167
    paired_cases,paired_entries,paired_ranks,paired_orders = 0,0,0,0
    for first,second in product(tuple(all_downsets(3)),repeat=2):
        cuts,metadata = coefficient_matrices(first,second)
        blocks = column_blocks(first,second)
        orders = sorted(diff.bit_count()+2*s.bit_count() for diff,_,outside in blocks for s in outside)
        assert orders == sorted(x.bit_count()+u.bit_count() for x,u in product(first,second))
        paired_orders += len(orders)
        for i,x in enumerate(first):
            for j,y in enumerate(second):
                for a,u in enumerate(second):
                    for b,v in enumerate(first):
                        polynomial = sum(((matrix[i*len(second)+j] >> (a*len(first)+b))&1) << k
                                         for k,matrix in enumerate(cuts[0]))
                        assert polynomial == literal_entry(x,y,u,v,3)
                        paired_entries += 1
        for mask in range(1 << (metadata['degree']+1)):
            assert reverse_rank(combine(cuts[1],mask)) == normalized_column_rank(blocks,mask)
            paired_ranks += 1
        paired_cases += 1
    assert paired_cases == 361
    assert not is_downset((1,)) and literal_entry(1,1,1,1,1) == 1
    return {'status':'pass','word_bits':4,'complete_nonempty_downsets':cases,
            'literal_entry_polynomials':entries,'constant_basis_rank_comparisons':ranks,
            'local_polynomial_identities':identities,'smith_orders_checked':smith,
            'distinct_downset_fibers_checked':len(local_checked),
            'paired_word_bits':3,'complete_nonempty_downset_pairs':paired_cases,
            'paired_literal_entry_polynomials':paired_entries,
            'paired_constant_basis_rank_comparisons':paired_ranks,
            'paired_smith_orders_checked':paired_orders,
            'non_downset_counterexample':{'code':[1],'actual_determinant':1,'invalid_predicted_order':2},
            'scope':'Exact algebra controls for the all-downset characteristic-two normal-form proof; not a factoring gain or an all-field claim.',
            'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('downset_crossing.py','check_downset_crossing_normal_form.py',
                                           'search_auxiliary_coefficient_span.py')}}


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
