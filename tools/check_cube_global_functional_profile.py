#!/usr/bin/env python3
"""Finite support controls for the all-B global-orientation separator."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_cube_favorable_reweighting import report as single_report
from check_two_axis_global_pairing import check_assignment_dependent_pairing

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/cube-global-functional-profile.json'


def tensor_product(first,second):
    """Two binary truth vectors, with first-block assignment the low byte."""
    result=0
    while second:
        bit=second & -second;index=bit.bit_length()-1
        result^=first << (256*index);second^=bit
    return result


def report():
    single=single_report()
    assert single==json.loads((ROOT/'evidence/cube-favorable-reweighting.json').read_text())
    original=check_assignment_dependent_pairing()['pieces']
    generators={(g['piece'],g['hidden_axis']):g for g in single['generators']}
    dual=(1 << 240) | (1 << 255)
    endpoint=int(single['expected_truth_bits'],16)
    supports=[]
    for hidden in range(3):
        nonzero=[i for i in range(22) if
            (int(generators[i,hidden]['truth_bits'],16) & dual).bit_count()%2]
        assert len(nonzero)==1
        i=nonzero[0];piece=original[i]
        parent_axis=next(a for a in range(3) if a!=hidden)
        assert piece['ranks']==[16 if a==parent_axis else 1 for a in range(3)]
        supports.append(dict(hidden_axis=hidden,retained_axis=parent_axis,
            unique_piece=i,exponent=list(piece['exponent']),ranks=piece['ranks'],
            widths=piece['widths'],all_22_values=[
                (int(generators[j,hidden]['truth_bits'],16) & dual).bit_count()%2
                for j in range(22)]))
    assert (dual & endpoint).bit_count()%2==1
    # Independent full-vector two-block span computation, not only local products of L values.
    target=tensor_product(endpoint,endpoint);dual2=tensor_product(dual,dual)
    pivots={};favorable=0;all_cases=0;nonzero_cases=[]
    for i in range(22):
        for j in range(22):
            widths=[original[i]['widths'][a]+original[j]['widths'][a] for a in range(3)]
            for hidden in range(3):
                first=int(generators[i,hidden]['truth_bits'],16)
                second=int(generators[j,hidden]['truth_bits'],16)
                tensor=tensor_product(first,second);all_cases+=1
                value=(tensor & dual2).bit_count()%2
                if value:
                    assert i==j==supports[hidden]['unique_piece']
                    nonzero_cases.append(dict(pieces=[i,j],hidden_axis=hidden,widths=widths))
                if widths[hidden]!=max(widths):continue
                favorable+=1;assert value==0
                while tensor:
                    k=tensor.bit_length()-1
                    if k in pivots:tensor^=pivots[k]
                    else:pivots[k]=tensor;break
    residual=target
    for k in sorted(pivots,reverse=True):
        if (residual >> k)&1:residual^=pivots[k]
    assert residual and (target & dual2).bit_count()%2==1
    assert all_cases==1452 and favorable==624 and len(pivots)==110
    # Each fixed orientation remains an unrestricted positive endpoint identity.
    for hidden in range(3):
        total=0
        for i in range(22):total^=int(generators[i,hidden]['truth_bits'],16)
        assert total==endpoint and tensor_product(total,total)==target
    paths=['tools/check_cube_global_functional_profile.py',
           'tools/check_cube_favorable_reweighting.py','tools/check_two_axis_global_pairing.py',
           'tools/cube_homogeneous_gaussian.py','evidence/cube-favorable-reweighting.json']
    return dict(scope='An all-B ordinary proof uses these exact single-block support data. Fixed m5/weight1/two-word base family, common orientation, literal products, maximum hidden cut width/rank. Larger native alphabets and recombined products are not covered.',
        factoring_exponent_improved=False,endpoint_functional_value=1,
        functional_masks=[240,255],orientation_supports=supports,
        two_block=dict(all_orientations=all_cases,favorable_generators=favorable,
            favorable_span_dimension=len(pivots),endpoint_in_favorable_span=False,
            endpoint_functional_value=1,nonzero_functional_cases=nonzero_cases,
            unrestricted_endpoint_identities=3),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=report()
    if args.write:RECEIPT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(RECEIPT.read_text())==data
    print('PASS: all66 functional values,1452 two-block orientations; favorable span624 generators/dimension110 excludes E4^tensor2.')
    print('The all-B conclusion uses the stated support lemma, not finite extrapolation; no score improvement.')


if __name__=='__main__':main()
