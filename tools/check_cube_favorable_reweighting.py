#!/usr/bin/env python3
"""Exact finite endpoint-span test for favorable cube orientations."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_two_axis_global_pairing import (AXES,pullback,
                                         check_assignment_dependent_pairing)

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/cube-favorable-reweighting.json'


def span(generators,target):
    pivots={}
    for i,tensor in enumerate(generators):
        coefficient=1 << i
        while tensor:
            k=tensor.bit_length()-1
            if k in pivots:
                tensor^=pivots[k][0];coefficient^=pivots[k][1]
            else:
                pivots[k]=(tensor,coefficient);break
    def reduce(tensor):
        coefficient=0
        for k in sorted(pivots,reverse=True):
            if (tensor >> k)&1:
                tensor^=pivots[k][0];coefficient^=pivots[k][1]
        return tensor,coefficient
    residual,coefficient=reduce(target)
    if residual:
        k=residual.bit_length()-1
        dual=sum(((reduce(1 << i)[0] >> k)&1) << i for i in range(256))
        assert all((dual & v).bit_count()%2==0 for v in generators)
        assert (dual & target).bit_count()%2==1
        solution=None
    else:
        dual=None;solution=[i for i in range(len(generators)) if (coefficient >> i)&1]
        total=0
        for i in solution:total^=generators[i]
        assert total==target
    return dict(generators=len(generators),dimension=len(pivots),endpoint_in_span=not residual,
                dual_truth_bits=None if dual is None else hex(dual),solution=solution)


def report():
    source=check_assignment_dependent_pairing()
    endpoint=int(source['expected_truth_bits'],16)
    permutations={axes:[sum(b << j for j,b in enumerate(
        pullback(tuple((mask >> i)&1 for i in range(8)),axes)))
        for mask in range(256)] for axes in AXES}
    generators=[]
    for i,piece in enumerate(source['pieces']):
        tensor=int(piece['original_truth_bits'],16)
        for hidden in range(3):
            axes=tuple(a for a in range(3) if a!=hidden)
            changed=sum(((tensor >> old)&1) << m
                        for m,old in enumerate(permutations[axes]))
            generators.append(dict(piece=i,exponent=list(piece['exponent']),hidden_axis=hidden,
                widths=piece['widths'],favorable=piece['widths'][hidden]==max(piece['widths']),
                first_choice=hidden==piece['selected_hidden_axis'],truth_bits=hex(changed)))
    modes={}
    for mode in ('first_choice','favorable','unrestricted'):
        selected=[i for i,g in enumerate(generators) if mode=='unrestricted' or g[mode]]
        result=span([int(generators[i]['truth_bits'],16) for i in selected],endpoint)
        result['generator_indices']=selected;modes[mode]=result
    assert [(modes[m]['generators'],modes[m]['dimension'],modes[m]['endpoint_in_span'])
            for m in modes]==[(22,8,False),(30,8,False),(66,14,True)]
    # An explicit two-entry functional is easier to audit than elimination.
    dual=(1 << 240) | (1 << 255)
    assert modes['favorable']['dual_truth_bits']==hex(dual)
    for g in generators:
        if g['favorable']:assert (int(g['truth_bits'],16) & dual).bit_count()%2==0
    assert (endpoint & dual).bit_count()%2==1
    paths=['tools/check_cube_favorable_reweighting.py',
           'tools/check_two_axis_global_pairing.py','tools/cube_homogeneous_gaussian.py']
    return dict(scope='Only the fixed m5/weight1/two-word Gaussian coefficient family. All scalar reweightings of its favorable single-block orientations fail, over every characteristic-two extension. No asymptotic or multi-block exclusion.',
        factoring_exponent_improved=False,word_length=5,weight=1,alphabet=[1,2],
        generators=generators,modes=modes,expected_truth_bits=hex(endpoint),
        explicit_dual=dict(masks=[240,255],endpoint_values=[(endpoint >> m)&1 for m in (240,255)]),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=report()
    if args.write:RECEIPT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(RECEIPT.read_text())==data
    print('PASS: favorable 30-generator span has dimension8 and excludes E4; unrestricted dimension14 contains it.')
    print('Exact dual: T[240]+T[255]. Fixed single-block family only; complete score unchanged.')


if __name__=='__main__':main()
