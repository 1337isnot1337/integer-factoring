#!/usr/bin/env python3
"""All-field uniform-marginal entry witnesses for both pencil endpoints."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from check_two_place_determinant import determinant
from joint_energy_conjugations import N, binary_product, model
from two_place_f4_census import build_input


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    text,construction = build_input()
    tokens = text.decode().split()
    entries = [tuple(map(int,tokens[1+4*i:1+4*(i+1)])) for i in range(625)]
    parameters = list(map(F,tokens[2501:]))
    witnesses = json.loads((root/'evidence/pencil-uniform-transport-witnesses.json').read_text())
    assert [w['place'] for w in witnesses] == ['zero','infinity']
    records = []
    for witness in witnesses:
        side = int(witness['place']=='infinity')
        marginals = [[F(0)]*5 for _ in range(4)]
        total,cost = F(0),F(0)
        for index,power,probability in witness['plan']:
            probability = F(probability)
            assert probability>0 and 0<=index<625 and power>=0
            low,high,k,l = entries[index]
            assert (low|high)>>power&1
            # q has degree <=11. Neither perturbation can change this
            # coefficient over any characteristic-two field.
            assert not (k and 0<=power<=11)
            assert not (l and 4<=power<=15)
            total += probability
            cost += probability*(-power if side else power)
            for leg in range(4):
                marginals[leg][index//5**leg%5] += probability
        assert total==1 and marginals==[[F(1,5)]*5 for _ in range(4)]
        expected = F(-66,5) if side else F(14,5)
        assert cost == F(witness['value']) == expected
        # The fixed natural gauges attain the upper bound for every q.
        native_side = 0 if side else 1
        def height(index):
            labels = [index//5**leg%5 for leg in range(4)]
            return sum(F(0) if label==0 else parameters[native_side*16+4*leg+label-1]
                       for leg,label in enumerate(labels))
        source_min,perturb_min = None,None
        for index,(low,high,k,l) in enumerate(entries):
            packed = low|high
            if packed:
                valuation = -(packed.bit_length()-1) if side else (packed&-packed).bit_length()-1
                value = valuation+height(index)
                source_min = value if source_min is None else min(source_min,value)
            for shift,present in ((0,k),(4,l)):
                if present:
                    value = height(index)-(11+shift) if side else height(index)+shift
                    perturb_min = value if perturb_min is None else min(perturb_min,value)
        minimum = F(-8) if side else F(8)
        assert source_min==minimum and perturb_min>=minimum
        assert entries[0]==(1<<8,0,0,0) and height(0)==0
        average = sum(parameters[native_side*16:(native_side+1)*16])/5
        assert minimum-average == cost
        records.append({'place':witness['place'],'uniform_transport_cost':str(cost),
                        'positive_support_entries':len(witness['plan']),
                        'uniform_marginals':[[str(v) for v in m] for m in marginals],
                        'natural_source_minimum':str(source_min),
                        'natural_perturbation_lower_bound':str(perturb_min),
                        'natural_uniform_height':str(average),'all_field_optimum':str(cost)})
    data = model()
    generator = sum(1<<(N*r+c) for pair in ((4,20),(7,11),(8,16)) for r in pair for c in pair)
    assert binary_product(generator,generator)==0
    coefficients = [sum(((data['operator'][r][c][0]>>degree)&1)<<(N*r+c)
                        for r in range(N) for c in range(N)) for degree in range(19)]
    changes = [binary_product(generator,a)^binary_product(a,generator) for a in coefficients]
    assert all(binary_product(binary_product(generator,a),generator)==0 for a in coefficients)
    assert changes[12]==changes[8]^changes[16]
    assert all(not v for i,v in enumerate(changes) if i not in (8,12,16))
    # [G,A]=t^8*(1+t)^4*(K+t^4*L); this and GAG=G^2=0
    # prove rational conjugation and determinant invariance for every q.
    row_determinant = determinant([[pair[0] for pair in row] for row in data['operator']])
    p0 = sum(1<<i for i in (0,4,6,10,12,14,16))
    assert row_determinant==p0<<314
    earlier = json.loads((root/'evidence/two-place-determinant.json').read_text())
    assert earlier['minimum_opposite_determinant_degree']==78
    for name,digest in earlier['source_sha256'].items():
        assert sha256((root/name).read_bytes()).hexdigest()==digest
    # The opposite-cut degree >=78 is the earlier all-field theorem;
    # the present replay freshly computes the row determinant and transport.
    gaps = {'zero_minimum_vs_infinity_mean':F(14,5)-F(78,25),
            'infinity_minimum_vs_zero_mean':F(-66,5)+F(314,25)}
    assert gaps=={'zero_minimum_vs_infinity_mean':F(-8,25),
                  'infinity_minimum_vs_zero_mean':F(-16,25)}
    names = ['tools/check_pencil_uniform_transport.py','tools/check_two_place_determinant.py',
             'tools/two_place_f4_census.py','tools/joint_energy_conjugations.py',
             'tools/search_f4_recentered_crossings.py','tools/tensor_pascal_box_crossings.py',
             'evidence/pencil-uniform-transport-witnesses.json','evidence/two-place-determinant.json']
    result = {'status':'pass','transport_certificates':records,'fresh_row_determinant':row_determinant,
              'row_determinant_formula':'t^314*(1+t^4+t^6+t^10+t^12+t^14+t^16)',
              'rational_conjugation_identity':'A_q=(I+fG)A(I+fG), f=q/(t^8*(1+t)^4)',
              'construction':construction,'opposite_degree_bound_reused':78,
              'all_field_endpoint_mean_gaps':{k:str(v) for k,v in gaps.items()},
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'All characteristic-two fields, full degree<=11 pencil and arbitrary '
              'physical gauges for the two uniform-type endpoint scalings. '
              'No all-angle support exclusion or stronger factoring exponent.'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'pass','all_field_optima':[r['all_field_optimum'] for r in records],
                      'all_field_endpoint_mean_gaps':result['all_field_endpoint_mean_gaps']}))


if __name__ == '__main__':
    main()
