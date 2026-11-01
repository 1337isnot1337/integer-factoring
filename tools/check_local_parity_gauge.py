#!/usr/bin/env python3
"""Replay constructive local parity gauges and exact inconsistency witnesses."""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
from random import Random

import local_parity_gauge as gauge
import cube_f4_bivariate_smith as smith
from check_cube_f4_bivariate_smith import seed_tensor

PAIRS=tuple((i,i^7) for i in range(4))
ROOT=Path(__file__).resolve().parents[1]


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def serialized(tensor):
    return [{'word':i,'terms':[[*e,c] for e,c in sorted(polynomial.items())]} for i,polynomial in enumerate(tensor)]


def load_tensor(path):
    record=json.loads(path.read_text())
    return tuple({tuple(term[:3]):term[3] for term in row['terms']} for row in record['entries']),record


def f4_inverse(matrix):
    a,b=matrix[0];c,d=matrix[1]
    determinant=gauge.F4_MUL[a][d]^gauge.F4_MUL[b][c]
    if not determinant:raise ValueError('Singular F4 matrix.')
    scale=gauge.F4_INV[determinant]
    return [[gauge.F4_MUL[scale][d],gauge.F4_MUL[scale][b]],
            [gauge.F4_MUL[scale][c],gauge.F4_MUL[scale][a]]]


def random_unit(rng):
    while True:
        matrix=[[rng.randrange(4) for _ in range(2)] for _ in range(2)]
        try:f4_inverse(matrix)
        except ValueError:continue
        return matrix


def extension_control(rng):
    # The paired global operator N satisfies N^2=N. Its first factor has
    # determinant alpha, whose quadratic does not split over F4.
    parameters=[1,2,0,0,0,0,0,0]
    monomials=((6,0,0),(0,6,0),(0,0,6));anchor=(2,2,2)
    result=[{anchor:1} if all((mask>>i&1)==(mask>>j&1) for i,j in PAIRS) else {} for mask in range(256)]
    for exponent in monomials:
        vector=[rng.randrange(4) for _ in range(256)]
        for mask in range(256):
            value=vector[mask] if not mask.bit_count()%2 else 0
            for p,(i,j) in enumerate(PAIRS):
                b,c=parameters[2*p:2*p+2]
                value^=gauge.F4_MUL[b if not(mask>>i&1) else c][vector[mask^(1<<i)]]
                value^=gauge.F4_MUL[c if not(mask>>j&1) else b][vector[mask^(1<<j)]]
            if value:
                result[mask][exponent]=value
                result[mask][anchor]=result[mask].get(anchor,0)^value
                if not result[mask][anchor]:del result[mask][anchor]
    return tuple(result),parameters


def check():
    rng=Random(2026110121);seed=seed_tensor();positive=[];equation_count=0
    for paired in (False,True):
        for trial in range(8):
            maps=[None]*8
            if paired:
                for i,j in PAIRS:
                    maps[i]=random_unit(rng);maps[j]=gauge.transpose(f4_inverse(maps[i]))
            else:maps=[random_unit(rng) for _ in range(8)]
            tensor=gauge.apply_maps(seed,maps,gauge.F4_MUL);pairs=PAIRS if paired else None
            record=gauge.find_gauge(tensor,pairs);assert record['equivalent']
            checked=gauge.verify_gauge(tensor,record,pairs)
            equation_count+=record['tested_equations']
            positive.append(dict(trial=trial,paired=paired,tensor_sha256=digest(serialized(tensor)),
                                 record=record,verification=checked))
    single=tuple({(1,):1} if mask==1 else {} for mask in range(4))
    single_record=gauge.find_gauge(single);assert single_record['equivalent']
    single_checked=gauge.verify_gauge(single,single_record)
    equation_count+=single_record['tested_equations']

    extension,parameters=extension_control(rng)
    extension_record=gauge.find_gauge(extension,PAIRS)
    assert extension_record['equivalent'] and extension_record['rank']==8
    assert extension_record['solution']==parameters
    assert not any(gauge.F4_MUL[x][x]^x==2 for x in range(4))
    assert any(c not in gauge.F4_EMBED for matrix in extension_record['gauges_f16'] for row in matrix for c in row)
    extension_checked=gauge.verify_gauge(extension,extension_record,PAIRS)
    general_extension=gauge.find_gauge(extension);assert general_extension['equivalent']
    general_extension_checked=gauge.verify_gauge(extension,general_extension)
    equation_count+=extension_record['tested_equations']+general_extension['tested_equations']

    reference_path=ROOT/'evidence/local-parity-nongauge-reference.json'
    reference,fixture=load_tensor(reference_path);negative=[]
    for paired in (False,True):
        pairs=PAIRS if paired else None;record=gauge.find_gauge(reference,pairs)
        assert not record['equivalent']
        checked=gauge.verify_obstruction(reference,record['obstruction'],pairs)
        negative.append(dict(paired=paired,record=record,verification=checked))
        equation_count+=record['tested_equations']
    assert len(negative[1]['record']['obstruction'])==2
    profiles=[]
    for axis in range(3):
        record=smith.profile_tensor(reference,axis);assert record['orders']==fixture['expected_orders'][axis]
        profiles.append(dict(axis=axis,orders=record['orders'],verification=smith.verify_tensor(reference,axis,record)))

    rejected=0
    bad=copy.deepcopy(negative[1]['record']['obstruction']);bad.pop()
    try:gauge.verify_obstruction(reference,bad,PAIRS)
    except AssertionError:rejected+=1
    else:raise AssertionError('A deleted certificate equation was accepted.')
    bad=copy.deepcopy(extension_record);bad['gauges_f16']=[[[1,0],[0,1]] for _ in range(8)]
    try:gauge.verify_gauge(extension,bad,PAIRS)
    except AssertionError:rejected+=1
    else:raise AssertionError('Identity maps were accepted on a tensor with odd support.')
    invalid=copy.deepcopy(seed);invalid[0][(2,2,2)]=2
    try:gauge.find_gauge(invalid,PAIRS)
    except ValueError:rejected+=1
    else:raise AssertionError('An invalid EQ endpoint was accepted.')
    try:gauge.find_gauge(tuple({} for _ in range(256)))
    except ValueError:rejected+=1
    else:raise AssertionError('The excluded zero tensor was accepted.')
    paths=[Path(__file__),Path(gauge.__file__),Path(smith.__file__),Path(smith.__file__).with_suffix('.cpp'),
           Path(__file__).with_name('check_cube_f4_bivariate_smith.py'),reference_path]
    return dict(schema=1,scope='Exact constant local parity-gauge equivalence over coefficient-field extensions; no new factoring exponent',
       coefficient_field='F4',gauge_field='F16',positive_controls=positive,
       odd_support_positive_control=dict(record=single_record,verification=single_checked),
       necessary_quadratic_extension_control=dict(tensor_sha256=digest(serialized(extension)),unique_paired_parameters=parameters,
        record=extension_record,verification=extension_checked,general_verification=general_extension_checked),
       nongauge_reference=dict(tensor_sha256=digest(serialized(reference)),obstructions=negative,complete_generic_profiles=profiles),
       positive_gauge_count=len(positive)+3,independent_obstruction_count=2,
       linear_equations_checked=equation_count,rejected_controls=rejected,
       source_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'evidence/local-parity-gauge.json')
    args=parser.parse_args();record=check();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:record[key] for key in ('positive_gauge_count','independent_obstruction_count','linear_equations_checked','rejected_controls')}))


if __name__=='__main__':main()
