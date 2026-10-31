#!/usr/bin/env python3
"""Literal controls for finite representatives of all regular local frames."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import hashlib
import json
from pathlib import Path
import random

from f4_lattice_linear import MUL
from formal_row_flags import dimensions,normal_form,row_catalogue
from search_binary_polynomial_flags import field_terms,lower_jet_frames


def add(a,b):
    result=dict(a)
    for power,value in b.items():
        result[power]=result.get(power,0)^value
        if not result[power]:del result[power]
    return result


def multiply(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():
            value=result.get(i+j,0)^MUL[x][y]
            if value:result[i+j]=value
            elif i+j in result:del result[i+j]
    return result


def matrix_product(a,b):
    size=len(a);result=[[{} for _ in range(size)] for _ in range(size)]
    for i,j,k in product(range(size),repeat=3):
        result[i][j]=add(result[i][j],multiply(a[i][k],b[k][j]))
    return result


def identity(size=4):
    return [[{0:1} if i==j else {} for j in range(size)] for i in range(size)]


def condition(matrix,targets):
    return all(not entry or min(entry)>=targets[i] for i,row in enumerate(matrix) for entry in row)


def key(matrix):
    return tuple(tuple(tuple(sorted(entry.items())) for entry in row) for row in matrix)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args();rng=random.Random(311700)
    normalized=conditions=0
    for trial in range(100):
        targets=tuple(sorted(rng.randrange(9) for _ in range(4)))
        frame=identity()
        for i in range(4):
            for j in range(i):
                frame[i][j]={k:value for k in range(1,targets[i]-targets[j])
                             if (value:=rng.randrange(4))}
        stabilizer=identity()
        for _ in range(12):
            i,j=rng.sample(range(4),2)
            minimum=max(1,targets[i]-targets[j])
            entry={k:value for k in range(minimum,minimum+3) if (value:=rng.randrange(4))}
            elementary=identity();elementary[i][j]=entry
            stabilizer=matrix_product(elementary,stabilizer)
        combined=matrix_product(stabilizer,frame)
        assert normal_form(combined,targets)==frame
        assert normal_form(combined,targets,precision=max(1,max(targets)-min(targets))+4)==frame
        normalized+=1
        # Test the defining row-condition equivalence on matrices with an
        # intentional exact-boundary failure and on matrices satisfying it.
        for bad in (False,True):
            source=[[{targets[i]+rng.randrange(3):rng.randrange(1,4)}
                     for j in range(4)] for i in range(4)]
            if bad:
                candidates=[i for i in range(4) if targets[i]>0]
                if not candidates:continue
                i=rng.choice(candidates);source[i][0][targets[i]-1]=1
            assert condition(source,targets)==condition(matrix_product(stabilizer,source),targets)
            conditions+=1
    targets=(0,1,1,2);classes=Counter()
    for frame in lower_jet_frames(1,4):
        matrix=[[dict(field_terms(entry)) for entry in row] for row in frame]
        classes[key(normal_form(matrix,targets))]+=1
    assert len(classes)==4 and set(classes.values())=={1024}
    finite=0;targets=(0,1,2,3)
    rows=[tuple(row_catalogue(targets,i)) for i in range(4)]
    assert [len(row) for row in rows]==[4**d for d in dimensions(targets)]
    for frame in product(*rows):
        matrix=[[dict(field_terms(entry)) for entry in row] for row in frame]
        assert normal_form(matrix,targets)==matrix
        finite+=1
    assert finite==256
    result={'status':'Formal row-frame classification controls passed.',
        'random_polynomial_stabilizer_normalizations':normalized,
        'literal_row_condition_equivalences':conditions,
        'complete_first_jet_inputs':4096,'distinct_canonical_classes':4,
        'uniform_class_multiplicity':1024,'complete_small_cell_representatives':finite,
        'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in ['formal_row_flags.py',Path(__file__).name]}}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
