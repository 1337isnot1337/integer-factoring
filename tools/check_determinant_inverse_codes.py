#!/usr/bin/env python3
"""Exact finite acquisition controls for the determinant physical codes."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from coherent_f8_maps import dot,mul
from determinant_inverse_codes import (INV,code,character,cauchy_matrix,
                                       minor_vector,scalar_minor,serializable)

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/determinant-inverse-code-controls.json'


def derivative(vector,deleted):
    return {w ^ deleted:c for w,c in vector.items() if w & deleted==deleted}


def scalar_det_gaussian(matrix):
    matrix=[row[:] for row in matrix];result=1
    for col in range(len(matrix)):
        pivot=next((i for i in range(col,len(matrix)) if matrix[i][col]),None)
        if pivot is None:return 0
        matrix[col],matrix[pivot]=matrix[pivot],matrix[col]
        value=matrix[col][col];result=mul(result,value)
        for i in range(col+1,len(matrix)):
            ratio=mul(matrix[i][col],INV[value])
            for j in range(col,len(matrix)):matrix[i][j]^=mul(ratio,matrix[col][j])
    return result


def report():
    cases=[];counts=Counter()
    for k in (2,3):
        expected_character=character(k)
        for kind in ('diagonal','antidiagonal','cauchy'):
            data=code(k,kind);p=data['embedding'];left=data['retraction'];labels=data['labels']
            actual=Counter((x['charge'],x['height']) for x in labels)
            assert actual==Counter({(r,h):n for r,hs in expected_character.items() for h,n in hs.items()})
            by_label={(tuple(x['rows']),tuple(x['columns'])):i for i,x in enumerate(labels)}
            for index,label in enumerate(labels):
                rows,cols=label['rows'],label['columns']
                for bit in range(k*k):
                    i,j=divmod(bit,k)
                    key=tuple(a for a in rows if a!=i),tuple(b for b in cols if b!=j)
                    target=by_label.get(key)
                    expected=p[target] if i in rows and j in cols else {}
                    assert derivative(p[index],1 << bit)==expected
                    counts['embedding_single_deletions']+=1
                    got=derivative(left[index],1 << bit)
                    if got:
                        assert target is not None
                        factor=next(iter(got.values()))
                        lead=next(iter(got))
                        factor=mul(factor,INV[left[target].get(lead,0)])
                        assert got=={w:mul(factor,c) for w,c in left[target].items()}
                    counts['retraction_single_deletions']+=1
            for i,l in enumerate(left):
                for j,v in enumerate(p):
                    assert dot(l,v)==int(i==j);counts['pairing_entries']+=1
            if kind=='cauchy':
                a=cauchy_matrix(k)
                for label in labels:
                    rows,cols=label['rows'],label['columns']
                    value=scalar_minor(a,rows,cols)
                    assert scalar_det_gaussian([[a[i][j] for j in cols] for i in rows])==value!=0
                    counts['independent_cauchy_minors']+=1
                assert any(c>1 for vector in left for c in vector.values())
            cases.append(data)
        # Every squarefree derivative of both top forms is a scalar minor or zero.
        for weights in (None,cauchy_matrix(k)):
            top=minor_vector(k,list(range(k)),list(range(k)),weights)
            for deleted in range(1 << (k*k)):
                selected=[divmod(bit,k) for bit in range(k*k) if deleted >> bit & 1]
                rows={i for i,j in selected};cols={j for i,j in selected}
                got=derivative(top,deleted)
                if len(rows)!=len(selected) or len(cols)!=len(selected):assert not got
                else:
                    factor=1
                    if weights is not None:
                        for i,j in selected:factor=mul(factor,weights[i][j])
                    expected=minor_vector(k,[i for i in range(k) if i not in rows],
                        [j for j in range(k) if j not in cols],weights)
                    assert got=={w:mul(factor,c) for w,c in expected.items()}
                counts['complete_top_derivatives']+=1
    characters=[]
    for k in range(1,13):
        c=character(k);flat=Counter()
        for r,values in c.items():
            assert sum(values.values())==comb(k,r)**2
            for h,n in values.items():
                assert c[k-r].get(k**3-h)==n
                flat[h]+=n
        assert sum(h*n for h,n in flat.items())*2==k**3*comb(2*k,k)
        characters.append(dict(k=k,q=comb(2*k,k),top_height=k**3,by_charge=c))
    paths=['tools/determinant_inverse_codes.py','tools/check_determinant_inverse_codes.py',
           'tools/coherent_f8_maps.py','tools/coherent_hot_compounds.py']
    return serializable(dict(scope='Exact finite determinant-minor map acquisition and general height-character controls; no crossing pressure or improved factoring score.',
        factoring_exponent_improved=False,cases=cases,counts=dict(counts),characters=characters,
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=report()
    if args.write:RECEIPT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(RECEIPT.read_text())==data
    print(json.dumps(dict(status='PASS',cases=len(data['cases']),counts=data['counts'],
                          character_sizes=len(data['characters']),factoring_exponent_improved=False)))


if __name__=='__main__':main()
