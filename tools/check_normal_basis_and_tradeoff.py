#!/usr/bin/env python3
"""Complete fixed small normal-basis controls and subthreshold support tests."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_carry_fourier_elimination import BinaryField
from normal_basis_and_tradeoff import (binary_rank,normal_basis_data,
    synchronized_and_rows,evaluate_rows)

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/normal-basis-and-tradeoff.json'
MODELS={2:0b111,3:0b1011,4:0b10011}


def support_controls(width,t):
    columns=width//t
    below=1 << (t//2)
    missing=t//2
    tested=0;bilinear=0;offsets=set()
    modulus=(1 << width)-1
    powers={1 << i:i for i in range(width)}
    for counts in product(range(below+1),repeat=2*columns):
        if sum(counts)>below:continue
        for row in range(t):
            tested+=1
            a=sum(counts[j] << (row+t*j) for j in range(columns))
            b=sum(counts[columns+j] << (row+t*j) for j in range(columns))
            a=0 if a==0 else 1+(a-1)%modulus
            b=0 if b==0 else 1+(b-1)%modulus
            if a in powers and b in powers:
                bilinear+=1
                offsets.add((powers[a]-powers[b])%t)
    assert missing not in offsets,(width,t,missing,offsets)
    return dict(width=width,t=t,degree_cutoff=below,tested_row_monomials=tested,
                bilinear_monomials=bilinear,realized_offsets=sorted(offsets),
                missing_antipodal_offset=missing)


def report():
    bases=[];supports=[];pairs=0;row_pairs=0
    for width,modulus in MODELS.items():
        field=BinaryField(modulus)
        divisors=[t for t in range(1,width+1) if width%t==0]
        supports.extend(support_controls(width,t) for t in divisors)
        for beta in range(1 << width):
            basis=[field.power(beta,1 << j) for j in range(width)]
            if binary_rank(basis)!=width:continue
            data=normal_basis_data(field,width,beta)
            encoding={}
            for mask in range(1 << width):
                value=0
                for j,x in enumerate(basis):
                    if (mask>>j)&1:value^=x
                encoding[mask]=value
            assert len(set(encoding.values()))==1 << width
            rows=[]
            for t in divisors:
                representation,invoice=synchronized_and_rows(width,t,data['coefficients'])
                rows.append(dict(invoice=invoice,rows=representation))
                for left,right in product(range(1 << width),repeat=2):
                    expected=encoding[left&right]
                    actual=evaluate_rows(field,width,t,representation,encoding[left],encoding[right])
                    assert actual==expected
                    row_pairs+=1
            for left,right in product(range(1 << width),repeat=2):
                x,y=encoding[left],encoding[right]
                value=0
                for i,row in enumerate(data['coefficients']):
                    for j,c in enumerate(row):
                        if c:value^=field.mul(field.power(x,1 << i),field.power(y,1 << j))
                assert value==encoding[left&right]
                pairs+=1
            bases.append(dict(width=width,modulus=modulus,data=data,representations=rows))
    sources=['tools/normal_basis_and_tradeoff.py','tools/check_normal_basis_and_tradeoff.py',
             'tools/check_carry_fourier_elimination.py']
    return json.loads(json.dumps(dict(
        scope='Exact direct normal-basis AND degree/relative-column tradeoff; no lower bound for factoring or arbitrary structured contraction.',
        bases=bases,subthreshold_support=supports,
        counts=dict(normal_bases=len(bases),bilinear_word_pairs=pairs,
                    synchronized_word_pairs=row_pairs,
                    subthreshold_row_monomials=sum(x['tested_row_monomials'] for x in supports)),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],counts=result['counts'])))


if __name__=='__main__':main()
