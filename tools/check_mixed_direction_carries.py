#!/usr/bin/env python3
"""Independent point sums and the registered normal-basis AND diagnostic."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_carry_fourier_elimination import BinaryField
from check_teichmuller_cyclic_carries import POLYNOMIALS
from mixed_direction_carries import polynomial_sum_mixed
from normal_basis_and_tradeoff import normal_basis_data,binary_rank
from teichmuller_cyclic_carries import UnramifiedRing,polynomial_sum

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/mixed-direction-carry-controls.json'


def evaluate(ring,width,factor,point,directions):
    row,terms=factor
    answer=0
    for c,counts in terms:
        value=c
        for x,n,sign in zip(point,counts,directions):
            value=ring.mul(value,ring.power(x,n << ((sign*row)%width)))
        answer=ring.add(answer,value)
    return answer


def cases():
    return [dict(name='exact full domain',H=5,w=2,d=2,v=2,signs=[1,-1],factors=[]),
        dict(name='reverse distinct from forward',H=3,w=3,d=3,v=2,signs=[1,-1],
             factors=[(0,[(1,(1,1))]),(1,[(1,(3,5))])]),
        dict(name='forward distinct control',H=3,w=3,d=3,v=2,signs=[1,1],
             factors=[(0,[(1,(1,1))]),(1,[(1,(3,5))])]),
        dict(name='mixed Frobenius product',H=5,w=3,d=3,v=2,signs=[1,-1],
             factors=[(r,[(1,(0,0)),(1,(1,1))]) for r in range(3)]),
        dict(name='mixed punctures',H=4,w=3,d=3,v=2,signs=[1,-1],mask=[True,False],
             factors=[(r,[(1,(0,0)),(3,(1,1))]) for r in range(3)]),
        dict(name='all reverse',H=4,w=3,d=3,v=2,signs=[-1,-1],
             factors=[(r,[(1,(0,0)),(1,(1,1))]) for r in range(3)]),
        dict(name='even width',H=3,w=4,d=4,v=2,signs=[1,-1],
             factors=[(0,[(1,(1,2)),(3,(0,0))]),(1,[(1,(3,1))]),
                      (3,[(2,(0,0)),(1,(1,1))])]),
        dict(name='proper summation subfield',H=3,w=2,d=4,v=2,signs=[1,-1],
             factors=[(0,[(11,(1,1)),(2,(0,0))]),(1,[(5,(1,1)),(1,(0,0))])]),
        dict(name='three signs',H=3,w=2,d=2,v=3,signs=[1,-1,1],
             factors=[(r,[(1,(0,0,0)),(1,(1,1,1))]) for r in range(2)]),
        dict(name='zero factor',H=3,w=3,d=3,v=2,signs=[1,-1],factors=[(2,[])]),
        dict(name='width one',H=3,w=1,d=1,v=2,signs=[1,-1],
             factors=[(0,[(1,(1,0)),(1,(0,1))])])]


def contraction_controls():
    records=[];assignments=0;uniform=0
    for case in cases():
        H,w,d,v=(case[k] for k in ('H','w','d','v'))
        ring=UnramifiedRing(H,POLYNOMIALS[d]);signs=case['signs'];factors=case['factors']
        mask=case.get('mask',[False]*v)
        domain=[ring.teichmuller(x) for x in range(ring.field_size)]
        domain=[x for x in domain if ring.power(x,1 << w)==x]
        expected=0
        for point in product(*([x for x in domain if x or not nz] for nz in mask)):
            value=1
            for factor in factors:value=ring.mul(value,evaluate(ring,w,factor,point,signs))
            expected=ring.add(expected,value);assignments+=1
        actual,invoice=polynomial_sum_mixed(ring,w,v,factors,signs,mask)
        assert actual==expected,(case['name'],actual,expected)
        assert invoice['peak_states']<=invoice['carry_box_size']
        assert all(b<=invoice['degree_budget'] for b in invoice['carry_bounds'])
        if len(set(signs))==1:
            relabeled=[((r*signs[0])%w,terms) for r,terms in factors]
            assert actual==polynomial_sum(ring,w,w,v,relabeled,mask)
            uniform+=1
        records.append(dict(input=case,value=actual,invoice=invoice))
    assert records[1]['value']==1 and records[2]['value']==0
    return dict(records=records,direct_assignments=assignments,canonical_uniform_comparisons=uniform)


def reverse_and_rows(width,coefficients):
    support={(i+j)%width for i in range(width) for j in range(width) if coefficients[i][j]}
    choices=[]
    for shift in range(width):
        maximum=max((s-shift)%width for s in support)
        degree=(1 << (maximum//2))+(1 << ((maximum+1)//2))
        choices.append((degree,shift))
    degree,shift=min(choices)
    rows=[[] for _ in range(width)]
    for i in range(width):
        for j in range(width):
            if not coefficients[i][j]:continue
            s=(i+j-shift)%width
            a,b=s//2,(s+1)//2
            row=(i-a)%width
            assert (shift-row+b)%width==j
            rows[row].append((1,(1 << a,1 << b)))
    assert max(sum(n) for terms in rows for c,n in terms)==degree
    return rows,dict(reverse_degree=degree,forward_minimum=(1 << (width//2))+1,
                     retime_Y=shift,antidiagonal_support=sorted(support))


def and_controls():
    records=[];pairs=0
    for width in (2,3,4):
        modulus=sum((c&1) << i for i,c in enumerate(POLYNOMIALS[width]))
        field=BinaryField(modulus)
        for beta in range(1 << width):
            basis=[field.power(beta,1 << i) for i in range(width)]
            if binary_rank(basis)!=width:continue
            data=normal_basis_data(field,width,beta)
            rows,metadata=reverse_and_rows(width,data['coefficients'])
            encoded=[]
            for mask in range(1 << width):
                x=0
                for j,b in enumerate(basis):
                    if (mask>>j)&1:x ^= b
                encoded.append(x)
            for left,right in product(range(1 << width),repeat=2):
                x=encoded[left];y=field.power(encoded[right],1 << metadata['retime_Y'])
                value=0
                for row,terms in enumerate(rows):
                    for c,(a,b) in terms:
                        value ^= field.mul(field.power(x,a << row),
                                           field.power(y,b << ((-row)%width)))
                assert value==encoded[left&right]
                pairs+=1
            if width==4:assert metadata['reverse_degree']==4<metadata['forward_minimum']
            records.append(dict(width=width,beta=beta,metadata=metadata,rows=rows))
    return dict(records=records,word_pairs=pairs,
                scope='Fixed registered small panel only; no asymptotic reverse-direction AND compression established')


def report():
    sources=['tools/mixed_direction_carries.py','tools/check_mixed_direction_carries.py',
             'tools/normal_basis_and_tradeoff.py','tools/teichmuller_cyclic_carries.py',
             'tools/check_carry_fourier_elimination.py','tools/check_teichmuller_cyclic_carries.py']
    ring=UnramifiedRing(2,POLYNOMIALS[2]);rejected=[]
    for name,call in [('sign',lambda:polynomial_sum_mixed(ring,2,1,[],[0])),
                      ('row',lambda:polynomial_sum_mixed(ring,2,1,[(2,[(1,(1,))])],[1])),
                      ('negative exponent',lambda:polynomial_sum_mixed(ring,2,1,[(0,[(1,(-1,))])],[1])),
                      ('mask',lambda:polynomial_sum_mixed(ring,2,1,[],[1],[True,False]))]:
        try:call()
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    return json.loads(json.dumps(dict(
        scope='Polynomial mixed-direction Frobenius sums and a direct AND diagnostic; scalar inverse/trace-projector closure and complete integer arithmetic remain open.',
        contractions=contraction_controls(),AND=and_controls(),rejected=rejected,
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],queries=len(result['contractions']['records']),
        direct_assignments=result['contractions']['direct_assignments'],
        canonical_comparisons=result['contractions']['canonical_uniform_comparisons'],
        AND_bases=len(result['AND']['records']),AND_word_pairs=result['AND']['word_pairs'])))


if __name__=='__main__':main()
