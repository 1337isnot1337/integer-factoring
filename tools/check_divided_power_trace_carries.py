#!/usr/bin/env python3
"""Independent point sums and integer algebra controls for trace lifting."""
import argparse
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

from teichmuller_cyclic_carries import UnramifiedRing
from divided_power_trace_carries import (DividedPowerRing,residue_bits,
    acquire_trace_factors,trace_constrained_sum,affine_trace_sum)
from check_teichmuller_cyclic_carries import evaluate,POLYNOMIALS
from check_carry_fourier_elimination import BinaryField
from square_zero_trace_contraction import square_zero_constrained_trace

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/divided-power-trace-controls.json'


def domain(ring,width):
    values=[ring.teichmuller(a) for a in range(ring.field_size)]
    result=[a for a in values if ring.power(a,1 << width)==a]
    assert len(result)==1 << width
    return result


def field_model(ring):
    return BinaryField(sum((c & 1) << i for i,c in enumerate(ring.polynomial)))


def field_trace_predicate(ring,field,width,terms,bit,point):
    value=0
    for coefficient,exponents in terms:
        term=residue_bits(ring,coefficient)
        for x,e in zip(point,exponents):term=field.mul(term,field.power(residue_bits(ring,x),e))
        value^=term
    trace=0
    for j in range(width):trace^=field.power(value,1 << j)
    assert trace in (0,1)
    return trace==bit


def evaluate_gamma(algebra,variables,t,factor,point):
    row,terms=factor;columns=len(terms[0][1])//variables if terms else 1
    answer=0
    for coefficient,counts in terms:
        value=coefficient
        for i,x in enumerate(point):
            exponent=sum(counts[i*columns+j] << (row+t*j) for j in range(columns))
            value=algebra.mul(value,algebra.base.power(x,exponent))
        answer=algebra.add(answer,value)
    return answer


def algebra_controls():
    base=UnramifiedRing(2,POLYNOMIALS[2]);associativity=0;exponentials=0
    for rank,limit in ((1,3),(2,3)):
        algebra=DividedPowerRing(base,rank,limit)
        basis=[algebra.tag(i) for i in range(algebra.dimension)]
        for a,b,c in product(basis,repeat=3):
            assert algebra.mul(algebra.mul(a,b),c)==algebra.mul(a,algebra.mul(b,c))
            associativity+=1
    algebra=DividedPowerRing(base,1,3)
    gamma=algebra.tag(1)
    assert algebra.mul(gamma,gamma)==algebra.tag(2,base.integer(2))!=0
    assert algebra.frobenius(gamma)!=algebra.power(gamma,2)
    try:algebra.unit_inverse(gamma)
    except ValueError:pass
    else:raise AssertionError('A nilpotent was inverted')
    def exponential(a):return sum(algebra.tag(k,base.power(a,k)) for k in range(4))
    for a,b in product(range(16),repeat=2):
        assert algebra.mul(exponential(a),exponential(b))==exponential(base.add(a,b))
        exponentials+=1
    return dict(associativity_basis_triples=associativity,exponential_identities=exponentials,
                divided_square_coefficient=2,lifted_frobenius_is_not_squaring=True,
                general_nilpotent_inverse_rejected=True)


def cases():
    result=[]
    def add(name,H,d,w,t,v,predicates,mask=False,numerators=(),denominators=()):
        result.append(dict(name=name,H=H,d=d,w=w,t=t,v=v,predicates=predicates,
                           mask=mask,numerators=numerators,denominators=denominators))
    linear=[([(1,(1,))],1)]
    add('empty domain count',3,2,2,2,1,[])
    add('precision-one trace',1,2,2,2,1,linear)
    add('two relative columns',2,2,2,1,1,linear)
    add('exact one-word count',3,2,2,2,1,linear)
    ring=UnramifiedRing(3,POLYNOMIALS[2]);omega=ring.teichmuller(2)
    add('two predicates one word',3,2,2,2,1,
        [([(1,(1,))],0),([(1,(2,)),(omega,(0,))],1)])
    ring=UnramifiedRing(2,POLYNOMIALS[2]);omega=ring.teichmuller(2)
    coupled=[([(1,(1,1)),(omega,(1,0))],1),([(1,(1,0)),(1,(0,2))],0)]
    add('coupled trace polynomials',2,2,2,2,2,coupled)
    add('mixed punctured domains',2,2,2,2,2,coupled,[True,False])
    ring=UnramifiedRing(2,POLYNOMIALS[4]);omega=ring.power(ring.teichmuller(2),5)
    coefficient=ring.add(omega,ring.mul(ring.integer(2),ring.generator))
    assert ring.frobenius(coefficient,2)!=coefficient
    assert ring.mod_two(ring.frobenius(coefficient,2))==ring.mod_two(coefficient)
    add('higher digits outside subring',2,4,2,2,1,[([(coefficient,(1,))],1)])
    add('nonlinear cubic trace',4,3,3,3,1,[([(1,(3,))],0)])
    add('selected scalar pole',2,2,2,2,1,[([(1,(1,))],0)],
        denominators=[(0,[(1,(1,)),(1,(0,))])])
    add('constant nonunit denominator',2,2,2,2,1,[([(1,(1,))],0)],
        denominators=[(0,[(2,(0,))])])
    add('inconsistent polynomial predicates',2,2,2,2,1,
        [([(1,(1,))],0),([(1,(1,))],1)])
    add('even-width target constant',3,2,2,2,1,[([(1,(0,))],1)])
    return result


def run_cases():
    records=[];assignments=0;pointwise=0
    for case in cases():
        name,H,d,w,t,v=(case[k] for k in ('name','H','d','w','t','v'))
        ring=UnramifiedRing(H,POLYNOMIALS[d]);field=field_model(ring)
        predicates=case['predicates'];mask=case['mask']
        flags=mask if isinstance(mask,list) else [mask]*v
        choices=[[x for x in domain(ring,w) if x or not flag] for flag in flags]
        algebra,factors,metadata=acquire_trace_factors(ring,w,t,v,predicates)
        expected=0;integer_count=0;points=[]
        for point in product(*choices):
            accepted=all(field_trace_predicate(ring,field,w,terms,bit,point)
                         for terms,bit in predicates)
            encoded=1
            for factor in factors:encoded=algebra.mul(encoded,evaluate_gamma(algebra,v,t,factor,point))
            assert algebra.indicator_functional(encoded)==int(accepted)
            pointwise+=1
            value=1
            for factor in case['numerators']:value=ring.mul(value,evaluate(ring,v,t,factor,point))
            for factor in case['denominators']:
                den=evaluate(ring,v,t,factor,point)
                value=ring.mul(value,ring.unit_inverse(den)) if ring.mod_two(den) else 0
            if not accepted:value=0
            else:integer_count+=1
            expected=ring.add(expected,value);assignments+=1
            points.append(dict(point=point,accepted=accepted,weighted_value=value))
        got,invoice=trace_constrained_sum(ring,w,t,v,predicates,
                                         case['numerators'],case['denominators'],mask)
        assert got==expected,(name,got,expected)
        assert invoice['algebra_dimension']==(2*H)**len(predicates)
        assert invoice['actual_backend']['max_degree_budget']<=invoice['max_degree_budget']
        exact=(H>v*w and not case['numerators'] and not case['denominators'])
        if exact:assert got==integer_count
        records.append(dict(input=case,result=got,invoice=invoice,points=points,
                            direct_accepted_count=integer_count,exact_count_recovered=exact))
    return dict(queries=len(records),direct_assignments=assignments,
                pointwise_indicator_controls=pointwise,cases=records)


def affine_controls():
    records=[]
    for H in (1,3):
        ring=UnramifiedRing(H,POLYNOMIALS[2]);omega=ring.teichmuller(2)
        systems=[([([1],1),([1],1),([0],0)],2,1),
                 ([([1],0),([omega],1)],1,2),
                 ([([1],0),([1],1)],0,None)]
        for rows,count,rank in systems:
            got,invoice=affine_trace_sum(ring,2,2,1,rows)
            assert got==count%ring.modulus and invoice['rank']==rank
            if H==1:
                field=field_model(ring)
                old,old_invoice=square_zero_constrained_trace(field,2,2,1,[],[],rows)
                assert old==got
            records.append(dict(precision=H,rows=rows,result=got,invoice=invoice))
    return records


def report():
    sources=['tools/divided_power_trace_carries.py','tools/check_divided_power_trace_carries.py',
             'tools/precision_rational_carries.py','tools/teichmuller_cyclic_carries.py',
             'tools/trace_constraint_contraction.py','tools/square_zero_trace_contraction.py',
             'tools/check_teichmuller_cyclic_carries.py','tools/check_carry_fourier_elimination.py']
    result=dict(scope='Growing-precision polynomial trace predicates through a charged divided-power algebra; exact finite point counts, no integer arithmetic representation or factoring improvement.',
                algebra=algebra_controls(),contractions=run_cases(),affine=affine_controls(),
                source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})
    return json.loads(json.dumps(result))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],algebra=result['algebra'],
                          contraction_queries=result['contractions']['queries'],
                          direct_assignments=result['contractions']['direct_assignments'],
                          affine_controls=len(result['affine']))))


if __name__=='__main__':main()
