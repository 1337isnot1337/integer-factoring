#!/usr/bin/env python3
"""Direct finite checks of collected inverse and packed-system compilers."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from teichmuller_cyclic_carries import UnramifiedRing
from precision_rational_carries import rational_sum, packed_system_count, conjugate_factor
from check_teichmuller_cyclic_carries import evaluate, POLYNOMIALS

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/precision-rational-carry-controls.json'


def domain(ring,width):
    values=[ring.teichmuller(x) for x in range(ring.field_size)]
    values=[x for x in values if ring.power(x,1 << width)==x]
    assert len(values)==1 << width
    return values


def rational_controls():
    results=[];assignments=0;conjugations=0
    for L in range(1,5):
        ring=UnramifiedRing(L,POLYNOMIALS[2]);omega=ring.teichmuller(2)
        for t in (1,2):
            columns=2//t;zero=(0,)*columns;first=(1,)+(0,)*(columns-1)
            last=(0,)*(columns-1)+(1,)
            numerator=(0,[(1,first),(ring.integer(2),zero)])
            denominators=[(0,[(1,first),(1,zero)]),
                          (t-1,[(omega,first),(1,zero)]),
                          (0,[(1,last),(omega,zero)])]
            for factor in [numerator]+denominators:
                for shift in range(4):
                    changed=conjugate_factor(ring,2,t,1,factor,shift)
                    for x in domain(ring,2):
                        expected=ring.frobenius(evaluate(ring,1,t,factor,[x]),shift)
                        assert evaluate(ring,1,t,changed,[x])==expected
                        conjugations+=1
            for punctured in (False,True):
                expected=0;values=[]
                for x in domain(ring,2):
                    if punctured and not x: continue
                    value=evaluate(ring,1,t,numerator,[x]);den_values=[]
                    for factor in denominators:
                        den=evaluate(ring,1,t,factor,[x]);den_values.append(den)
                        value=ring.mul(value,ring.unit_inverse(den)) if ring.mod_two(den) else 0
                    expected=ring.add(expected,value);assignments+=1
                    values.append(dict(input=x,denominators=den_values,value=value))
                got,invoice=rational_sum(ring,2,t,1,[numerator],denominators,punctured)
                assert got==expected
                assert invoice['queries']==L+1
                assert invoice['max_degree_budget']==1+((2*L+1)*2-1)*3
                results.append(dict(precision=L,row_modulus=t,punctured=punctured,
                                    numerator=numerator,denominators=denominators,
                                    direct_values=values,result=got,invoice=invoice))
    ring=UnramifiedRing(3,POLYNOMIALS[2])
    for denominators,expected in (([],4),([(0,[(ring.integer(2),(0,))])],0)):
        got,invoice=rational_sum(ring,2,2,1,[],denominators)
        assert got==expected
        results.append(dict(precision=3,numerators=[],denominators=denominators,
                            result=got,invoice=invoice,constant_control=True))
    return dict(queries=len(results),direct_assignments=assignments,
                literal_frobenius_factor_comparisons=conjugations,cases=results)


def root_controls():
    results=[];assignments=0
    configurations=[(1,1,1,1),(2,1,2,2),(3,1,1,3),
                    (4,2,1,2),(4,2,2,2)]
    for L in (1,2,3):
        for degree,width,variables,count in configurations:
            ring=UnramifiedRing(L,POLYNOMIALS[degree])
            if variables==1:
                equations=[[(1,(1,)),(ring.integer(i%2),(0,))] for i in range(count)]
            elif width==1:
                equations=[[(1,(1,0))],[(1,(0,1))]]
            else:
                omega=ring.power(ring.teichmuller(2),5)
                assert ring.frobenius(omega,2)==omega and omega not in (0,1)
                equations=[[(1,(1,0)),(omega,(0,1)),(1,(0,0))],
                           [(1,(0,2)),(1,(0,1))]]
            for punctured in (False,True):
                expected=0;values=[]
                choices=[x for x in domain(ring,width) if x or not punctured]
                for point in product(choices,repeat=variables):
                    outputs=[evaluate(ring,variables,width,(0,eq),point) for eq in equations]
                    is_root=all(ring.mod_two(value)==0 for value in outputs)
                    expected+=is_root;assignments+=1
                    values.append(dict(point=point,equation_values=outputs,is_root=is_root))
                got,invoice=packed_system_count(ring,width,variables,equations,punctured)
                assert got==expected%ring.modulus
                assert invoice['variables']==variables and invoice['precision']==L
                assert invoice['queries']==L+1 and invoice['relative_columns']==1
                results.append(dict(precision=L,coefficient_degree=degree,width=width,
                                    variables=variables,punctured=punctured,equations=equations,
                                    direct_count=expected,result=got,invoice=invoice,direct_values=values))
    # A precision exceeding the entire four-element one-variable domain.
    ring=UnramifiedRing(5,POLYNOMIALS[4])
    equations=[[(1,(1,)),(1,(0,))],[(1,(2,)),(1,(0,))]]
    got,invoice=packed_system_count(ring,2,1,equations)
    expected=sum(all(not ring.mod_two(evaluate(ring,1,2,(0,eq),[x])) for eq in equations)
                 for x in domain(ring,2))
    assert got==expected==1 and ring.precision>2
    assignments+=4
    results.append(dict(precision=5,coefficient_degree=4,width=2,variables=1,
                        equations=equations,result=got,direct_count=expected,invoice=invoice,
                        exact_count_precision_exceeds_domain_bits=True))
    for equations,punctured,expected in (([],False,4),([],True,1),([[]],False,4),
                                         ([[(1,(0,0))]],False,0)):
        binary=UnramifiedRing(3,POLYNOMIALS[2])
        got,invoice=packed_system_count(binary,1,2,equations,punctured)
        assert got==expected
        results.append(dict(precision=3,coefficient_degree=2,width=1,variables=2,
                            equations=equations,punctured=punctured,result=got,
                            direct_count=expected,invoice=invoice,constant_control=True))
    try:
        packed_system_count(ring,2,1,[[(ring.generator,(1,))]])
    except ValueError:
        pass
    else:
        raise AssertionError('Coefficient outside the subfield ring was accepted')
    return dict(queries=len(results),direct_assignments=assignments,cases=results)


def report():
    sources=['tools/precision_rational_carries.py','tools/check_precision_rational_carries.py',
             'tools/teichmuller_cyclic_carries.py','tools/check_teichmuller_cyclic_carries.py']
    result=dict(scope='Implemented collected total unit inverses and packed residue systems; unchanged summation variables and precision, charged extension/degree; no complete factoring gain.',
                rational=rational_controls(),packed_systems=root_controls(),
                source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})
    return json.loads(json.dumps(result))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],rational_queries=result['rational']['queries'],
                          rational_direct_assignments=result['rational']['direct_assignments'],
                          frobenius_factor_checks=result['rational']['literal_frobenius_factor_comparisons'],
                          packed_queries=result['packed_systems']['queries'],
                          packed_direct_assignments=result['packed_systems']['direct_assignments'])))


if __name__=='__main__':main()
