#!/usr/bin/env python3
"""Exact small characters, point counts, precision and acquisition controls."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from artin_hasse_trace_carries import (splitting_coefficients,
    acquire_character_factors, trace_character_sum, residue_bits)
from check_carry_fourier_elimination import BinaryField
from check_teichmuller_cyclic_carries import POLYNOMIALS, evaluate
from divided_power_trace_carries import trace_constrained_sum
from precision_rational_carries import conjugate_factor, factor_degree
from teichmuller_cyclic_carries import UnramifiedRing
from witt_artin_hasse import artin_hasse_coefficients

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/artin-hasse-trace-controls.json'


def field_model(ring):
    return BinaryField(sum((c & 1) << i for i, c in enumerate(ring.polynomial)))


def field_trace(field, width, value):
    result = 0
    for j in range(width):
        result ^= field.power(value, 1 << j)
    assert result in (0, 1)
    return result


def predicate_bit(ring, field, width, polynomial, point):
    result = 0
    for c, exponents in polynomial:
        term = residue_bits(ring, c)
        for x, e in zip(point, exponents):
            term = field.mul(term, field.power(residue_bits(ring, x), e))
        result ^= term
    return field_trace(field, width, result)


def character_controls():
    records = []
    for P in (2, 3, 6, 9):
        pi, theta = splitting_coefficients(P)
        _, higher = splitting_coefficients(P+2)
        assert theta == tuple(c % (1 << P) for c in higher[:P])
        assert pi % 4 == 2
        assert all(c % (1 << n) == 0 for n, c in enumerate(theta))
        assert all(c.denominator % 2 for c in artin_hasse_coefficients(P-1))
        for d in (1, 2, 3, 4):
            ring = UnramifiedRing(P, POLYNOMIALS[d])
            field = field_model(ring)
            for residue in range(1 << d):
                x = ring.teichmuller(residue)
                value = 1
                for j in range(d):
                    value = ring.mul(value, ring.evaluate(theta, ring.frobenius(x, j)))
                bit = field_trace(field, d, residue)
                assert value == ring.integer((-1)**bit)
                records.append(dict(P=P, d=d, residue=residue, value=value, trace=bit))
    ring = UnramifiedRing(5, POLYNOMIALS[2])
    pi, theta = splitting_coefficients(5)
    # 1+1 is not Teichmuller. Its residue trace is zero, but direct splitting
    # of this sum gives 25 rather than 1. Monomialwise splitting gives 1.
    unsplit = ring.power(ring.evaluate(theta, ring.integer(2)), 2)
    assert unsplit == 25 != 1
    proper = []
    ring = UnramifiedRing(5, POLYNOMIALS[4])
    for residue in range(16):
        x = ring.teichmuller(residue)
        if ring.frobenius(x, 2) == x:
            value = ring.mul(ring.evaluate(theta, x),
                             ring.evaluate(theta, ring.frobenius(x)))
            bit = field_trace(field_model(ring), 2, residue)
            assert value == ring.integer((-1)**bit)
            proper.append(dict(residue=residue, value=value, trace=bit))
    assert len(proper) == 4
    return dict(full_field_cases=records, proper_subfield_cases=proper,
                unsplit_counterexample=dict(P=5, width=2, pi=pi,
                    input_value=2, residue_trace=0, wrong_value=unsplit, correct_value=1))


def cases():
    result = []
    def add(name, H, d, w, t, v, predicates, nonzero=False, numerators=(), denominators=()):
        result.append(dict(name=name,H=H,d=d,w=w,t=t,v=v,predicates=predicates,
                           nonzero=nonzero,numerators=numerators,denominators=denominators))
    add('no constraints',3,2,2,2,1,[])
    add('exact affine count',3,2,2,2,1,[([(1,(1,))],1)])
    add('relative columns',2,2,2,1,1,[([(1,(1,))],1)])
    ring = UnramifiedRing(5,POLYNOMIALS[2]);omega=ring.teichmuller(2)
    add('two independent predicates',3,2,2,2,1,
        [([(1,(1,))],0), ([(omega,(1,))],1)])
    add('nonlinear cubic',4,3,3,3,1,[([(1,(3,))],0)])
    # All supplied coefficient encodings use the case's INPUT precision P=H+r.
    ring = UnramifiedRing(4,POLYNOMIALS[2]);omega=ring.teichmuller(2)
    coupled=[([(1,(1,1)),(omega,(1,0))],1), ([(1,(1,0)),(1,(0,2))],0)]
    add('coupled polynomial traces',2,2,2,2,2,coupled)
    add('mixed punctures',2,2,2,2,2,coupled,[True,False])
    add('exact binary nonlinear count',3,1,1,1,2,[([(1,(1,1))],1)])
    ring = UnramifiedRing(4,POLYNOMIALS[4]);sub=ring.power(ring.teichmuller(2),5)
    coefficient=ring.add(sub,ring.mul(ring.integer(2),ring.generator))
    assert ring.frobenius(coefficient,2) != coefficient
    add('higher digits outside subring',3,4,2,2,1,[([(coefficient,(1,))],1)])
    add('scalar pole',3,2,2,2,1,[([(1,(1,))],0)],
        denominators=[(0,[(1,(1,)),(1,(0,))])])
    add('nonunit denominator',2,2,2,2,1,[([(1,(1,))],0)],
        denominators=[(0,[(2,(0,))])])
    add('contradictory traces',2,2,2,2,1,
        [([(1,(1,))],0), ([(1,(1,))],1)])
    add('even width constant target',3,2,2,2,1,[([(1,(0,))],1)])
    ring=UnramifiedRing(5,POLYNOMIALS[2]);omega=ring.teichmuller(2)
    add('coefficientwise division with nonscalar weight',3,2,2,2,1,
        [([(1,(1,))],0), ([(omega,(1,))],1)],
        numerators=[(0,[(ring.generator,(0,))])])
    return result


def full_row_caps(ring, case, groups, bound):
    factors = list(case['numerators'])+[f for group in groups for f in group]
    for denominator in case['denominators']:
        orbit=[conjugate_factor(ring,case['w'],case['t'],case['v'],denominator,s)
               for s in range(ring.degree)]
        factors.extend(orbit* (2*ring.precision))
        factors.extend(orbit[1:])
    degrees=[sum(factor_degree(f) for f in factors if f[0]==row)
             for row in range(case['w'])]
    carries=[]
    for a in range(case['v']):
        rotated=sum(max([0]+[counts[a] << row for c,counts in terms if c])
                    for row,terms in factors)
        carries.append(rotated//((1 << case['w'])-1))
    assert max(degrees+[0]) <= bound
    assert max(carries+[0]) <= bound
    return dict(row_degrees=degrees,initial_carry_caps=carries,bound=bound)


def contraction_controls():
    records=[];assignments=0;point_characters=0;dp_comparisons=0
    for case in cases():
        H,d,w,t,v=(case[k] for k in ('H','d','w','t','v'))
        predicates=case['predicates'];P=H+len(predicates)
        ring=UnramifiedRing(P,POLYNOMIALS[d]);out=UnramifiedRing(H,POLYNOMIALS[d])
        field=field_model(ring)
        mask=case['nonzero']
        if isinstance(mask,bool):mask=[mask]*v
        domain=[ring.teichmuller(x) for x in range(1 << d)]
        domain=[x for x in domain if ring.frobenius(x,w)==x]
        assert len(domain)==1 << w
        groups=[acquire_character_factors(ring,w,t,v,f)[0] for f,b in predicates]
        expected=0;count=0;points=[]
        for point in product(*([x for x in domain if x or not nz] for nz in mask)):
            assignments+=1
            bits=[predicate_bit(ring,field,w,f,point) for f,b in predicates]
            for group,bit in zip(groups,bits):
                character=1
                for row,terms in group:
                    character=ring.mul(character,evaluate(ring,v,t,(row,terms),point))
                assert character==ring.integer((-1)**bit)
                point_characters+=1
            accepted=all(bit==target for bit,(f,target) in zip(bits,predicates))
            weight=1
            for factor in case['numerators']:
                weight=ring.mul(weight,evaluate(ring,v,t,factor,point))
            for factor in case['denominators']:
                value=evaluate(ring,v,t,factor,point)
                weight=ring.mul(weight,ring.unit_inverse(value) if ring.mod_two(value) else 0)
            if accepted:
                count+=1
                expected=ring.add(expected,weight)
            points.append(dict(point=point,bits=bits,accepted=accepted,weight=weight))
        expected=out.from_polynomial(ring.coefficients(expected))
        actual,invoice=trace_character_sum(ring,w,t,v,predicates,
                            case['numerators'],case['denominators'],mask)
        assert actual==expected,(case['name'],actual,expected)
        exact=H>v*w and not case['numerators'] and not case['denominators']
        if exact:assert actual==count
        def reduce_factor(factor):
            row,terms=factor
            return row,[(out.from_polynomial(ring.coefficients(c)),k) for c,k in terms]
        reduced=[([(out.from_polynomial(ring.coefficients(c)),k) for c,k in f],b)
                 for f,b in predicates]
        # Independent coefficient-algebra path, not used by the character compiler.
        dp,dp_invoice=trace_constrained_sum(out,w,t,v,reduced,
            list(map(reduce_factor,case['numerators'])),
            list(map(reduce_factor,case['denominators'])),mask)
        assert dp==actual
        dp_comparisons+=1
        caps=full_row_caps(ring,case,groups,invoice['full_row_degree_and_carry_bound']) if t==w else None
        records.append(dict(input=case,result=actual,invoice=invoice,points=points,
                            direct_count=count,exact_integer_recovery=exact,full_row_caps=caps))
    return dict(cases=records,direct_assignments=assignments,
                pointwise_characters=point_characters,dp_comparisons=dp_comparisons)


def rejection_controls():
    ring=UnramifiedRing(3,POLYNOMIALS[4]);rejected=[]
    calls=[('precision',lambda:trace_character_sum(ring,2,2,1,[([],0)]*3)),
           ('coefficient subfield',lambda:acquire_character_factors(ring,2,2,1,[(ring.generator,(1,))])),
           ('negative exponent',lambda:acquire_character_factors(ring,2,2,1,[(1,(-1,))])),
           ('target',lambda:trace_character_sum(ring,2,2,1,[([],2)])),
           ('dimensions',lambda:acquire_character_factors(ring,3,3,1,[(1,(1,))]))]
    for name,call in calls:
        try:call()
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    return rejected


def report():
    sources=['tools/artin_hasse_trace_carries.py','tools/check_artin_hasse_trace_carries.py',
             'tools/witt_artin_hasse.py','tools/precision_rational_carries.py',
             'tools/teichmuller_cyclic_carries.py','tools/divided_power_trace_carries.py',
             'tools/check_teichmuller_cyclic_carries.py','tools/check_carry_fourier_elimination.py']
    return json.loads(json.dumps(dict(
        scope='Classical Artin-Hasse binary characters plus charged Fourier projection; exact represented trace counts, no native digit evaluator or factoring improvement.',
        characters=character_controls(),contractions=contraction_controls(),
        rejected=rejection_controls(),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],full_character_cases=len(result['characters']['full_field_cases']),
        proper_subfield_cases=len(result['characters']['proper_subfield_cases']),
        queries=len(result['contractions']['cases']),direct_assignments=result['contractions']['direct_assignments'],
        pointwise_characters=result['contractions']['pointwise_characters'],
        dp_comparisons=result['contractions']['dp_comparisons'],rejected=result['rejected'])))


if __name__=='__main__':main()
