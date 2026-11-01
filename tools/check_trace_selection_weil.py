#!/usr/bin/env python3
"""Literal trace characters, formal AS classes, exceptional fibres and density."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_carry_fourier_elimination import BinaryField
from trace_artin_schreier import (field_trace,degree,reduce_polynomial,
                                  reduce_system,slice_witness)

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/trace-selection-weil-controls.json'
MODELS={2:0b111,3:0b1011,4:0b10011}


def evaluate(field,polynomial,point):
    answer=0
    for c,e in polynomial:
        value=c
        for x,k in zip(point,e):value=field.mul(value,field.power(x,k))
        answer ^= value
    return answer


def domain(field,width):
    result=[x for x in range(field.size) if field.power(x,1 << width)==x]
    assert len(result)==1 << width
    return result


def character_control(field,width,variables,polynomial):
    canonical,constant=reduce_polynomial(field,width,variables,polynomial)
    assert canonical
    witness=slice_witness(field,width,variables,canonical)
    values=domain(field,width);q=1 << width
    axis=witness['axis'];M=witness['odd_degree'];K=witness['clearing_frobenius']
    exceptional=witness['exceptional_polynomial'];exception_degree=witness['exceptional_degree']
    total=0;bad=0;fibres=[];point_checks=0
    for others in product(values,repeat=variables-1):
        coefficients={}
        for c,e in canonical:
            value=c
            for x,k in zip(others,e[:axis]+e[axis+1:]):value=field.mul(value,field.power(x,k))
            coefficients[e[axis]]=coefficients.get(e[axis],0) ^ value
        fibre=[(c,(j,)) for j,c in sorted(coefficients.items()) if c]
        reduced_fibre,fibre_constant=reduce_polynomial(field,width,1,fibre)
        leading=next((c for c,e in reduced_fibre if e[0]==M),0)
        exceptional_value=evaluate(field,exceptional,others)
        assert field.power(leading,1 << K)==exceptional_value
        fibre_sum=0
        for y in values:
            point=others[:axis]+(y,)+others[axis:]
            original_bit=field_trace(field,width,evaluate(field,polynomial,point))
            canonical_bit=field_trace(field,width,evaluate(field,canonical,point)) ^ constant
            assert original_bit==canonical_bit
            fibre_sum+=(-1)**original_bit
            point_checks+=1
        is_bad=exceptional_value==0
        if is_bad:bad+=1
        else:
            assert degree(reduced_fibre)==M and M%2==1
            assert fibre_sum**2 <= (M-1)**2*q
        total+=fibre_sum
        fibres.append(dict(others=others,exceptional_value=exceptional_value,
                           odd_leading_coefficient=leading,bad=is_bad,character_sum=fibre_sum))
    assert bad*q<=exception_degree*q**(variables-1)
    allowance=exception_degree*q**(variables-1)
    excess=max(0,abs(total)-allowance)
    assert excess**2 <= (M-1)**2*q**(2*variables-1)
    D=degree(canonical)
    excess_general=max(0,abs(total)-D*D*q**(variables-1))
    assert excess_general**2 <= (D-1)**2*q**(2*variables-1)
    return dict(width=width,ambient_degree=field.degree,variables=variables,
                original=polynomial,canonical=canonical,constant_trace=constant,
                witness=witness,character_sum=total,bad_fibres=bad,fibres=fibres,
                pointwise_reduction_checks=point_checks)


def character_panel():
    records=[]
    for width,modulus in MODELS.items():
        field=BinaryField(modulus);a=2
        cases=[(1,[(1,(1,))]),(1,[(1,(3,))]),(1,[(1,(5,)),(a,(1,))]),
               (1,[(a,(8,)),(1,(2,))]),
               (2,[(1,(1,1))]),(2,[(1,(1,2)),(1,(2,1))]),
               (2,[(1,(2,3)),(1,(3,2))]),
               (3,[(1,(1,1,0)),(1,(0,0,3))])]
        for variables,polynomial in cases:
            records.append(character_control(field,width,variables,polynomial))
    field=BinaryField(MODELS[4]);sub=next(x for x in domain(field,2) if x>1)
    records.append(character_control(field,2,2,[(sub,(1,2)),(1,(2,1))]))
    special=next(r for r in records if r['width']==2 and r['ambient_degree']==2
                 and r['variables']==2 and r['original']==[(1,(1,2)),(1,(2,1))])
    assert special['bad_fibres']==4 and special['character_sum']==16
    return records


def system_controls():
    records=[];points_checked=0;fourier_checks=0
    for width in (2,3):
        field=BinaryField(MODELS[width]);a=2
        f=[(a,(1,1)),(1,(3,0))]
        square=[(field.mul(c,c),tuple(2*k for k in e)) for c,e in f]
        coboundary=f+square
        systems=[('nonlinear redundancy',[(f,1),(square,1),(coboundary,0)]),
                 ('nonlinear contradiction',[(f,1),(square,0)]),
                 ('constant target outside trace',[([(1,(0,0))],field_trace(field,width,1))]),
                 ('constant contradiction',[([(1,(0,0))],1 ^ field_trace(field,width,1))]),
                 ('two formal classes',[(f,0), ([(1,(1,2)),(1,(0,1))],1)]),
                 ('small-field functional relation',[([(1,(1,2)),(1,(2,1))],1)])]
        values=domain(field,width);q=1 << width
        for name,predicates in systems:
            reduced,invoice=reduce_system(field,width,2,predicates)
            count=0;truth=[]
            for point in product(values,repeat=2):
                original=all(field_trace(field,width,evaluate(field,f,point))==b for f,b in predicates)
                after=False if reduced is None else all(
                    field_trace(field,width,evaluate(field,f,point))==b for f,b in reduced)
                assert original==after
                truth.append(original);count+=original;points_checked+=1
            if name=='nonlinear redundancy':assert invoice['formal_rank']==1
            if name in ('nonlinear contradiction','constant contradiction'):assert reduced is None
            fourier=None
            if reduced is not None:
                rank=len(reduced);fourier=0
                for mask in range(1 << rank):
                    for point in product(values,repeat=2):
                        bit=0
                        for i,(f,b) in enumerate(reduced):
                            if (mask >> i)&1:bit ^= field_trace(field,width,evaluate(field,f,point)) ^ b
                        fourier+=(-1)**bit
                assert fourier==(1 << rank)*count
                fourier_checks+=1
                if count:
                    D=max(1,max((degree(f) for f,b in reduced),default=0))
                    for h in range(width//2+1):
                        if count*(1 << h)<=q*q:
                            assert (1 << rank)*D*(D+1)>=1 << h
            records.append(dict(name=name,width=width,predicates=predicates,
                reduced=reduced,invoice=invoice,solutions=count,fourier_numerator=fourier))
    field=BinaryField(MODELS[3]);predicates=[([(4,(1,))],1), ([(7,(1,))],0)]
    solutions=[x for x in range(8) if all(field_trace(field,3,evaluate(field,f,(x,)))==b
                                         for f,b in predicates)]
    assert solutions==[3,6]
    old_false=dict(width=3,variables=1,degree=1,predicates=predicates,solutions=solutions,
                   old_lhs=1,old_rhs=2,old_bound_rejected=True)
    return dict(records=records,pointwise_equivalences=points_checked,
                integer_fourier_checks=fourier_checks,withdrawn_ANF_counterexample=old_false)


def rejection_controls():
    field=BinaryField(MODELS[4]);rejected=[]
    for name,call in [('subfield',lambda:reduce_polynomial(field,2,1,[(2,(1,))])),
                      ('negative exponent',lambda:reduce_polynomial(field,4,1,[(1,(-1,))])),
                      ('target',lambda:reduce_system(field,4,1,[([],2)]))]:
        try:call()
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    return rejected


def rectangle_controls():
    records=[];polynomials=0
    for width,A,B in ((2,[0,1],[0,1,2]),(3,[1,2,3],[0,1]),
                      (2,[0,1,2],[1,2,3])):
        field=BinaryField(MODELS[width]);threshold=min(len(A),len(B))
        exponents=[(i,j) for i in range(threshold) for j in range(threshold-i)]
        vanishing=0
        for coefficients in product(range(1 << width),repeat=len(exponents)):
            polynomial=[(c,e) for c,e in zip(coefficients,exponents) if c]
            zeros=all(evaluate(field,polynomial,(x,y))==0 for x,y in product(A,B))
            if zeros:
                vanishing+=1
                assert not polynomial
            polynomials+=1
        # The degree threshold is sharp for this root-count argument.
        axis=0 if len(A)<=len(B) else 1
        root_set=A if axis==0 else B
        boundary={(0,0):1}
        for a in root_set:
            updated={}
            for e,c in boundary.items():
                shifted=list(e);shifted[axis]+=1;shifted=tuple(shifted)
                updated[shifted]=updated.get(shifted,0) ^ c
                updated[e]=updated.get(e,0) ^ field.mul(c,a)
            boundary={e:c for e,c in updated.items() if c}
        boundary=[(c,e) for e,c in sorted(boundary.items())]
        assert boundary and degree(boundary)==threshold
        assert all(evaluate(field,boundary,(x,y))==0 for x,y in product(A,B))
        records.append(dict(width=width,A=A,B=B,strict_degree_cutoff=threshold,
                            tested_polynomials=(1 << width)**len(exponents),
                            vanishing_polynomials=vanishing,sharp_boundary=boundary))
    return dict(records=records,tested_polynomials=polynomials)


def report():
    sources=['tools/trace_artin_schreier.py','tools/check_trace_selection_weil.py',
             'tools/check_carry_fourier_elimination.py']
    characters=character_panel();systems=system_controls();rectangles=rectangle_controls()
    return json.loads(json.dumps(dict(
        scope='Formal AS reduction, exceptional-fibre character bound and thin pointwise trace selection; no runtime lower bound or factoring improvement.',
        characters=characters,systems=systems,rectangles=rectangles,
        rejected=rejection_controls(),
        counts=dict(character_cases=len(characters),
                    pointwise_reduction_checks=sum(r['pointwise_reduction_checks'] for r in characters),
                    fibres=sum(len(r['fibres']) for r in characters),
                    system_cases=len(systems['records']),
                    pointwise_system_equivalences=systems['pointwise_equivalences'],
                    integer_fourier_checks=systems['integer_fourier_checks'],
                    subthreshold_rectangle_polynomials=rectangles['tested_polynomials']),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],counts=result['counts'],rejected=result['rejected'])))


if __name__=='__main__':main()
