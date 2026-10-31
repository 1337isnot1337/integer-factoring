#!/usr/bin/env python3
"""Exact nonlinear-charge polytope and typed synthesis controls."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_typed_charge_lattice_synthesis import chi_order,coefficient_table_control
from check_reflection_quartic_tiles import two_copy_control
from nonlinear_symbol_charge_synthesis import SymbolChargeLattice,cube_lipschitz,four_symbol_coordinates
from typed_charge_lattice_synthesis import natural_heights,CompressedFullAlphabetFiber
from reflection_dual_rail_gaussian import tensor


def polytope_certificate():
    values=tensor(4)['values']
    basis=((2,2,4,0,0,0),(2,-2,0,0,0,0),(2,2,0,0,0,0),
           (0,0,0,2,2,4),(0,0,0,2,-2,0),(0,0,0,2,2,0))
    facets={};checks=0;orders=Counter();nonzero=0
    for index,(value,high) in enumerate(values):
        assert not high
        if not value:continue
        nonzero+=1;valuation=chi_order(value);orders[valuation]+=1
        x,y,u,v=(index//4**i%4 for i in range(4))
        row=[0]*6
        for letter,offset,sign in ((x,0,-1),(v,0,1),(y,3,-1),(u,3,1)):
            if letter:row[offset+letter-1]+=sign
        coefficients=tuple(sum(a*b for a,b in zip(row,direction)) for direction in basis)
        for coefficient in coefficients:
            assert abs(coefficient)<=valuation;checks+=2
        if valuation and all(abs(c)==valuation for c in coefficients):
            signs=tuple(c//valuation for c in coefficients)
            facets[signs]={'entry':[x,y,u,v],'chi_order':valuation,'constraint':row}
    assert set(facets)==set(product((-1,1),repeat=6))
    assert all(record['chi_order']==2 for record in facets.values())
    for direction in basis:
        x=(0,)+direction[:3];y=(0,)+direction[3:]
        assert cube_lipschitz(x)+cube_lipschitz(y)==2
        coordinates=four_symbol_coordinates(x)+four_symbol_coordinates(y)
        assert sum(map(abs,coordinates))==1
    return {'nonzero_original_coefficients':nonzero,'valuation_histogram':dict(orders),
            'signed_vertices':12,'exact_vertex_inequality_checks':checks,
            'all_sign_facets':64,'basis':basis,
            'facet_witnesses':[{'signs':signs,**record} for signs,record in sorted(facets.items())]}


def nonlinear_coefficient_controls():
    checked=neighbors=0;records=[]
    for M in (4,8):
        m=M//2;q=1<<m;source=tensor(M)['values']
        valuations={index:chi_order(value) for index,(value,high) in enumerate(source) if value}
        functions=[([2*(state.bit_count()%2) for state in range(q)],[0]*q),
                   ([state.bit_count() for state in range(q)],[(state.bit_count()%2) for state in range(q)]),
                   ([min(state.bit_count(),(state^((1<<m)-1)).bit_count()) for state in range(q)],
                    [state.bit_count()%2 for state in range(q)])]
        for fx,fy in functions:
            assert cube_lipschitz(fx)+cube_lipschitz(fy)<=2
            for index,valuation in valuations.items():
                x,y,u,v=(index//q**i%q for i in range(4))
                assert valuation+fx[x]+fy[y]-fy[u]-fx[v]>=0;checked+=1
        for x,y,i,j in product(range(q),range(q),range(m),range(m)):
            v=x^(1<<i);u=y^(1<<j)
            index=x+q*y+q*q*u+q**3*v
            assert valuations[index]==2;neighbors+=1
        records.append({'interval':M,'alphabet':q,'nonzero_coefficients':len(valuations)})
    return checked,neighbors,records


def verify():
    rng=random.Random(927361);polytope=polytope_certificate()
    coefficients,neighbors,intervals=nonlinear_coefficient_controls()
    heights,_=natural_heights(4)
    choices=[(([0,2,2,0],[0,0,0,0]),False),
             (([0,2,2,0],[0,0,0,0]),True),
             (([0,1,1,0],[0,1,1,0]),True),
             (([0,1,1,1],[0,1,1,0]),False),
             (([0,Fraction(1,100),Fraction(4,100),Fraction(1,10)],
               [0,Fraction(1,4),Fraction(1,2),Fraction(4,5)]),False)]
    records=[];affine=typed=profiles=0
    for weights,boundary in choices:
        lattice=SymbolChargeLattice(4,heights,raise_cap=3,charge_weights=weights,negative_boundary=boundary)
        low=rng.getrandbits(lattice.variables);high=rng.getrandbits(lattice.variables)
        affine+=coefficient_table_control(lattice,low,high)
        values=lattice.seed(low,high)
        typed+=two_copy_control(lattice.as_seed(values))
        actual=lattice.profiles(values);profiles+=2
        floor=tuple(max(0,order-3) for order in lattice.orders)
        fiber=CompressedFullAlphabetFiber(lattice,floor)
        pivots=fiber.solve(lattice.orders);assert pivots is not None
        witness=fiber.sample(pivots,rng)
        lattice.check_witness(witness,lattice.orders);profiles+=2
        typed+=two_copy_control(lattice.as_seed(witness))
        records.append({'potentials':[[str(v) for v in row] for row in lattice.symbol_weights],
            'lipschitz_constants':[str(cube_lipschitz(row)) for row in lattice.symbol_weights],
            'negative_boundary':boundary,'variables':lattice.variables,'common_pole':lattice.pole,
            'numerator_degree':lattice.degree,'effective_degree':lattice.degree-lattice.pole,
            'random_profiles_infinity_then_zero':actual,'fiber_variables':fiber.variables})
    assert records[-1]['variables']==6264 and records[-1]['common_pole']==0
    rejected=0
    try:SymbolChargeLattice(4,heights,charge_weights=([0,3,3,0],[0,0,0,0]))
    except ValueError:rejected+=1
    else:raise AssertionError('An inadmissible nonlinear charge was accepted.')
    return {'status':'Exact nonlinear symbol-charge controls passed','polytope':polytope,
            'all_coefficient_charge_checks':coefficients,'exact_neighbor_valuations':neighbors,
            'literal_intervals':intervals,'affine_coefficient_comparisons':affine,
            'two_copy_typed_coefficients':typed,'complete_local_profiles':profiles,
            'rejected_invalid_charge':rejected,'affine_cases':records}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('evidence/nonlinear-symbol-charges.json'))
    args=p.parse_args();result=verify()
    names=[Path(__file__).name,'nonlinear_symbol_charge_synthesis.py','typed_charge_lattice_synthesis.py',
           'check_typed_charge_lattice_synthesis.py','check_reflection_quartic_tiles.py',
           'reflection_dual_rail_gaussian.py','full_alphabet_two_cut_lattices.py',
           'native_matrix_smith_f4.py','native_matrix_smith_f4.cpp']
    result['source_sha256']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    note=Path(__file__).resolve().parents[1]/'research'/'nonlinear-symbol-charge-lattices.md'
    result['note_sha256']=sha256(note.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in
        ('polytope','literal_intervals','affine_cases','source_sha256','note_sha256')}))


if __name__=='__main__':main()
