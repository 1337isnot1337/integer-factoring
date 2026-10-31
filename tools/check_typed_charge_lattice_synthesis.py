#!/usr/bin/env python3
"""Exact controls for complete-entry typed charge lattices."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

from check_rectangular_defect_support import pdivmod
from check_reflection_quartic_tiles import two_copy_control,exact_fock_profile,polynomial_product,endpoint
from reflection_dual_rail_gaussian import tensor,laurent_kernel
from reflection_quartic_tiles import full_rank_certificate
from typed_charge_lattice_synthesis import (TypedChargeLattice,CompressedFullAlphabetFiber,
    natural_heights,chi_power,cut_location,matrix_multiply,polynomial_degrees)


def chi_order(value):
    assert value
    result=0
    while True:
        quotient,remainder=pdivmod(value,3)
        if remainder:return result
        value=quotient;result+=1


def edge_controls():
    result=0
    for M in (4,8,16,32):
        m=M//2;kernel=laurent_kernel(M)
        choices=[tuple(Fraction(i,m-1) for _ in range(2) for i in range(m)),
                 (Fraction(1),)*(2*m),(Fraction(2),)*m+(Fraction(0),)*m,
                 tuple(Fraction((-1)**i,2) for i in range(m))+
                 tuple(Fraction(3*(-1)**(i+1),2) for i in range(m))]
        for weights in choices:
            charge=weights[:m]+weights[m:]+tuple(-v for v in weights[m:])+tuple(-v for v in weights[:m])
            assert max(map(abs,weights[:m]))+max(map(abs,weights[m:]))<=2
            for i in range(2*M):
                for j in range(i+1,2*M):
                    entry=kernel[i][j]
                    if not entry:continue
                    low=min(entry);value=sum(1<<(k-low) for k in entry)
                    assert chi_order(value)+charge[i]+charge[j]>=0
                    result+=1
    return result


def coefficient_table_control(lattice,low,high):
    values=lattice.seed(low,high);checked=0
    for side,(left,right) in enumerate(lattice.bases):
        h=lattice.heights[4*side:4*side+4];scalar=lattice.scalar_shifts[side]
        for plane,assignment in enumerate((low,high)):
            matrix=[[0]*lattice.size for _ in range(lattice.size)]
            for entry,pair in enumerate(values):
                row,col=cut_location(lattice.labels[entry],lattice.q,side)
                weight=sum(v[b] for v,b in zip(h,lattice.labels[entry]))
                for power in polynomial_degrees(pair[plane]):
                    exponent=scalar+weight+(-power if side==0 else power)
                    assert exponent>=0;matrix[row][col]^=1<<exponent
            transformed=matrix_multiply(matrix_multiply(left,matrix),right)
            for row in range(lattice.size):
                index=side*lattice.size+row;expected={}
                for (power,equation,rhs),col in zip(lattice.equations[index],lattice.equation_columns[index]):
                    value=((equation&assignment).bit_count()&1)^(rhs if plane==0 else 0)
                    expected[power,col]=value
                for power in range(lattice.caps[index]):
                    for col in range(lattice.size):
                        assert ((transformed[row][col]>>power)&1)==expected.get((power,col),0)
                        checked+=1
    return checked


def verify():
    rng=random.Random(4398531);edges=edge_controls();heights,_=natural_heights(4)
    choices=[([0,0,0,0],False),([0,1,0,1],False),([1,1,1,1],False),
             ([2,2,0,0],False),([1,-1,1,-1],False),
             ([1,1,1,1],True),([2,2,0,0],True),
             ([Fraction(1,2)]*4,True)]
    records=[];table_coefficients=typed_entries=base_coefficients=profile_controls=0
    independent_profiles=negative_controls=alternate_frames=0
    for case,(weights,boundary) in enumerate(choices):
        h=[list(row) for row in heights]
        if case%3==1:
            h=[[value+rng.randrange(-2,3) for value in row] for row in h]
        lattice=TypedChargeLattice(4,h,raise_cap=4,charge_weights=weights,negative_boundary=boundary)
        for value,energy in zip(lattice.source,lattice.charges):
            if value:
                valuation=chi_order(value)
                assert valuation+energy>=0
                if valuation+energy==0 and energy:
                    assert energy<0
                base_coefficients+=1
        low=rng.getrandbits(lattice.variables);high=rng.getrandbits(lattice.variables)
        table_coefficients+=coefficient_table_control(lattice,low,high)
        arbitrary=lattice.seed(low,high);seed=lattice.as_seed(arbitrary)
        typed_entries+=two_copy_control(seed)
        rank_witness=full_rank_certificate(seed)
        arbitrary_profiles=lattice.profiles(arbitrary);profile_controls+=2
        floors=tuple(max(0,order-4) for order in lattice.orders)
        fiber=CompressedFullAlphabetFiber(lattice,floors)
        pivots=fiber.solve(lattice.orders);assert pivots is not None
        witness=fiber.sample(pivots,rng)
        profiles=lattice.check_witness(witness,lattice.orders);profile_controls+=2
        typed_entries+=two_copy_control(lattice.as_seed(witness))
        exact_sides=[]
        # An independent rational DVR elimination checks dense binary entries,
        # including the restored common endpoint denominator.
        if case in (1,5):
            sparse_low=0
            for entry in rng.sample(range(256),12):
                sparse_low|=1<<(lattice.offsets[entry]+rng.randrange(lattice.widths[entry]))
            binary=lattice.seed(sparse_low);binary_seed=lattice.as_seed(binary)
            actual=lattice.profiles(binary)
            for side in (0,1):
                original_side=1-side
                independent=exact_fock_profile(binary_seed,original_side,h[4*side:4*side+4])
                assert tuple(independent)==actual[side]
                exact_sides.append(side);independent_profiles+=1
            if case==5:
                alternate=TypedChargeLattice(4,h,raise_cap=3,charge_weights=weights,
                    negative_boundary=boundary,reference_values=binary)
                assert alternate.check_witness(binary,alternate.orders)==actual
                table_coefficients+=coefficient_table_control(alternate,low,high)
                alternate_floors=tuple(max(0,v-2) for v in alternate.orders)
                alternate_fiber=CompressedFullAlphabetFiber(alternate,alternate_floors)
                alternate_pivots=alternate_fiber.solve(alternate.orders)
                assert alternate_pivots is not None
                altered=alternate_fiber.sample(alternate_pivots,rng)
                alternate.check_witness(altered,alternate.orders)
                typed_entries+=two_copy_control(alternate.as_seed(altered))
                profile_controls+=4;alternate_frames+=2
        records.append({'weights':[str(v) for v in lattice.charge_weights],
            'negative_boundary':boundary,'variables':lattice.variables,'pole':lattice.pole,
            'numerator_degree':lattice.degree,'effective_degree':lattice.degree-lattice.pole,
            'fiber_variables':fiber.variables,'reference_solution_dimension':fiber.variables-len(pivots),
            'random_profiles_infinity_then_zero':arbitrary_profiles,
            'feasible_profiles_infinity_then_zero':profiles,'rank_specializations':rank_witness,
            'independent_rational_sides':exact_sides})
    assert records[0]['variables']==records[1]['variables']==records[2]['variables']==6144
    assert records[5]['variables']==6237
    # Equality at positive charge can create a non-crossing typed endpoint.
    base=tensor(4)['values'];q=4
    at=lambda x,y,u,v:x+q*y+q*q*u+q**3*v
    positive=at(1,1,0,0);negative=at(0,0,1,1)
    assert chi_order(base[negative][0])==2
    bad_numerator=[(value[0]<<0,value[1]) for value in base]
    # Common denominator chi^2: add one only to the positive-charge entry.
    from weyl_lattice_filter_synthesis import polynomial_multiply
    bad_numerator=[(polynomial_multiply(a,chi_power(2)),polynomial_multiply(b,chi_power(2)))
                   for a,b in bad_numerator]
    a,b=bad_numerator[positive];bad_numerator[positive]=(a^1,b)
    product=polynomial_product(bad_numerator[positive],bad_numerator[negative])
    quotient=[]
    for plane in product:
        value,remainder=pdivmod(plane,chi_power(4));assert not remainder;quotient.append(value)
    assert endpoint(quotient)==1;negative_controls+=1
    try:TypedChargeLattice(4,heights,charge_weights=[2,2,1,1])
    except ValueError:negative_controls+=1
    else:raise AssertionError('Invalid edge budget was accepted.')
    return {'status':'Exact typed charge-lattice synthesis controls passed',
        'kernel_edge_controls':edges,'original_coefficient_charge_controls':base_coefficients,
        'independent_affine_coefficient_comparisons':table_coefficients,
        'two_copy_typed_coefficients':typed_entries,'full_local_profiles':profile_controls,
        'independent_rational_profiles':independent_profiles,'alternate_reference_frames':alternate_frames,
        'negative_controls':negative_controls,
        'cases':records}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('evidence/typed-charge-lattice-synthesis.json'))
    args=parser.parse_args();result=verify()
    names=[Path(__file__).name,'typed_charge_lattice_synthesis.py','full_alphabet_two_cut_lattices.py',
           'weyl_lattice_filter_synthesis.py','reflection_quartic_tiles.py',
           'check_reflection_quartic_tiles.py','reflection_dual_rail_gaussian.py',
           'native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','finite_precision_smith_frames.py']
    result['source_sha256']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    note=Path(__file__).resolve().parents[1]/'research'/'typed-charge-lattice-synthesis.md'
    if note.exists():result['note_sha256']=sha256(note.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('cases','source_sha256','note_sha256')}))


if __name__=='__main__':main()
