#!/usr/bin/env python3
"""Put regular physical core CAR into coordinate bases for fast real gauges.

Vacuum-projector pivot minors are units, certified independently modulo z.
The full matrix is L K B; both L and B are regular with regular inverses.
Only the mathematically zero y-both rows and u-both columns are removed.
"""
from __future__ import annotations

import argparse
import ctypes
import json
from pathlib import Path
import time

from check_rectangular_defect_core_basis import core_basis_data,data_for
from rectangular_defect_profiles import Program,Expr,library,lcm2,quotient2,MASK,precision_sizes
from native_matrix_smith_f4 import MatrixOracle

def coordinate_matrix(data,alpha=1,beta=2,quartic_terms=None,repair_spec=None,profile_precision=None,input_spec=None,repair_kind='partial',identity_spec=None,extra_quartics=None,full_alphabet_gamma=None,extra_monomials=None,residual_monomials=None):
    start=time.monotonic();bases=[core_basis_data(data,side) for side in ('left','right')]
    extended=dict(data);extended['active_images']=dict(data['active_images']);names=[]
    for side,basis in zip(('left','right'),bases):
        prefix='basis_'+side+'_';names.append([prefix+name for name in basis['core_names']])
        extended['active_images'].update({prefix+name:value for name,value in basis['regular_images'].items()})
    if full_alphabet_gamma is None:
        if residual_monomials is not None:raise ValueError('A general residual correction requires the full physical alphabet.')
        program=Program(extended)
    else:
        from full_alphabet_quartic_profiles import FullAlphabetProgram
        program=FullAlphabetProgram(extended,full_alphabet_gamma)
    k=program.k;n=1<<k;left,right=bases
    left_names,right_names=names;left_states=1<<len(left_names);right_states=1<<len(right_names)
    diagonal=[sum(s for j,s in enumerate(data['active_orders']) if mask>>j&1) for mask in range(n)]
    def vacuum(a,core_names):
        for name in core_names:a=program.add(a,program.number(a,name))
        return a
    right_vacuum=vacuum(Expr(0),right_names)
    right_vectors={0:right_vacuum};outputs={}
    for state in range(right_states):
        if state:
            bit=state&-state;index=bit.bit_length()-1
            right_vectors[state]=program.car(right_vectors[state^bit],right_names[index]+'_create')
        initial=program.diagonal(right_vectors[state],diagonal)
        arguments=(alpha,beta,quartic_terms,initial,repair_spec,input_spec,repair_kind,identity_spec,extra_quartics,extra_monomials)
        value=program.build(*arguments,residual_monomials=residual_monomials) if residual_monomials is not None else program.build(*arguments)
        annihilated={0:value}
        for mask in range(left_states):
            if mask:
                bit=mask&-mask;index=bit.bit_length()-1
                annihilated[mask]=program.car(annihilated[mask^bit],left_names[index]+'_annihilate')
            outputs[mask,state]=vacuum(annihilated[mask],left_names)
    base=min(v.base for v in outputs.values());denominator=1
    for value in outputs.values():denominator=lcm2(denominator,value.denominator)
    output_ids=[]
    for mask in range(left_states):
        for state in range(right_states):
            value=outputs[mask,state]
            factor=quotient2(denominator,value.denominator)<<(value.base-base)
            value=program.polynomial(value,factor);output_ids.append(value.register)
    bound=max(program.degree_bounds[i] for i in output_ids)
    window,builder=precision_sizes(program,base,data['M'],profile_precision);words=builder//64;output_words=window//64
    instructions=(ctypes.c_int*(5*len(program.instructions)))(*(v for row in program.instructions for v in row))
    packed=(ctypes.c_uint64*(2*words*len(program.coefficients)))(*(value>>(64*w)&MASK for pair in program.coefficients for value in pair for w in range(words)))
    terms=[];offsets=[0]
    for car in program.cars:terms.extend(v for term in car for v in term);offsets.append(len(terms)//3)
    cc=(ctypes.c_int*len(terms))(*terms);co=(ctypes.c_int*len(offsets))(*offsets)
    dd=(ctypes.c_int*n)(*(v-min(diagonal) for v in diagonal));ii=(ctypes.c_int*len(left['vacuum_rows']))(*left['vacuum_rows']);jj=(ctypes.c_int*len(right['vacuum_columns']))(*right['vacuum_columns'])
    rr=(ctypes.c_int*len(output_ids))(*output_ids);output=(ctypes.c_uint64*(2*output_words*n*n))();shifts=(ctypes.c_int64*(n*n))();meta=(ctypes.c_int64*4)()
    lib=library(builder);I=ctypes.POINTER(ctypes.c_int);U=ctypes.POINTER(ctypes.c_uint64);L=ctypes.POINTER(ctypes.c_int64)
    lib.fock_f4_coordinates_v2.argtypes=[ctypes.c_int,ctypes.c_int,I,ctypes.c_int,U,ctypes.c_int,I,I,I,ctypes.c_int,I,ctypes.c_int,I,I,ctypes.c_int,ctypes.c_int,U,L,L]
    lib.fock_f4_coordinates_v2.restype=ctypes.c_int
    status=lib.fock_f4_coordinates_v2(k,len(program.instructions),instructions,len(program.coefficients),packed,len(program.cars),co,cc,dd,left['multiplicity'],ii,right['multiplicity'],jj,rr,base,output_words,output,shifts,meta)
    if status:raise ArithmeticError(f'Core-coordinate program failed or lacked common precision: {status}; metadata={list(meta)}.')
    infinity=data['place']=='infinity'
    row_states=[i//left['multiplicity'] for i in range(n)];column_states=[j//right['multiplicity'] for j in range(n)]
    kept_rows=[i for i,state in enumerate(row_states) if full_alphabet_gamma is not None or ((state>>2)&3 if infinity else state&3)!=3]
    kept_columns=[j for j,state in enumerate(column_states) if full_alphabet_gamma is not None or (state&3 if infinity else (state>>2)&3)!=3]
    assert len(kept_rows)==len(kept_columns)==(n if full_alphabet_gamma is not None else 3*n//4)
    keep_row_set=set(kept_rows);keep_column_set=set(kept_columns)
    for i in range(n):
        for j in range(n):
            if i not in keep_row_set or j not in keep_column_set:
                index=i*n+j
                if any(output[(2*index)*output_words+w] or output[(2*index+1)*output_words+w] for w in range(output_words)):
                    raise ArithmeticError('A theoretically zero physical-core row or column is nonzero.')
    matrix=[];shift_matrix=[]
    for i in kept_rows:
        row=[];shift_row=[]
        for j in kept_columns:
            index=i*n+j;row.append(tuple(sum(int(output[(2*index+plane)*output_words+w])<<(64*w) for w in range(output_words)) for plane in range(2)));shift_row.append(int(shifts[index]))
        matrix.append(row);shift_matrix.append(shift_row)
    metadata={'place':data['place'],'M':data['M'],'matrix_order':len(matrix),'padded_matrix_order':n,'registers':len(program.instructions),'coefficient_count':len(program.coefficients),'unit_denominator_degree':denominator.bit_length()-1,'unit_denominator_hex':hex(denominator),'program_degree_bound':bound,'maximum_normalized_degree':int(meta[3]),'exact_polynomial_matrix':bound<builder and int(meta[3])<window,'minimum_valuation':int(meta[0]),'known_below':int(meta[1]),'builder_precision':builder,'profile_precision':window,'matrix_seconds':time.monotonic()-start,'row_core_states':[row_states[i] for i in kept_rows],'column_core_states':[column_states[j] for j in kept_columns],'row_charges':[left['coordinate_row_charges'][i] for i in kept_rows],'column_charges':[right['coordinate_column_charges'][j] for j in kept_columns]}
    metadata['row_core_letters']=left['core_letters'];metadata['column_core_letters']=right['core_letters'];metadata['full_core']=bool(data.get('full_core',False))
    metadata['full_alphabet_gamma']=full_alphabet_gamma
    if residual_monomials is not None:metadata['residual_monomials']=residual_monomials
    return matrix,shift_matrix,metadata

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--M',type=int,default=4);parser.add_argument('--place',choices=('infinity','zero'),default='infinity');parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2);parser.add_argument('--quartic-power',type=int,default=0);parser.add_argument('--output',type=Path)
    args=parser.parse_args();data=data_for(args.place,args.M);matrix,shifts,metadata=coordinate_matrix(data,args.alpha,args.beta,[(1,args.quartic_power)])
    start=time.monotonic();result=MatrixOracle(matrix,len(matrix),shifts,precision=metadata['profile_precision']).profile()
    if not metadata['exact_polynomial_matrix'] and max(result['orders'])>=metadata['known_below']:raise ArithmeticError('Profile exceeds known polynomial precision.')
    metadata['profile_seconds']=time.monotonic()-start;metadata['orders']=result['orders'];metadata['alpha']=args.alpha;metadata['beta']=args.beta;metadata['quartic_power']=args.quartic_power
    for key in ('row_core_states','column_core_states','row_charges','column_charges'):metadata.pop(key)
    rendered=json.dumps(metadata,indent=2)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered)

if __name__=='__main__':main()
