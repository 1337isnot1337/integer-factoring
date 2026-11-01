#!/usr/bin/env python3
"""Complete degree-six cube coefficient fibers in exact generic row frames.

The 256 endpoint equations leave 6912 field coefficients. There is no
physical parity, translation, or axis symmetry restriction.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import cube_f4_bivariate_smith as exact
from cube_homogeneous_gaussian import axis_parts

ANCHOR=(2,2,2)

def serialize_tensor(tensor):
    return [[[a,b,c,value] for (a,b,c),value in sorted(polynomial.items()) if value] for polynomial in tensor]

def parse_tensor(record):
    rows=record.get('tensor',record.get('entries')) if isinstance(record,dict) else record
    if len(rows)!=256:raise ValueError('Expected all 256 physical coefficients.')
    result=[]
    for word,row in enumerate(rows):
        if isinstance(row,dict):
            if row.get('word',word)!=word:raise ValueError('Physical words are not in canonical order.')
            row=row.get('coefficients',row.get('monomials',row.get('terms')))
        polynomial={}
        for term in row:
            if len(term) not in (3,4):raise ValueError('Expected three exponents and an optional field coefficient.')
            exponent=tuple(term[:3]);coefficient=term[3] if len(term)==4 else 1
            polynomial[exponent]=polynomial.get(exponent,0)^coefficient
        result.append({e:c for e,c in polynomial.items() if c})
    exact.validate_tensor(result)
    return tuple(result)

def tensor_digest(tensor):
    return hashlib.sha256(json.dumps(serialize_tensor(tensor),separators=(',',':')).encode()).hexdigest()

def coefficient_basis(field):
    if field not in (2,4):raise ValueError('Supported coefficient fields have order 2 or 4.')
    directions=[];symbolic=[{} for _ in range(256)];components=(1,) if field==2 else (1,2)
    for mask in range(256):
        for x in range(7):
            for y in range(7-x):
                exponent=(x,y,6-x-y)
                if exponent==ANCHOR:continue
                bits=[0,0]
                for component,coefficient in enumerate(components):
                    bits[component]=1<<len(directions);directions.append((mask,exponent,coefficient))
                symbolic[mask][exponent]=tuple(bits)
                a,b=symbolic[mask].get(ANCHOR,(0,0));symbolic[mask][ANCHOR]=a^bits[0],b^bits[1]
    assert len(directions)==6912*(1 if field==2 else 2)
    return directions,symbolic

def tensor_from_assignment(seed,directions,value):
    if value<0 or value>>len(directions):raise ValueError('Coefficient assignment is outside the declared space.')
    result=[dict(polynomial) for polynomial in seed]
    while value:
        bit=value&-value;value^=bit;mask,exponent,coefficient=directions[bit.bit_length()-1]
        for term in (exponent,ANCHOR):result[mask][term]=result[mask].get(term,0)^coefficient
    return tuple({e:c for e,c in polynomial.items() if c} for polynomial in result)

def profiles(tensor,verify=True):
    result=[]
    for axis in range(3):
        record=exact.profile_tensor(tensor,axis)
        if verify:exact.verify_matrix(exact.tensor_matrix(tensor,axis),record)
        result.append(record)
    return result

def insert_equation(pivots,row,mask):
    while row&mask:
        pivot=(row&mask).bit_length()-1
        if pivot not in pivots:pivots[pivot]=row;return True
        row^=pivots[pivot]
    return row==0

def solve(system,targets):
    if len(targets)!=48 or any(not isinstance(v,int) or not 0<=v<=cap for v,cap in zip(targets,system['caps'])):
        raise ValueError('Invalid row-floor target.')
    pivots={};mask=(1<<system['variables'])-1
    for row,target in enumerate(targets):
        for power in range(target):
            for equation in system['equations'][row][power]:
                if not insert_equation(pivots,equation,mask):return None
    return pivots

def complete_assignment(pivots,free_value,variables):
    if free_value<0 or free_value>>variables or any(free_value>>p&1 for p in pivots):
        raise ValueError('The assignment does not use only free coordinates.')
    value=free_value
    for pivot,row in sorted(pivots.items()):value|=(((row>>variables)&1)^((row&value).bit_count()&1))<<pivot
    return value

def build_system(seed,field=2,extra_order=6):
    exact.validate_tensor(seed)
    if field not in (2,4) or any(c>=field for polynomial in seed for c in polynomial.values()):
        raise ValueError('Seed coefficients are outside the declared field.')
    if not isinstance(extra_order,int) or extra_order<1:raise ValueError('The extra-order window must be positive.')
    directions,symbolic=coefficient_basis(field);variables=len(directions);frames=profiles(seed)
    reference=[order for record in frames for order in record['orders']]
    caps=[order+extra_order for order in reference];table={};expanded_terms=0
    for axis in range(3):
        first_vertices,second_vertices=axis_parts(axis);other=next(i for i in range(3) if i!=axis)
        frame=frames[axis]['frames']
        for word in range(256):
            row_word=sum(((word>>position)&1)<<i for i,position in enumerate(first_vertices))
            column_word=sum(((word>>position)&1)<<i for i,position in enumerate(second_vertices))
            polynomial=dict(symbolic[word])
            for exponent,coefficient in seed[word].items():
                a,b=polynomial.get(exponent,(0,0))
                if coefficient&1:a^=1<<variables
                if coefficient&2:b^=1<<variables
                polynomial[exponent]=a,b
            for i in range(16):
                row=16*axis+i
                for fx,fy,coefficient in exact.terms(frame[i][row_word]):
                    for exponent,(a,b) in polynomial.items():
                        x=fx+exponent[axis]
                        if x>=caps[row]:continue
                        real=(a if coefficient&1 else 0)^(b if coefficient&2 else 0)
                        alpha=(a if coefficient&2 else 0)^(b if bool(coefficient&1)!=bool(coefficient&2) else 0)
                        for component,value in enumerate((real,alpha)):
                            if not value:continue
                            key=row,x,column_word,fy+exponent[other],component
                            table[key]=table.get(key,0)^value;expanded_terms+=1
    equations=[[set() for _ in range(cap)] for cap in caps]
    for (row,x,column,y,component),value in table.items():
        if value:equations[row][x].add(value)
    return dict(field=field,variables=variables,directions=directions,frames=frames,reference=reference,caps=caps,
                expanded_terms=expanded_terms,equations=[[tuple(sorted(values)) for values in row] for row in equations],
                seed_digest=tensor_digest(seed),extra_order=extra_order)

def verify_row_floors(tensor,system,targets):
    """Independent polynomial products against the fixed reference frames."""
    if len(targets)!=48 or any(not isinstance(t,int) or t<0 for t in targets):
        raise ValueError('Expected 48 nonnegative row-floor targets.')
    exact.validate_tensor(tensor);entries=0;actual=[]
    for axis in range(3):
        matrix=exact.tensor_matrix(tensor,axis);frame=system['frames'][axis]['frames']
        for row in range(16):
            minimum=None
            for column in range(16):
                value=exact.ZERO
                for k in range(16):value=exact.add(value,exact.multiply(frame[row][k],matrix[k][column]))
                powers=[x for x,y,c in exact.terms(value)]
                if powers:
                    at=min(powers);minimum=at if minimum is None else min(minimum,at)
                    assert at>=targets[16*axis+row]
                entries+=1
            assert minimum is not None;actual.append(minimum)
    return dict(polynomial_product_entries=entries,row_orders=actual)

def pressure(orders):
    """Numerical common-degree-six pressure; never a proof certificate."""
    if len(orders)!=3 or any(len(p)!=16 for p in orders):
        raise ValueError('Expected three full profiles.')
    def evaluate(u):
        log_sum=weighted=0.0
        for profile in orders:
            minimum=min(profile);weights=[math.exp(-u*(s-minimum)) for s in profile]
            total=sum(weights)
            log_sum+=math.log(total)-u*minimum
            weighted+=sum(w*s for w,s in zip(weights,profile))/total
        return (log_sum/3+2*u)/(4*math.log(2)),weighted/3
    if evaluate(0)[1]<=2:return dict(rho=1.0,u=0.0,base=1.0,numerical=True)
    if sum(min(p) for p in orders)>=6:
        minimum_sum=sum(min(p) for p in orders)
        limit=sum(math.log(p.count(min(p))) for p in orders)/(12*math.log(2))
        return dict(rho=limit if minimum_sum==6 else -math.inf,u=None,base=0.0,numerical=True)
    left,right=0.0,1.0
    while evaluate(right)[1]>2:right*=2
    for _ in range(64):
        middle=(left+right)/2
        if evaluate(middle)[1]>2:left=middle
        else:right=middle
    u=(left+right)/2
    return dict(rho=evaluate(u)[0],u=u,base=math.exp(-u),numerical=True)

def certify_pressure(orders,numerator,denominator,base_numerator,base_denominator):
    """Exact rational-base upper certificate for a declared full profile."""
    if (len(orders)!=3 or any(len(p)!=16 for p in orders)
            or any(not isinstance(s,int) or s<0 for p in orders for s in p)
            or not 0<numerator<=denominator or not 0<base_numerator<base_denominator):
        raise ValueError('Invalid profile or rational pressure parameters.')
    a,b=base_numerator,base_denominator;product=1;degree=0
    for profile in orders:
        maximum=max(profile);degree+=maximum
        product*=sum(a**s*b**(maximum-s) for s in profile)
    left=(product*b**6)**denominator
    right=2**(12*numerator)*(b**degree*a**6)**denominator
    delta=abs(right-left)
    return dict(passed=left<right,rho=[numerator,denominator],base=[a,b],
                left_bits=left.bit_length(),right_bits=right.bit_length(),
                difference_sha256=hashlib.sha256(delta.to_bytes((delta.bit_length()+7)//8,'big')).hexdigest(),
                upper_bound_below_cube_threshold=3*numerator*numerator<denominator*denominator)

def source_hashes():
    names=('cube_flat_coefficient_fibers.py','cube_f4_bivariate_smith.py',
           'cube_f4_bivariate_smith.cpp','cube_homogeneous_gaussian.py')
    return {name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
