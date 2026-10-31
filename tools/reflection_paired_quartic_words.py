#!/usr/bin/env python3
"""Growing fixed-xor quartic words with amplified removable endpoint poles."""
from itertools import combinations

from check_rectangular_defect_support import pdivmod
from reflection_dual_rail_gaussian import tensor,validate
from reflection_gaussian_physical_frames import norm,add,mul
from reflection_quartic_tiles import gate_data,chi_power


def validate_word(M,gates):
    validate(M);m=M//2;data=[];used=set();orientation=None;flip=None
    for gate in gates:
        modes,kinds,powers,pole,numerator,_,_=gate_data(M,gate)
        created=[mode for mode,kind in zip(modes,kinds) if kind=='create']
        removed=[mode for mode,kind in zip(modes,kinds) if kind=='annihilate']
        if len(created)!=2 or len(removed)!=2:
            raise ValueError('Each paired quartet has two creations and two annihilations.')
        positive={mode//m for mode in created};negative={mode//m for mode in removed}
        if len(positive)!=1 or len(negative)!=1:raise ValueError('Each operator kind stays on one physical letter.')
        pair=(next(iter(positive)),next(iter(negative)))
        if pair not in ((0,3),(3,0),(1,2),(2,1)):
            raise ValueError('The physical letters must be a matched crossing pair.')
        indices=tuple(sorted(mode%m for mode in created))
        if tuple(sorted(mode%m for mode in removed))!=indices:
            raise ValueError('Creation and annihilation use the same two logical indices.')
        current_flip=indices[0]^indices[1]
        if current_flip&(current_flip-1):raise ValueError('Use a fixed single-bit xor pairing.')
        if orientation is None:orientation=pair;flip=current_flip
        if pair!=orientation or current_flip!=flip:
            raise ValueError('A word uses one oriented physical pair and one fixed xor bit.')
        if any(index in used for index in indices):raise ValueError('The logical index pairs must be disjoint.')
        used.update(indices)
        gain=4*flip
        if pole>7+gain:raise ValueError('The pole exceeds the proved Pfaffian row-factor budget.')
        data.append({'modes':modes,'kinds':kinds,'powers':powers,'pole':pole,
                     'numerator':numerator,'indices':indices,'pfaffian_gain':gain})
    return {'orientation':orientation,'xor_bit':flip,'entries':data,
            'pole_bound':0 if flip is None else 7+4*flip}


def chi_order(value):
    assert value
    order=0
    while True:
        quotient,remainder=pdivmod(value,3)
        if remainder:return order
        value=quotient;order+=1


def construct(M,gates):
    """Literal M4/M8 coefficients; the validated family is uniform in M."""
    certificate=validate_word(M,gates)
    values=[norm(a,b,0) for a,b in tensor(M)['values']]
    raw_pole=0
    for data in certificate['entries']:
        mask=sum(1<<i for i in data['modes'])
        occupied=sum(1<<i for i,kind in zip(data['modes'],data['kinds']) if kind=='annihilate')
        denominator=(chi_power(data['pole']),0,0)
        output=[mul(denominator,value) for value in values]
        for state,value in enumerate(values):
            if state&mask==occupied:
                output[state^mask]=add(output[state^mask],mul(data['numerator'],value))
        values=output;raw_pole+=data['pole']
    clearing=max(0,-min(shift for a,b,shift in values if a|b))
    polynomial=[(a<<(shift+clearing),b<<(shift+clearing)) for a,b,shift in values]
    content=min(chi_order(plane) for value in polynomial for plane in value if plane)
    cancelled=min(content,raw_pole)
    if cancelled:
        divisor=chi_power(cancelled);reduced=[]
        for value in polynomial:
            pair=[]
            for plane in value:
                quotient,remainder=pdivmod(plane,divisor);assert not remainder;pair.append(quotient)
            reduced.append(tuple(pair))
        polynomial=reduced
    pole=raw_pole-cancelled
    assert pole<=M-1
    degree=max(max(value).bit_length()-1 for value in polynomial)
    effective_bound=3*M*M//2+sum(max(0,-min(row['powers']))+
        max(0,max(row['powers'])-row['pole']) for row in certificate['entries'])
    assert degree-pole<=effective_bound
    q=1<<(M//2)
    return {'values':polynomial,'dimensions':[q]*4,'q_x':q,'q_u':q,'interval':M,
        'degree':degree,'common_w_clearing':clearing,'endpoint_pole_order':pole,
        'raw_endpoint_pole_order':raw_pole,'cancelled_common_chi_order':cancelled,
        'effective_cleared_degree':degree-pole,'effective_degree_bound':effective_bound,
        'strict_endpoint':False,'type_projection_required':True,'gates':gates,
        'oriented_physical_pair':certificate['orientation'],'xor_bit':certificate['xor_bit'],
        'per_gate_pole_bound':certificate['pole_bound']}


def paired_gates(M,xor_bit,*,orientation=(0,3),pole=None,power=0,gamma=1,indices=None):
    """A complete or selected subword of the fixed xor pairing."""
    validate(M);m=M//2
    if xor_bit<1 or xor_bit>=m or xor_bit&(xor_bit-1):raise ValueError('Use a logical mode bit.')
    if indices is None:indices=[i for i in range(m) if not i&xor_bit]
    if pole is None:pole=7+4*xor_bit
    positive,negative=orientation
    return [{'modes':[positive*m+i,positive*m+(i^xor_bit),negative*m+i,negative*m+(i^xor_bit)],
             'kinds':['create','create','annihilate','annihilate'],
             'powers':[power],'pole':pole,'gamma':gamma} for i in indices]
