#!/usr/bin/env python3
"""Mixed directions on disjoint fixed-xor physical quartic pairs."""
from check_rectangular_defect_support import pdivmod
from reflection_dual_rail_gaussian import tensor,validate
from reflection_gaussian_physical_frames import norm,add,mul
from reflection_quartic_tiles import gate_data,chi_power
from reflection_paired_quartic_words import chi_order


def validate_mixed_word(M,gates,base_orientation=(0,3)):
    validate(M);m=M//2
    base_orientation=tuple(base_orientation)
    if base_orientation not in ((0,3),(3,0),(1,2),(2,1)):
        raise ValueError('Choose one oriented matched physical pair for the common charge.')
    used=set();flip=None;data=[]
    for gate in gates:
        modes,kinds,powers,pole,numerator,_,_=gate_data(M,gate)
        created=[v for v,k in zip(modes,kinds) if k=='create']
        removed=[v for v,k in zip(modes,kinds) if k=='annihilate']
        if len(created)!=2 or len(removed)!=2:raise ValueError('Use two creations and two annihilations.')
        p={v//m for v in created};n={v//m for v in removed}
        if len(p)!=1 or len(n)!=1:raise ValueError('Each operator kind acts on one physical letter.')
        orientation=(next(iter(p)),next(iter(n)))
        if orientation==base_orientation:direction=1
        elif orientation==base_orientation[::-1]:direction=-1
        else:raise ValueError('Every gate uses the same unordered matched physical pair.')
        indices=tuple(sorted(v%m for v in created))
        if tuple(sorted(v%m for v in removed))!=indices:raise ValueError('Both letters use the same logical pair.')
        current=indices[0]^indices[1]
        if current&(current-1):raise ValueError('Use a single-bit xor pairing.')
        if flip is None:flip=current
        if current!=flip or any(i in used for i in indices):
            raise ValueError('Logical pairs must be disjoint and use one fixed xor bit.')
        used.update(indices);gain=4*flip;bound=gain+(7 if direction==1 else -1)
        if pole>bound:raise ValueError('The pole exceeds this direction\'s proved margin.')
        data.append({'modes':modes,'kinds':kinds,'powers':powers,'pole':pole,
            'numerator':numerator,'indices':indices,'direction':direction,
            'correction_grade':gain+(8 if direction==1 else 0),'pole_bound':bound})
    return {'orientation':base_orientation,'xor_bit':flip,'entries':data}


def construct(M,gates,base_orientation=(0,3)):
    certificate=validate_mixed_word(M,gates,base_orientation)
    values=[norm(a,b,0) for a,b in tensor(M)['values']];raw_pole=0
    for data in certificate['entries']:
        mask=sum(1<<i for i in data['modes'])
        occupied=sum(1<<i for i,kind in zip(data['modes'],data['kinds']) if kind=='annihilate')
        output=[mul((chi_power(data['pole']),0,0),value) for value in values]
        for state,value in enumerate(values):
            if state&mask==occupied:
                output[state^mask]=add(output[state^mask],mul(data['numerator'],value))
        values=output;raw_pole+=data['pole']
    clearing=max(0,-min(shift for a,b,shift in values if a|b))
    polynomial=[(a<<(shift+clearing),b<<(shift+clearing)) for a,b,shift in values]
    content=min(chi_order(plane) for value in polynomial for plane in value if plane)
    cancelled=min(content,raw_pole)
    if cancelled:
        divisor=chi_power(cancelled);output=[]
        for value in polynomial:
            pair=[]
            for plane in value:
                quotient,remainder=pdivmod(plane,divisor);assert not remainder;pair.append(quotient)
            output.append(tuple(pair))
        polynomial=output
    pole=raw_pole-cancelled;assert pole<=M-1
    degree=max(max(value).bit_length()-1 for value in polynomial)
    bound=3*M*M//2+sum(max(0,-min(row['powers']))+max(0,max(row['powers'])-row['pole'])
                       for row in certificate['entries'])
    assert degree-pole<=bound
    q=1<<(M//2)
    return {'values':polynomial,'dimensions':[q]*4,'q_x':q,'q_u':q,'interval':M,
        'degree':degree,'common_w_clearing':clearing,'endpoint_pole_order':pole,
        'raw_endpoint_pole_order':raw_pole,'cancelled_common_chi_order':cancelled,
        'effective_cleared_degree':degree-pole,'effective_degree_bound':bound,
        'strict_endpoint':False,'type_projection_required':True,'gates':gates,
        'base_orientation':certificate['orientation'],'xor_bit':certificate['xor_bit'],
        'directions':[row['direction'] for row in certificate['entries']],
        'per_gate_margins':[row['correction_grade'] for row in certificate['entries']]}
