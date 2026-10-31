#!/usr/bin/env python3
"""Literal strict and type-projected quartic words on reflection alphabets.

The direct coefficient constructor uses the existing M<=8 control tensor.
The companion note proves the family and degree bounds for every interval.
"""
from itertools import product

from check_rectangular_defect_support import pmul
from reflection_dual_rail_gaussian import tensor
from reflection_gaussian_physical_frames import norm, add, mul
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import flatten


def chi_power(power):
    """The binary polynomial (1+w)^power, by the exact submask identity."""
    assert isinstance(power,int) and power>=0
    result=0
    sub=power
    while True:
        result|=1<<sub
        if sub==0:
            return result
        sub=(sub-1)&power


def charge_vector(M):
    m=M//2
    return tuple((1 if leg<2 else -1)*i for leg in range(4) for i in range(m))


def gate_data(M,gate):
    modes=tuple(gate['modes'])
    kinds=tuple(gate.get('kinds',('create',)*4))
    powers=tuple(gate['powers'])
    gamma=gate.get('gamma',1)
    pole=gate.get('pole',0)
    if (len(modes)!=4 or len(set(modes))!=4 or any(not 0<=i<2*M for i in modes)
            or len(kinds)!=4 or any(kind not in ('create','annihilate') for kind in kinds)):
        raise ValueError('Use four distinct physical modes and four operator kinds.')
    if (not powers or len(set(powers))!=len(powers) or gamma not in (1,2,3)
            or not isinstance(pole,int) or pole<0):
        raise ValueError('Use distinct numerator powers, a nonzero F4 color, and a nonnegative pole order.')
    low=min(powers)
    polynomial=sum(1<<(power-low) for power in powers)
    numerator=(polynomial if gamma&1 else 0,polynomial if gamma&2 else 0,low)
    charge=charge_vector(M)
    increase=sum(charge[i]*(1 if kind=='create' else -1) for i,kind in zip(modes,kinds))
    endpoint=gamma if len(powers)%2 else 0
    return modes,kinds,powers,pole,numerator,increase,endpoint


def construct(M,gates):
    """Return a common polynomial numerator and its fully charged scalar data."""
    data=[gate_data(M,gate) for gate in gates]
    has_pole=any(row[3] for row in data)
    strict=not has_pole and all(row[6]==0 for row in data)
    loss=M//2-1
    if has_pole:
        if any(increase<=loss*pole for _,_,_,pole,_,increase,_ in data):
            raise ValueError('Every gate in a removable-pole word must have charge increase >(m-1)*pole.')
    elif not strict:
        if any(endpoint and increase<=0 for _,_,_,_,_,increase,endpoint in data):
            raise ValueError('Every nonvanishing endpoint gate must strictly increase paired charge.')
    values=[norm(low,high,0) for low,high in tensor(M)['values']]
    pole_total=0
    for modes,kinds,powers,pole,numerator,increase,endpoint in data:
        mask=sum(1<<i for i in modes)
        occupied=sum(1<<i for i,kind in zip(modes,kinds) if kind=='annihilate')
        denominator=(chi_power(pole),0,0)
        output=[mul(denominator,value) for value in values]
        for state in range(len(values)):
            if state&mask==occupied:
                output[state^mask]=add(output[state^mask],mul(numerator,values[state]))
        values=output
        pole_total+=pole
    clearing=max(0,-min(shift for low,high,shift in values if low|high))
    polynomial=[(low<<(shift+clearing),high<<(shift+clearing))
                for low,high,shift in values]
    degree=max(max(value).bit_length()-1 for value in polynomial)
    clearing_bound=sum(max(0,-min(row[2])) for row in data)
    effective_bound=3*M*M//2+clearing_bound+sum(max(0,max(row[2])-row[3]) for row in data)
    assert clearing<=clearing_bound and degree-pole_total<=effective_bound
    m=M//2;q=1<<m;charge=charge_vector(M)
    if not has_pole:
        for state,(low,high) in enumerate(polynomial):
            x,y,u,v=((state>>(leg*m))&(q-1) for leg in range(4))
            value=(low.bit_count()&1)+2*(high.bit_count()&1)
            energy=sum(charge[i] for i in range(2*M) if state>>i&1)
            if strict or energy<=0:
                assert value==int(x==v and y==u)
    return {'values':polynomial,'dimensions':[q]*4,'q_x':q,'q_u':q,'interval':M,
            'degree':degree,'common_w_clearing':clearing,'endpoint_pole_order':pole_total,
            'effective_cleared_degree':degree-pole_total,
            'effective_degree_bound':effective_bound,'clearing_bound':clearing_bound,
            'strict_endpoint':strict,'type_projection_required':not strict,
            'charge':charge,'gates':gates}


def complete_profiles(seed,heights):
    """Original rational Fock orders at zero and infinity, with every scalar restored."""
    q=seed['q_x']
    if len(heights)!=2 or any(len(side)!=4 or any(len(leg)!=q for leg in side) for side in heights):
        raise ValueError('Supply four full physical-symbol height lists at each of two places.')
    if seed['endpoint_pole_order'] and 'full_rank_witness' not in seed:
        seed['full_rank_witness']=full_rank_certificate(seed)
    profiles=[]
    for side in (0,1):
        matrix,rows,columns=flatten(seed,side)
        rowlegs=(0,1) if side==0 else (0,2)
        collegs=(2,3) if side==0 else (1,3)
        rw=[sum(heights[side][leg][letter] for leg,letter in zip(rowlegs,pair)) for pair in rows]
        cw=[sum(heights[side][leg][letter] for leg,letter in zip(collegs,pair)) for pair in columns]
        precision=max(512,1<<seed['degree'].bit_length())
        while True:
            oracle=MatrixOracle(matrix,len(matrix),precision=precision)
            try:
                record=oracle.profile(rw,cw,infinity=bool(side))
                break
            except ArithmeticError:
                if precision>=8192:
                    raise
                precision*=2
        scalar=(-seed['common_w_clearing'] if side==0 else
                seed['common_w_clearing']+seed['endpoint_pole_order'])
        assert all(value==round(value) for value in record['orders'])
        profiles.append([round(value)+scalar for value in record['orders']])
    return profiles


def field_multiply(a,b):
    """F16 in the basis defined by x^4+x+1."""
    result=0
    while b:
        if b&1:
            result^=a
        b>>=1
        a<<=1
        if a&16:
            a^=19
    return result


def field_power(a,power):
    result=1
    while power:
        if power&1:
            result=field_multiply(result,a)
        a=field_multiply(a,a)
        power>>=1
    return result


def evaluate_at(value,w):
    omega=field_power(2,5)
    low,high=value
    result=0;power=1
    while low|high:
        coefficient=(low&1)^(omega if high&1 else 0)
        result^=field_multiply(coefficient,power)
        power=field_multiply(power,w)
        low>>=1;high>>=1
    return result


def finite_field_rank(matrix):
    work=[list(row) for row in matrix]
    height=len(work);width=len(work[0]);rank=0
    for column in range(width):
        pivot=next((row for row in range(rank,height) if work[row][column]),None)
        if pivot is None:
            continue
        work[rank],work[pivot]=work[pivot],work[rank]
        inverse=field_power(work[rank][column],14)
        work[rank]=[field_multiply(value,inverse) for value in work[rank]]
        for row in range(rank+1,height):
            if work[row][column]:
                value=work[row][column]
                work[row]=[a^field_multiply(value,b) for a,b in zip(work[row],work[rank])]
        rank+=1
        if rank==height:
            break
    return rank


def full_rank_certificate(seed):
    """An independent finite specialization control of the full-rank theorem.

    All-length paired-type endpoint identities already prove full local rank.
    Failure to find a witness in F16 is an unresolved specialization control,
    not a rank-deficiency theorem.
    """
    result=[]
    for side in (0,1):
        matrix,_,_=flatten(seed,side)
        for value in range(2,16):
            rank=finite_field_rank([[evaluate_at(entry,value) for entry in row] for row in matrix])
            if rank==len(matrix):
                result.append({'side':side,'field_value':value,'rank':rank})
                break
        else:
            raise ArithmeticError('No full-rank F16 specialization was found.')
    return result
