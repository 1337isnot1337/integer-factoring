#!/usr/bin/env python3
"""Uniform reflection-pair compression of a strict Pascal Gaussian.

The formula is a Laurent Gaussian times one common monomial. Complete
coefficients are computed in a regular zero gauge and then restored;
restoration rejects any unexpected pole. For interval M, each physical
alphabet has 2^(M/2) letters, instead of 2^M.
"""
from __future__ import annotations
from functools import lru_cache
from check_schur_gaussian_reduction import polynomial_multiply


def validate(interval):
    if not isinstance(interval,int) or interval<4 or interval&(interval-1):
        raise ValueError('Use a power-of-two interval of length at least four.')


def reflection_letters(interval):
    validate(interval);m=interval//2
    return tuple(sum(1<<(interval-1-i if state>>i&1 else i) for i in range(m))
                 for state in range(1<<m))


def logical_heights(interval,place):
    validate(interval);m=interval//2
    if place not in ('zero','infinity'):raise ValueError('Choose zero or infinity.')
    active=(0,1) if place=='zero' else (2,3)
    return tuple(2*(interval-1-2*i) if leg in active else 0
                 for leg in range(4) for i in range(m))


def laurent_kernel(interval):
    """Symmetric alternating kernel in x,y,u,v variable order."""
    validate(interval);m=interval//2;size=4*m
    kernel=[[{} for _ in range(size)] for _ in range(size)]
    def put(leg1,i,leg2,j,exponents):
        a,b=leg1*m+i,leg2*m+j
        if a==b:raise ArithmeticError('The Gaussian kernel must be alternating.')
        terms={}
        for exponent in exponents:
            if exponent in terms:del terms[exponent]
            else:terms[exponent]=1
        if kernel[a][b]:raise ArithmeticError('A kernel pair was specified twice.')
        kernel[a][b]=terms;kernel[b][a]=terms.copy()
    for leg in range(4):
        second=2*interval if leg in (2,3) else -2*interval
        for i in range(m//2):put(leg,i,leg,m-1-i,(0,second))
    for i in range(m):
        for j in range(m):
            exponent=4*(j&~i)
            put(3,i,1,j,(exponent,exponent+2))
            put(2,i,0,j,(exponent,exponent+2))
            exponent=4*(m-1-(i|j))
            put(3,i,2,j,(exponent,exponent+2))
            exponent=-2*interval+4*(i&j)
            put(1,i,0,j,(exponent,exponent+2))
            if i==j:
                put(3,i,0,j,(0,));put(1,i,2,j,(0,))
    return kernel


def gauged_kernel(interval,place):
    kernel=laurent_kernel(interval);heights=logical_heights(interval,place)
    sign=1 if place=='zero' else -1
    output=[]
    for i,row in enumerate(kernel):
        values=[]
        for j,entry in enumerate(row):
            exponents=[sign*power+heights[i]+heights[j] for power in entry]
            if any(power<0 for power in exponents):
                raise ArithmeticError('The natural gauged kernel has an unexpected pole.')
            values.append(sum(1<<power for power in exponents))
        output.append(values)
    return output


def tensor(interval):
    """Construct every exact polynomial coefficient for M≤8."""
    validate(interval)
    if interval>8:raise ValueError('Use the implicit kernel and profiles above M=8.')
    kernel=gauged_kernel(interval,'zero');heights=logical_heights(interval,'zero')
    size=len(kernel);scalar=interval*interval//2
    @lru_cache(None)
    def pfaffian(mask):
        if not mask:return 1
        if mask.bit_count()%2:return 0
        first=mask&-mask;i=first.bit_length()-1;rest=mask^first;value=0;remaining=rest
        while remaining:
            chosen=remaining&-remaining;remaining^=chosen;j=chosen.bit_length()-1
            if kernel[i][j]:value^=polynomial_multiply(kernel[i][j],pfaffian(rest^chosen))
        return value
    values=[]
    for mask in range(1<<size):
        value=pfaffian(mask)<<scalar
        shift=sum(height for i,height in enumerate(heights) if mask>>i&1)
        if value&((1<<shift)-1):raise ArithmeticError('The restored compressed Gaussian has a pole.')
        values.append((value>>shift,0))
    q=1<<(interval//2)
    for mask,(value,_) in enumerate(values):
        x,y,u,v=((mask>>(interval//2*leg))&(q-1) for leg in range(4))
        assert value.bit_count()%2==int(x==v and y==u)
    return {'values':values,'dimensions':[q]*4,'q_x':q,'q_u':q,'denominator_degree':0,
            'strict_endpoint':True,'interval':interval,'physical_modes':interval//2,
            'degree':max(value.bit_length()-1 for value,_ in values),
            'scalar_w_power':scalar,'x_letters':reflection_letters(interval),
            'y_letters':reflection_letters(interval)}


def profile_factors(interval):
    validate(interval)
    return tuple(range(2,4*interval,4))


def physical_bit_energies(interval):
    validate(interval)
    return tuple(4*(interval-1-2*i) for i in range(interval//2))
