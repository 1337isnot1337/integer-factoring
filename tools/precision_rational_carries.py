"""Collected total unit inverses and packed systems for Teichmuller sums.

The polynomial evaluator, ring model and parity polynomial are credited to
teichmuller_cyclic_carries.py / PR71 b9cd96f. This module compiles additional
inputs to that evaluator; it never enumerates the field domain.
"""
from teichmuller_cyclic_carries import polynomial_sum, parity_lift_coefficients


def factor_degree(factor):
    return max((sum(counts) for coefficient,counts in factor[1] if coefficient),default=0)


def conjugate_factor(ring,width,row_modulus,variables,factor,shift):
    """Apply sigma^shift on values and wrap Teichmuller binary positions."""
    row,terms=factor
    if width<1 or row_modulus<1 or width%row_modulus or not 0<=row<row_modulus or shift<0:
        raise ValueError('Invalid synchronized row')
    columns=width//row_modulus
    new_row=(row+shift)%row_modulus
    rotation=(row+shift)//row_modulus
    result=[]
    for coefficient,counts in terms:
        if len(counts)!=variables*columns or any(x<0 for x in counts):
            raise ValueError('Invalid explicit count vector')
        changed=[0]*len(counts)
        for i in range(variables):
            for j in range(columns):
                changed[i*columns+(j+rotation)%columns]=counts[i*columns+j]
        result.append((ring.frobenius(coefficient,shift),tuple(changed)))
    return new_row,result


def rational_sum(ring,width,row_modulus,variables,numerators,denominators,nonzero=False):
    """Sum numerator products times total unit inverses of denominators.

    A nonunit denominator makes the value zero. This is not division by
    a nonzero zero divisor. Return (ring value, explicit resource invoice).
    """
    numerator_degree=sum(map(factor_degree,numerators))
    denominator_degree=sum(map(factor_degree,denominators))
    if not denominators:
        answer=polynomial_sum(ring,width,row_modulus,variables,numerators,nonzero)
        return answer,dict(queries=1,max_degree_budget=numerator_degree,
                           summation_width=width,coefficient_degree=ring.degree,
                           relative_columns=width//row_modulus)
    norms=[];adjugates=[]
    for factor in denominators:
        for shift in range(ring.degree):
            changed=conjugate_factor(ring,width,row_modulus,variables,factor,shift)
            norms.append(changed)
            if shift: adjugates.append(changed)
    answer=0;queries=0
    # Formal division of E_(L+1)(X) by X, before evaluating any ring value.
    inverse_coefficients=parity_lift_coefficients(ring.precision+1)[1:]
    for exponent,coefficient in enumerate(inverse_coefficients):
        if coefficient:
            factors=list(numerators)+adjugates+norms*exponent
            value=polynomial_sum(ring,width,row_modulus,variables,factors,nonzero)
            answer=ring.add(answer,ring.mul(ring.integer(coefficient),value))
            queries+=1
    budget=numerator_degree+((2*ring.precision+1)*ring.degree-1)*denominator_degree
    return answer,dict(queries=queries,max_degree_budget=budget,
                       summation_width=width,coefficient_degree=ring.degree,
                       relative_columns=width//row_modulus,
                       denominator_factors=len(denominators),
                       norm_factors=len(norms),adjugate_factors=len(adjugates))


def packed_system_count(ring,width,variables,equations,nonzero=False):
    """Count a system over the width-bit residue subfield modulo 2^L.

    The supplied basic irreducible quotient has degree >= width*len(equations)
    divisible by width. Coefficients must lie in the fixed subring of
    sigma^width. Powers of its generator give an independent subfield basis.
    Each equation uses ordinary nonnegative exponents, with no free inverse
    or arbitrary arithmetic-circuit input. No summation variables are added.
    """
    if width<1 or ring.degree%width or len(equations)>ring.degree//width:
        raise ValueError('Insufficient unramified coefficient extension')
    packed=[];degree=0
    for i,equation in enumerate(equations):
        basis=ring.power(ring.generator,i)
        for coefficient,counts in equation:
            if len(counts)!=variables or any(x<0 for x in counts):
                raise ValueError('Invalid ordinary polynomial monomial')
            if ring.frobenius(coefficient,width)!=coefficient:
                raise ValueError('Equation coefficient is outside the subfield ring')
            packed.append((ring.mul(basis,coefficient),tuple(counts)))
            if coefficient: degree=max(degree,sum(counts))
    constant=polynomial_sum(ring,width,width,variables,[],nonzero)
    if not equations:
        return constant,dict(queries=1,variables=variables,precision=ring.precision,
                              coefficient_degree=ring.degree,max_degree_budget=0)
    norms=[conjugate_factor(ring,width,width,variables,(0,packed),shift)
           for shift in range(ring.degree)]
    answer=constant;queries=1
    for exponent,coefficient in enumerate(parity_lift_coefficients(ring.precision)):
        if coefficient:
            moment=polynomial_sum(ring,width,width,variables,norms*exponent,nonzero)
            answer=ring.sub(answer,ring.mul(ring.integer(coefficient),moment))
            queries+=1
    assert not any(ring.coefficients(answer)[1:])
    return answer,dict(queries=queries,variables=variables,precision=ring.precision,
                       coefficient_degree=ring.degree,summation_width=width,
                       equations=len(equations),relative_columns=1,
                       max_degree_budget=(2*ring.precision-1)*ring.degree*degree)
