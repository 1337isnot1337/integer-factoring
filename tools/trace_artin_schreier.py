"""Formal Artin--Schreier reduction and slicing witnesses for trace systems.

The field model is supplied. Reduction is exact for its binary trace but
does not claim to remove every function relation on a small finite field.
No domain enumeration, character-sum evaluator or factoring counter occurs.
"""


def field_trace(field,width,value):
    answer=0
    for j in range(width):answer ^= field.power(value,1 << j)
    if answer not in (0,1):raise ValueError('Coefficient is outside the trace subfield')
    return answer


def degree(polynomial):
    return max((sum(exponents) for c,exponents in polynomial if c),default=0)


def reduce_polynomial(field,width,variables,polynomial):
    """Return (nonconstant formal AS-reduced terms, constant trace BIT).

    Every nonconstant exponent vector is divided by its common power of
    two; the coefficient receives the inverse Frobenius. Targets stay
    outside the trace, including at even widths.
    """
    if width<1 or field.degree%width or variables<1:
        raise ValueError('Incompatible field, trace width or variable count')
    terms={};constant=0
    for c,exponents in polynomial:
        if not 0<=c<field.size or len(exponents)!=variables or any(e<0 for e in exponents):
            raise ValueError('Invalid ordinary field polynomial')
        if field.power(c,1 << width)!=c:
            raise ValueError('Coefficient is outside the summation subfield')
        if not c:continue
        exponents=tuple(exponents)
        if not any(exponents):
            constant ^= field_trace(field,width,c)
            continue
        twos=min((e & -e).bit_length()-1 for e in exponents if e)
        exponents=tuple(e >> twos for e in exponents)
        c=field.power(c,1 << ((-twos)%width))
        terms[exponents]=terms.get(exponents,0) ^ c
    result=[(c,e) for e,c in sorted(terms.items()) if c]
    assert all(any(e&1 for e in exponents) for c,exponents in result)
    return result,constant


def reduce_system(field,width,variables,predicates):
    """Equivalent echelon constraints in the FORMAL nonconstant AS classes.

    Returns (None, invoice) for a detected contradiction. Other relations
    among finite-field functions may remain; this is not their exact rank.
    """
    rows=[];monomials=set();constant_bits=[]
    for polynomial,bit in predicates:
        if bit not in (0,1):raise ValueError('Binary target required')
        terms,constant=reduce_polynomial(field,width,variables,polynomial)
        rows.append((terms,bit ^ constant))
        constant_bits.append(constant)
        monomials.update(e for c,e in terms)
    monomials=sorted(monomials)
    positions={e:i for i,e in enumerate(monomials)}
    pivots={}
    for i,(terms,target) in enumerate(rows):
        vector=sum(c << (field.degree*positions[e]) for c,e in terms)
        origin=1 << i
        while vector:
            lead=vector.bit_length()-1
            if lead in pivots:
                old,old_target,old_origin=pivots[lead]
                vector ^= old;target ^= old_target;origin ^= old_origin
            else:
                pivots[lead]=(vector,target,origin)
                break
        if not vector and target:
            return None,dict(inconsistent=True,formal_rank=None,constant_bits=constant_bits,
                             contradiction_combination=origin,input_predicates=len(predicates))
    result=[];origins=[]
    mask=field.size-1
    for lead,(vector,target,origin) in sorted(pivots.items(),reverse=True):
        terms=[((vector >> (field.degree*i)) & mask,e) for i,e in enumerate(monomials)]
        result.append(([(c,e) for c,e in terms if c],target))
        origins.append(origin)
    return result,dict(inconsistent=False,formal_rank=len(result),constant_bits=constant_bits,
                       basis_combinations=origins,input_predicates=len(predicates),
                       input_degree=max((degree(f) for f,b in predicates),default=0),
                       reduced_degree=max((degree(f) for f,b in result),default=0))


def slice_witness(field,width,variables,polynomial):
    """A nonzero polynomial witnessing a good positive odd-degree fibre.

    Input must already be nonconstant and formally reduced. A may still
    vanish on every point of a small finite field; the degree/root bound
    remains valid and can then be uninformative.
    """
    reduced,constant=reduce_polynomial(field,width,variables,polynomial)
    assert reduced and reduced==list(polynomial) and constant==0
    axis=next(i for i in range(variables) if any(e[i] for c,e in reduced))
    groups={}
    for c,e in reduced:
        j=e[axis]
        if not j:continue
        twos=(j & -j).bit_length()-1
        odd=j >> twos
        groups.setdefault(odd,[]).append((twos,c,e[:axis]+e[axis+1:]))
    odd=max(groups)
    height=max(k for k,c,e in groups[odd])
    coefficients={}
    for k,c,e in groups[odd]:
        shift=height-k
        out=tuple(x << shift for x in e)
        coefficient=field.power(c,1 << shift)
        coefficients[out]=coefficients.get(out,0) ^ coefficient
    witness=[(c,e) for e,c in sorted(coefficients.items()) if c]
    assert witness
    witness_degree=degree(witness)
    total_degree=degree(reduced)
    assert witness_degree<=total_degree*(1 << height)<=total_degree**2
    return dict(axis=axis,odd_degree=odd,clearing_frobenius=height,
                exceptional_polynomial=witness,exceptional_degree=witness_degree,
                original_degree=total_degree)
