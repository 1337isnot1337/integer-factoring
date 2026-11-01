"""Growing-precision trace predicates in a divided-power coefficient algebra.

Credits: square_zero_trace_contraction.py and trace_constraint_contraction.py
for the packet/algebra construction; teichmuller_cyclic_carries.py and
precision_rational_carries.py for the ring carry and collected scalar inverse
interfaces. No inverse of a general divided-power element is exposed.
"""
from math import comb

from precision_rational_carries import rational_sum, factor_degree
from teichmuller_cyclic_carries import parity_lift_coefficients
from trace_constraint_contraction import independent_constraints


def residue_bits(ring,value):
    return sum((c & 1) << i for i,c in enumerate(ring.coefficients(value)))


class DividedPowerRing:
    """Free truncated divided powers over the supplied unramified base ring.

    degree is the BASE coefficient-field degree used by scalar moments.
    dimension is the separately charged algebra dimension, not a field size.
    Elements pack base-ring coefficients in fixed-width integer blocks.
    """
    def __init__(self,base,rank,limit):
        if rank<0 or limit<0:raise ValueError('Nonnegative algebra dimensions required')
        self.base=base;self.rank=rank;self.limit=limit
        self.precision=base.precision;self.degree=base.degree
        self.modulus=base.modulus;self.radix=limit+1
        self.dimension=self.radix**rank
        self.coefficient_bits=base.precision*base.degree
        self.mask=(1 << self.coefficient_bits)-1
        self.places=[self.radix**i for i in range(rank)]
        self.product_cache={};self.exponent_cache={}

    def terms(self,value):
        if value<0 or value.bit_length()>self.dimension*self.coefficient_bits:
            raise ValueError('Noncanonical coefficient-algebra element')
        while value:
            index=((value & -value).bit_length()-1)//self.coefficient_bits
            coefficient=(value >> (index*self.coefficient_bits)) & self.mask
            yield index,coefficient
            value^=coefficient << (index*self.coefficient_bits)

    def exponents(self,index):
        if index not in self.exponent_cache:
            self.exponent_cache[index]=tuple((index//place)%self.radix for place in self.places)
        return self.exponent_cache[index]

    def tag(self,index,coefficient=1):
        if not 0<=index<self.dimension or not 0<=coefficient<=self.mask:
            raise ValueError('Invalid divided-power basis coefficient')
        return coefficient << (index*self.coefficient_bits)

    def integer(self,value):return self.base.integer(value)

    def add(self,a,b):
        answer=a
        for index,c in self.terms(b):
            shift=index*self.coefficient_bits
            old=(answer >> shift) & self.mask
            new=self.base.add(old,c)
            answer^=(old^new) << shift
        return answer

    def neg(self,a):
        return sum(self.tag(index,self.base.neg(c)) for index,c in self.terms(a))

    def sub(self,a,b):return self.add(a,self.neg(b))

    def product_data(self,i,j):
        key=(min(i,j),max(i,j))
        if key not in self.product_cache:
            weight=1
            for a,b in zip(self.exponents(i),self.exponents(j)):
                if a+b>self.limit:
                    weight=0;break
                weight=weight*comb(a+b,a)%self.modulus
            self.product_cache[key]=(i+j,weight) if weight else None
        return self.product_cache[key]

    def mul(self,a,b):
        answer=0;left=list(self.terms(a));right=list(self.terms(b))
        for i,x in left:
            for j,y in right:
                data=self.product_data(i,j)
                if data is None:continue
                index,weight=data
                value=self.base.mul(self.base.mul(x,y),self.base.integer(weight))
                answer=self.add(answer,self.tag(index,value))
        return answer

    def power(self,a,exponent):
        if exponent<0:raise ValueError('No divided-power inverse is defined')
        answer=1
        while exponent:
            if exponent & 1:answer=self.mul(answer,a)
            a=self.mul(a,a);exponent >>=1
        return answer

    def frobenius(self,a,power=1):
        # This automorphism fixes the formal gamma basis. It is not a^2.
        return sum(self.tag(index,self.base.frobenius(c,power)) for index,c in self.terms(a))

    def mod_two(self,a):
        if a >> self.coefficient_bits:
            raise ValueError('Only embedded scalar denominators are supported')
        return self.base.mod_two(a)

    def unit_inverse(self,a):
        if a >> self.coefficient_bits:
            raise ValueError('No inverse of a general divided-power element is exposed')
        return self.base.unit_inverse(a)

    def indicator_functional(self,a):
        coeffs=[-c for c in parity_lift_coefficients(self.precision)]
        coeffs[0]+=1
        answer=0
        for index,c in self.terms(a):
            weight=1
            for e in self.exponents(index):weight=weight*coeffs[e]%self.modulus
            if weight:answer=self.base.add(answer,self.base.mul(c,self.base.integer(weight)))
        return answer


def polynomial_product(ring,left,right):
    answer={}
    for a,x in left.items():
        for b,y in right.items():
            counts=tuple(i+j for i,j in zip(a,b))
            value=ring.add(answer.get(counts,0),ring.mul(x,y))
            if value:answer[counts]=value
            else:answer.pop(counts,None)
    return answer


def exponential_factor(algebra,axis,polynomial,dimension):
    """Acquire sum_k gamma_(k e_axis) Q^k as an explicit row polynomial."""
    zero=(0,)*dimension;power={zero:1};answer={zero:1}
    acquisition_terms=1
    for k in range(1,algebra.limit+1):
        power=polynomial_product(algebra.base,power,polynomial)
        acquisition_terms+=len(power)
        if not power:break
        for counts,c in power.items():
            value=algebra.tag(k*algebra.places[axis],c)
            value=algebra.add(answer.get(counts,0),value)
            if value:answer[counts]=value
            else:answer.pop(counts,None)
    return sorted((c,counts) for counts,c in answer.items()),acquisition_terms


def acquire_trace_factors(ring,width,t,variables,predicates):
    """Predicates are (ordinary polynomial terms, bit target).

    Coefficients need only REDUCE into the width-bit subfield. The partial
    lifted trace may have higher digits outside the base integers; its
    residue is still Boolean, which is all the integral projector needs.
    """
    if width<1 or t<1 or width%t or ring.degree%width or variables<1:
        raise ValueError('Incompatible trace dimensions')
    columns=width//t;dimension=variables*columns
    algebra=DividedPowerRing(ring,len(predicates),2*ring.precision-1)
    factors=[];degrees=[];acquisition=[]
    for axis,(terms,bit) in enumerate(predicates):
        if bit not in (0,1):raise ValueError('A trace target must be a bit')
        degree=0
        for c,exponents in terms:
            if c<0 or c.bit_length()>ring.precision*ring.degree:
                raise ValueError('Trace coefficients must be base-ring scalars')
            if len(exponents)!=variables or any(e<0 for e in exponents):
                raise ValueError('Bad ordinary trace monomial')
            if ring.mod_two(ring.frobenius(c,width))!=ring.mod_two(c):
                raise ValueError('Trace coefficient residue lies outside the summation subfield')
            if c:degree=max(degree,sum(exponents))
        degrees.append(degree)
        for row in range(t):
            packet={}
            if row==0 and bit:packet[(0,)*dimension]=ring.integer(-bit)
            for j in range(columns):
                shift=row+t*j
                for c,exponents in terms:
                    counts=tuple(exponents[i] if column==j else 0
                                 for i in range(variables) for column in range(columns))
                    value=ring.add(packet.get(counts,0),ring.frobenius(c,shift))
                    if value:packet[counts]=value
                    else:packet.pop(counts,None)
            acquired,count=exponential_factor(algebra,axis,packet,dimension)
            factors.append((row,acquired));acquisition.append(count)
    return algebra,factors,dict(trace_degrees=degrees,acquisition_terms=acquisition,
                               trace_degree_budget=(2*ring.precision-1)*t*sum(degrees))


def trace_constrained_sum(ring,width,t,variables,predicates,
                          numerators=(),denominators=(),nonzero=False):
    """Exact ring sum under polynomial field-trace bit predicates.

    Numerator/denominator factor coefficients are embedded BASE R scalars.
    General coefficient-algebra denominators are deliberately rejected.
    """
    predicates=list(predicates);numerators=list(numerators);denominators=list(denominators)
    algebra,factors,metadata=acquire_trace_factors(ring,width,t,variables,predicates)
    for row,terms in list(numerators)+list(denominators):
        if any(c<0 or c.bit_length()>algebra.coefficient_bits for c,counts in terms):
            raise ValueError('Weights must have scalar base-ring coefficients')
    answer,backend=rational_sum(algebra,width,t,variables,list(numerators)+factors,
                                list(denominators),nonzero)
    value=algebra.indicator_functional(answer)
    budget=(sum(map(factor_degree,numerators))+
            ((2*ring.precision+1)*ring.degree-1)*sum(map(factor_degree,denominators))+
            metadata['trace_degree_budget'])
    assert backend['max_degree_budget']<=budget
    invoice=dict(precision=ring.precision,coefficient_degree=ring.degree,width=width,
                 row_modulus=t,variables=variables,relative_columns=width//t,
                 predicates=len(predicates),algebra_dimension=algebra.dimension,
                 divided_power_limit=algebra.limit,max_degree_budget=budget,
                 actual_backend=backend,**metadata)
    return value,invoice


class ResidueAdapter:
    def __init__(self,ring):self.ring=ring;self.degree=ring.degree;self.size=1 << ring.degree
    def power(self,a,exponent):
        lift=self.ring.from_polynomial([(a >> i)&1 for i in range(self.degree)])
        return residue_bits(self.ring,self.ring.power(lift,exponent))


def affine_trace_sum(ring,width,t,variables,constraints,
                     numerators=(),denominators=(),nonzero=False):
    """Reuse the earlier binary row reduction; r is its actual affine rank."""
    rows=[([residue_bits(ring,c) for c in coefficients],bit)
          for coefficients,bit in constraints]
    independent=independent_constraints(ResidueAdapter(ring),width,variables,rows)
    if independent is None:return 0,dict(inconsistent=True,rank=None,queries=0)
    predicates=[]
    for coefficients,bit in independent:
        terms=[]
        for i,c in enumerate(coefficients):
            if c:
                counts=tuple(int(i==j) for j in range(variables))
                terms.append((ring.teichmuller(c),counts))
        predicates.append((terms,bit))
    value,invoice=trace_constrained_sum(ring,width,t,variables,predicates,
                                       numerators,denominators,nonzero)
    invoice.update(inconsistent=False,rank=len(independent))
    return value,invoice
