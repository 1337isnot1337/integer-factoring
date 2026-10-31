#!/usr/bin/env python3
"""Exact complete-entry affine perturbations in the typed charge lattice.

For rational paired mode weights with max|a_i|+max|b_i|<=2, the baseline
satisfies ord_chi(T_lambda)+E(lambda)>=0, where chi=1+w. Perturbations of
order floor(-E)+1 preserve every paired-type polynomial crossing endpoint.
An optional extension also permits equality at strictly negative charge.
These are sufficient local conditions, not a classification of endpoints.
"""
from itertools import product
from fractions import Fraction
from reflection_dual_rail_gaussian import tensor
from weyl_lattice_filter_synthesis import (regular_smith_basis,polynomial_multiply,
    polynomial_degrees,matrix_multiply)
from full_alphabet_two_cut_lattices import (FullAlphabetTwoCutLattice,
    CompressedFullAlphabetFiber,cut_location)
from native_matrix_smith_f4 import MatrixOracle


def chi_power(power):
    result=0;sub=power
    while True:
        result|=1<<sub
        if sub==0:return result
        sub=(sub-1)&power


class TypedChargeLattice(FullAlphabetTwoCutLattice):
    def __init__(self,M,heights,raise_cap=10,extra_degree=0,charge_weights=None,
                 negative_boundary=False,reference_values=None):
        if M!=4:raise ValueError('The finite affine implementation currently supports M4; the theorem is uniform in M.')
        if (not isinstance(raise_cap,int) or raise_cap<1 or
                not isinstance(extra_degree,int) or extra_degree<0):
            raise ValueError('Use a positive target allowance and nonnegative extra polynomial degree.')
        self.M=M;self.m=M//2;self.q=1<<self.m;self.size=self.q**2
        self.clearing=0
        self.heights=tuple(map(tuple,heights))
        if len(heights)!=8 or any(len(row)!=self.q or any(not isinstance(v,int) for v in row) for row in heights):
            raise ValueError('Eight full integer height rows are required.')
        data=tensor(M);self.source=tuple(a for a,b in data['values'])
        assert all(b==0 for a,b in data['values'])
        self.source_degree=max(v.bit_length()-1 for v in self.source)
        self.labels=tuple(tuple((index//self.q**i)%self.q for i in range(4)) for index in range(self.q**4))
        if charge_weights is None:
            charge_weights=[Fraction(i,self.m-1) for _ in range(2) for i in range(self.m)]
        charge_weights=tuple(map(Fraction,charge_weights))
        if len(charge_weights)!=2*self.m or max(map(abs,charge_weights[:self.m]))+max(map(abs,charge_weights[self.m:]))>2:
            raise ValueError('Paired mode weights must satisfy max|a_i|+max|b_i|<=2.')
        self.charge_weights=charge_weights
        first=charge_weights[:self.m];second=charge_weights[self.m:]
        charge=first+second+tuple(-v for v in second)+tuple(-v for v in first)
        self.charges=[sum(charge[i] for i in range(2*M) if state>>i&1) for state in range(self.q**4)]
        self.orders_at_chi=[(-energy)//1+1 for energy in self.charges]
        self.negative_boundary=negative_boundary
        if negative_boundary:
            self.orders_at_chi=[k-int(energy<0 and energy.denominator==1)
                                for k,energy in zip(self.orders_at_chi,self.charges)]
        self.pole=max(0,-min(self.orders_at_chi))
        self.degree=self.source_degree+self.pole+extra_degree
        self.factors=[chi_power(self.pole+k) for k in self.orders_at_chi]
        self.widths=[max(0,self.degree-(self.pole+k)+1) for k in self.orders_at_chi]
        self.offsets=[];offset=0
        for width in self.widths:self.offsets.append(offset);offset+=width
        self.variables=offset
        self.chi=chi_power(self.pole)
        self.base=tuple(polynomial_multiply(v,self.chi) for v in self.source)
        self.reference_values=tuple((v,0) for v in self.base) if reference_values is None else tuple(map(tuple,reference_values))
        if len(self.reference_values)!=self.q**4 or any(len(v)!=2 or v[1] or v[0]<0 or
                v[0].bit_length()>self.degree+1 for v in self.reference_values):
            raise ValueError('A frame reference must be a complete binary numerator within the degree cap.')
        from check_rectangular_defect_support import pdivmod
        if any(pdivmod(v[0]^base,factor)[1] for v,base,factor in zip(self.reference_values,self.base,self.factors)):
            raise ValueError('The reference must belong to the same typed charge lattice.')
        self.bases=[];self.orders=[];self.scalar_shifts=[];self.equations=[];self.caps=[]
        self.equation_columns=[]
        for side in range(2):
            h=self.heights[side*4:side*4+4]
            scalar=-sum(min(row) for row in h)+(self.degree if side==0 else 0)
            base_scalar=-sum(min(row) for row in h)+(self.source_degree if side==0 else 0)
            self.scalar_shifts.append(scalar)
            weights=[sum(row[label] for row,label in zip(h,bits)) for bits in self.labels]
            matrix=[[0]*self.size for _ in range(self.size)]
            for entry,value in enumerate(self.source):
                i,j=cut_location(self.labels[entry],self.q,side)
                for power in polynomial_degrees(value):
                    exponent=base_scalar+weights[entry]+(-power if side==0 else power)
                    assert exponent>=0;matrix[i][j]^=1<<exponent
            if reference_values is None:
                orders,left,right,diagonal=regular_smith_basis(matrix)
                shift=extra_degree if side==0 else 0
                reference=[[polynomial_multiply(value,self.chi)<<shift for value in row] for row in diagonal]
                orders=tuple(v+shift for v in orders)
            else:
                from finite_precision_smith_frames import bounded_smith_frame
                frame_source=[[0]*self.size for _ in range(self.size)]
                origin=[[0]*self.size for _ in range(self.size)]
                for entry,pair in enumerate(self.reference_values):
                    i,j=cut_location(self.labels[entry],self.q,side)
                    for target,value in ((frame_source,pair[0]),(origin,self.base[entry])):
                        for power in polynomial_degrees(value):
                            exponent=scalar+weights[entry]+(-power if side==0 else power)
                            assert exponent>=0;target[i][j]^=1<<exponent
                orders,left,right,_=bounded_smith_frame(frame_source,raise_cap)
                reference=matrix_multiply(matrix_multiply(left,origin),right)
            self.bases.append((left,right));self.orders.extend(orders)
            caps=[v+raise_cap for v in orders];self.caps.extend(caps)
            tables=[[[0 for _ in range(cap)] for _ in range(self.size)] for cap in caps]
            for entry,bits in enumerate(self.labels):
                i,j=cut_location(bits,self.q,side);width=self.widths[entry]
                if not width:continue
                factor=self.factors[entry];factor_degree=factor.bit_length()-1
                # Each chi power is self-reciprocal.  Anchor the highest P
                # exponent at infinity and the lowest exponent at zero.
                anchor=scalar+weights[entry]-(factor_degree+width-1) if side==0 else scalar+weights[entry]
                assert anchor>=0
                for row in range(self.size):
                    if not left[row][i]:continue
                    cap=caps[row]
                    for col in range(self.size):
                        if not right[j][col]:continue
                        outer=polynomial_multiply(polynomial_multiply(left[row][i],right[j][col]),factor)<<anchor
                        for power in polynomial_degrees(outer&((1<<cap)-1)):
                            for step in range(min(width,cap-power)):
                                variable=self.offsets[entry]+(width-1-step if side==0 else step)
                                tables[row][col][power+step]^=1<<variable
            for row,cap in enumerate(caps):
                equations=[];columns=[]
                for power in range(cap):
                    for col in range(self.size):
                        value=tables[row][col][power];rhs=(reference[row][col]>>power)&1
                        if value or rhs:
                            equations.append((power,value,rhs));columns.append(col)
                self.equations.append(equations);self.equation_columns.append(columns)
        self.orders=tuple(self.orders);self.caps=tuple(self.caps)

    def seed(self,low,high=0):
        if any(not isinstance(v,int) or v<0 or v.bit_length()>self.variables for v in (low,high)):
            raise ValueError('Coefficient planes must fit the complete affine variable list.')
        result=[]
        for entry,base in enumerate(self.base):
            mask=(1<<self.widths[entry])-1;offset=self.offsets[entry]
            a=(low>>offset)&mask;b=(high>>offset)&mask;factor=self.factors[entry]
            result.append((base^polynomial_multiply(a,factor),polynomial_multiply(b,factor)))
        return tuple(result)

    def as_seed(self,values):
        return {'values':values,'dimensions':[self.q]*4,'q_x':self.q,'q_u':self.q,
                'degree':self.degree,'typed_denominator_degree':self.pole,'denominator_degree':self.pole,
                'clearing':self.clearing,'strict_endpoint':False,'typed_endpoint':True,
                'common_w_clearing':self.clearing,'endpoint_pole_order':self.pole,
                'effective_cleared_degree':self.degree-self.pole}

    def profiles(self,values):
        profiles=[]
        for side in range(2):
            matrix=[[(0,0)]*self.size for _ in range(self.size)]
            for entry,pair in enumerate(values):
                i,j=cut_location(self.labels[entry],self.q,side);matrix[i][j]=pair
            h=self.heights[side*4:side*4+4]
            if side==0:
                rw=[h[0][a]+h[2][b] for a,b in product(range(self.q),repeat=2)]
                cw=[h[1][a]+h[3][b] for a,b in product(range(self.q),repeat=2)]
            else:
                rw=[h[0][a]+h[1][b] for a,b in product(range(self.q),repeat=2)]
                cw=[h[2][a]+h[3][b] for a,b in product(range(self.q),repeat=2)]
            raw=MatrixOracle(matrix,self.size).profile(rw,cw,infinity=side==0)['orders']
            correction=self.clearing+self.pole if side==0 else -self.clearing
            assert all(abs(value-round(value))<1e-8 for value in raw)
            profiles.append(tuple(round(value)+correction for value in raw))
        return tuple(profiles)

    def check_row_witness(self,values,targets):
        for side,(left,right) in enumerate(self.bases):
            h=self.heights[4*side:4*side+4];scalar=self.scalar_shifts[side]
            for plane in range(2):
                matrix=[[0]*self.size for _ in range(self.size)]
                for entry,pair in enumerate(values):
                    i,j=cut_location(self.labels[entry],self.q,side)
                    weight=sum(row[label] for row,label in zip(h,self.labels[entry]))
                    for power in polynomial_degrees(pair[plane]):
                        exponent=scalar+weight+(-power if side==0 else power)
                        assert exponent>=0
                        matrix[i][j]^=1<<exponent
                transformed=matrix_multiply(matrix_multiply(left,matrix),right)
                for i,row in enumerate(transformed):
                    assert all(value&((1<<targets[side*self.size+i])-1)==0 for value in row)
        profiles=self.profiles(values)
        for side in range(2):
            correction=self.clearing+self.pole if side==0 else -self.clearing
            bounds=sorted(targets[side*self.size+i]-self.scalar_shifts[side]+correction for i in range(self.size))
            for k in range(1,self.size+1):assert sum(profiles[side][:k])>=sum(bounds[:k])
        return profiles

    def check_witness(self,values,targets):
        # The prescribed entrywise chi divisibility is checked independently
        # of the linear coefficient tables before testing the two row lattices.
        from check_rectangular_defect_support import pdivmod
        for entry,(low,high) in enumerate(values):
            assert not pdivmod(low^self.base[entry],self.factors[entry])[1]
            assert not pdivmod(high,self.factors[entry])[1]
            assert (self.orders_at_chi[entry]+self.charges[entry]>0 or
                    self.negative_boundary and self.charges[entry]<0 and
                    self.orders_at_chi[entry]+self.charges[entry]==0)
            assert max(low,high).bit_length()<=self.degree+1
        return self.check_row_witness(values,targets)



def natural_heights(M):
    m=M//2;q=1<<m
    energy=[sum(4*(M-1-2*i) for i in range(m) if state>>i&1) for state in range(q)]
    half=[v//2 for v in energy];zero=[0]*q
    return [zero,zero,half,half,half,half,zero,zero],energy
