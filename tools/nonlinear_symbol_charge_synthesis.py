#!/usr/bin/env python3
"""Exact nonlinear symbol charges and their complete-entry affine family."""
from fractions import Fraction
from typed_charge_lattice_synthesis import (TypedChargeLattice,tensor,chi_power,
    regular_smith_basis,polynomial_multiply,polynomial_degrees,matrix_multiply,cut_location)


def cube_lipschitz(values):
    values=tuple(map(Fraction,values));q=len(values)
    if q<2 or q&(q-1):raise ValueError('Use a complete power-of-two occupation alphabet.')
    m=q.bit_length()-1
    return max(abs(values[state]-values[state^(1<<i)]) for state in range(q) for i in range(m))


def four_symbol_coordinates(values):
    values=tuple(map(Fraction,values))
    if len(values)!=4:raise ValueError('Use four physical-symbol values.')
    a,b,c=(value-values[0] for value in values[1:])
    return c/4,(a-b)/4,(a+b-c)/4


class SymbolChargeLattice(TypedChargeLattice):
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
            raise ValueError('Supply two complete physical-symbol potential lists.')
        if len(charge_weights)!=2 or any(len(row)!=self.q for row in charge_weights):
            raise ValueError('Supply q entries in each paired symbol potential.')
        fx,fy=(tuple(map(Fraction,row)) for row in charge_weights)
        fx=tuple(v-fx[0] for v in fx);fy=tuple(v-fy[0] for v in fy)
        if cube_lipschitz(fx)+cube_lipschitz(fy)>2:
            raise ValueError('The exact symbol-charge criterion is Lip(f)+Lip(g)<=2.')
        self.symbol_weights=(fx,fy);self.charge_weights=fx+fy
        self.charges=[fx[x]+fy[y]-fy[u]-fx[v] for x,y,u,v in self.labels]
        from check_rectangular_defect_support import pdivmod
        for value,energy in zip(self.source,self.charges):
            if not value:continue
            order=0
            while True:
                quotient,remainder=pdivmod(value,3)
                if remainder:break
                value=quotient;order+=1
            if order+energy<0:raise ValueError('The exact original coefficient charge inequality fails.')
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

