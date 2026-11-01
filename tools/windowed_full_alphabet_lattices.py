#!/usr/bin/env python3
"""Windowed constructor retaining the historical complete-affine source unchanged.

The existing class supplies solving, sampling and literal witness validation.
This constructor builds identical low-coefficient equations using bounded
regular frames. The original source has pinned historical search receipts.
"""
from itertools import product

from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice, cut_location, shifted
from windowed_binary_smith import windowed_regular_smith_basis
from weyl_lattice_filter_synthesis import polynomial_degrees, polynomial_multiply


class WindowedFullAlphabetTwoCutLattice(FullAlphabetTwoCutLattice):
    def __init__(self,values,q,degree,heights,raise_cap=8):
        if q<2 or len(values)!=q**4 or len(heights)!=8:
            raise ValueError('Supply all q^4 entries and eight physical height lists.')
        if any(len(row)!=q or any(not isinstance(value,int) for value in row) for row in heights):
            raise ValueError('Each physical height list has q integer entries.')
        if any(len(pair)!=2 or any(not isinstance(v,int) or v<0 for v in pair)
               or pair[1] or pair[0].bit_length()>degree+1 for pair in values):
            raise ValueError('The frame source must be binary and fit the polynomial degree.')
        if degree<1 or raise_cap<1:raise ValueError('Use a positive degree and row-target allowance.')
        self.frame_backend='windowed';self.frame_windows=[]
        self.q=q;self.size=q*q;self.degree=degree;self.variables=q**4*degree
        self.base=tuple(pair[0] for pair in values);self.heights=tuple(map(tuple,heights))
        self.labels=tuple(tuple((index//q**i)%q for i in range(4)) for index in range(q**4))
        for entry,value in enumerate(self.base):
            x,y,u,v=self.labels[entry]
            if value.bit_count()%2!=int(x==v and y==u):
                raise ValueError('The frame source must have the literal crossing endpoint.')
        self.bases=[];self.orders=[];self.scalar_shifts=[];self.equations=[];self.caps=[]
        variable_mask=(1<<degree)-1
        for side in range(2):
            height_rows=self.heights[4*side:4*side+4]
            scalar=-sum(min(row) for row in height_rows)+(degree if side==0 else 0)
            self.scalar_shifts.append(scalar)
            weights=[sum(row[label] for row,label in zip(height_rows,bits)) for bits in self.labels]
            source=[[0]*self.size for _ in range(self.size)]
            for entry,base in enumerate(self.base):
                i,j=cut_location(self.labels[entry],q,side)
                for power in polynomial_degrees(base):
                    exponent=scalar+weights[entry]+(-power if side==0 else power)
                    assert exponent>=0
                    source[i][j]^=1<<exponent
            orders,left,right,diagonal,window=windowed_regular_smith_basis(source,raise_cap)
            self.frame_windows.append(window)
            if len(orders)!=self.size:raise ArithmeticError('The strict endpoint must give full generic rank.')
            self.bases.append((left,right));self.orders.extend(orders)
            caps=[order+raise_cap for order in orders];self.caps.extend(caps)
            tables=[[[0 for _ in range(cap)] for _ in range(self.size)] for cap in caps]
            for entry,bits in enumerate(self.labels):
                i,j=cut_location(bits,q,side)
                offset=entry*degree
                for a in range(self.size):
                    if not left[a][i]:continue
                    cap=caps[a];truncation=(1<<cap)-1
                    for b in range(self.size):
                        if not right[j][b]:continue
                        outer=polynomial_multiply(left[a][i],right[j][b])
                        outer=(outer^(outer<<1))&truncation
                        if not outer:continue
                        if side==1:
                            width=outer.bit_length();reverse=int(format(outer,'0'+str(width)+'b')[::-1],2)
                            for power in range(cap):
                                at=power-scalar-weights[entry]
                                section=shifted(reverse,at-width+1)&variable_mask
                                tables[a][b][power]^=section<<offset
                        else:
                            for power in range(cap):
                                at=power-scalar-weights[entry]+1
                                section=shifted(outer,-at)&variable_mask
                                tables[a][b][power]^=section<<offset
            for a,cap in enumerate(caps):
                rows=[]
                for power in range(cap):
                    for b in range(self.size):
                        equation=tables[a][b][power];rhs=(diagonal[a][b]>>power)&1
                        if equation or rhs:rows.append((power,equation,rhs))
                self.equations.append(rows)
        self.orders=tuple(self.orders);self.caps=tuple(self.caps)
