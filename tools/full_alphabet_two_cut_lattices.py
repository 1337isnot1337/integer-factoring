#!/usr/bin/env python3
"""Complete-entry affine crossing synthesis in fixed regular binary frames.

The physical alphabet can have arbitrary finite size. Every coefficient
entry is free over F4, subject to its exact endpoint and degree bound.
The reference used to construct the regular frames must have binary
coefficients. Frames are sufficient row-lattice constraints, not a
parameterization of all tensors with a specified Smith profile.
"""
from __future__ import annotations

from itertools import product
import random

from native_matrix_smith_f4 import MatrixOracle
from weyl_lattice_filter_synthesis import (
    regular_smith_basis,polynomial_multiply,polynomial_degrees,matrix_multiply)


def index_of(bits,q):
    x,y,u,v=bits
    return x+q*y+q*q*u+q*q*q*v


def cut_location(bits,q,side):
    x,y,u,v=bits
    return (q*x+u,q*y+v) if side==0 else (q*x+y,q*u+v)


def shifted(value,shift):
    return value<<shift if shift>=0 else value>>(-shift)


class FullAlphabetTwoCutLattice:
    def __init__(self,values,q,degree,heights,raise_cap=8):
        if q<2 or len(values)!=q**4 or len(heights)!=8:
            raise ValueError('Supply all q^4 entries and eight physical height lists.')
        if any(len(row)!=q or any(not isinstance(value,int) for value in row) for row in heights):
            raise ValueError('Each physical height list has q integer entries.')
        if any(len(pair)!=2 or any(not isinstance(v,int) or v<0 for v in pair)
               or pair[1] or pair[0].bit_length()>degree+1 for pair in values):
            raise ValueError('The frame source must be binary and fit the polynomial degree.')
        if degree<1 or raise_cap<1:raise ValueError('Use a positive degree and row-target allowance.')
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
            orders,left,right,diagonal=regular_smith_basis(source)
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

    def solve(self,targets):
        if len(targets)!=2*self.size or any(t<0 or t>cap for t,cap in zip(targets,self.caps)):
            raise ValueError('The requested row targets exceed the completely constructed coefficient window.')
        pivots={}
        for target,equations in zip(targets,self.equations):
            for power,row,rhs in equations:
                if power>=target:break
                while row:
                    pivot=row.bit_length()-1
                    if pivot not in pivots:
                        pivots[pivot]=(row,rhs);break
                    old,value=pivots[pivot];row^=old;rhs^=value
                if not row and rhs:return None
        return pivots

    def assignment(self,pivots,free_bits=0,homogeneous=False):
        value=free_bits
        for pivot,(row,rhs) in sorted(pivots.items()):
            assert not (value>>pivot)&1
            value|=(((row&value).bit_count()&1)^(0 if homogeneous else rhs))<<pivot
        return value

    def seed(self,low,high=0):
        mask=(1<<self.degree)-1
        result=[]
        for entry,base in enumerate(self.base):
            a=(low>>(entry*self.degree))&mask;b=(high>>(entry*self.degree))&mask
            result.append((base^a^(a<<1),b^(b<<1)))
        return tuple(result)

    def sample(self,pivots,rng):
        free=[index for index in range(self.variables) if index not in pivots]
        planes=[sum(rng.randrange(2)<<index for index in free) for _ in range(2)]
        return self.seed(self.assignment(pivots,planes[0]),self.assignment(pivots,planes[1],True))

    def profiles(self,values):
        profiles=[]
        for side in range(2):
            matrix=[[(0,0)]*self.size for _ in range(self.size)]
            for entry,pair in enumerate(values):
                i,j=cut_location(self.labels[entry],self.q,side);matrix[i][j]=pair
            h=self.heights[4*side:4*side+4]
            if side==0:
                rw=[h[0][a]+h[2][b] for a,b in product(range(self.q),repeat=2)]
                cw=[h[1][a]+h[3][b] for a,b in product(range(self.q),repeat=2)]
            else:
                rw=[h[0][a]+h[1][b] for a,b in product(range(self.q),repeat=2)]
                cw=[h[2][a]+h[3][b] for a,b in product(range(self.q),repeat=2)]
            raw=MatrixOracle(matrix,self.size).profile(rw,cw,infinity=side==0)['orders']
            assert all(abs(value-round(value))<1e-8 for value in raw)
            profiles.append(tuple(round(value) for value in raw))
        return tuple(profiles)

    def check_witness(self,values,targets):
        for entry,pair in enumerate(values):
            x,y,u,v=self.labels[entry]
            assert ((pair[0].bit_count()&1)|((pair[1].bit_count()&1)<<1))==int(x==v and y==u)
        for side,(left,right) in enumerate(self.bases):
            h=self.heights[4*side:4*side+4];scalar=self.scalar_shifts[side]
            for plane in range(2):
                matrix=[[0]*self.size for _ in range(self.size)]
                for entry,pair in enumerate(values):
                    i,j=cut_location(self.labels[entry],self.q,side)
                    weight=sum(row[label] for row,label in zip(h,self.labels[entry]))
                    for power in polynomial_degrees(pair[plane]):
                        matrix[i][j]^=1<<(scalar+weight+(-power if side==0 else power))
                transformed=matrix_multiply(matrix_multiply(left,matrix),right)
                for i,row in enumerate(transformed):
                    assert all(value&((1<<targets[side*self.size+i])-1)==0 for value in row)
        profiles=self.profiles(values)
        for side in range(2):
            bounds=sorted(targets[side*self.size+i]-self.scalar_shifts[side] for i in range(self.size))
            for k in range(1,self.size+1):assert sum(profiles[side][:k])>=sum(bounds[:k])
        return profiles


class CompressedFullAlphabetFiber:
    """Exact coordinates on a common row-lattice fiber for repeated solves."""
    def __init__(self,lattice,floors):
        self.lattice=lattice;self.floors=tuple(floors)
        self.base_pivots=lattice.solve(floors)
        if self.base_pivots is None:raise ValueError('The proposed common fiber is empty.')
        self.free=tuple(i for i in range(lattice.variables) if i not in self.base_pivots)
        self.variables=len(self.free);positions={value:i for i,value in enumerate(self.free)}
        self.equations=[]
        for row in lattice.equations:
            result=[]
            for exponent,value,rhs in row:
                reduced=0
                while value:
                    pivot=value.bit_length()-1
                    if pivot in self.base_pivots:
                        old,constant=self.base_pivots[pivot];value^=old;rhs^=constant
                    else:
                        reduced^=1<<positions[pivot];value^=1<<pivot
                if reduced or rhs:result.append((exponent,reduced,rhs))
            self.equations.append(result)

    def solve(self,targets):
        if len(targets)!=len(self.floors) or any(t<f or t>cap for t,f,cap in
            zip(targets,self.floors,self.lattice.caps)):
            raise ValueError('Targets must stay inside the completely reduced fiber window.')
        pivots={}
        for target,equations in zip(targets,self.equations):
            for power,row,rhs in equations:
                if power>=target:break
                while row:
                    pivot=row.bit_length()-1
                    if pivot not in pivots:pivots[pivot]=(row,rhs);break
                    old,value=pivots[pivot];row^=old;rhs^=value
                if not row and rhs:return None
        return pivots

    def sample(self,pivots,rng):
        free=[i for i in range(self.variables) if i not in pivots]
        result=[]
        for plane in range(2):
            value=sum(rng.randrange(2)<<i for i in free)
            for pivot,(row,rhs) in sorted(pivots.items()):
                value|=(((row&value).bit_count()&1)^(0 if plane else rhs))<<pivot
            original=0
            while value:
                bit=value&-value;value^=bit;original|=1<<self.free[bit.bit_length()-1]
            result.append(self.lattice.assignment(self.base_pivots,original,homogeneous=bool(plane)))
        return self.lattice.seed(*result)
