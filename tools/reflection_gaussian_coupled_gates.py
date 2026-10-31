#!/usr/bin/env python3
"""Growing cross-leg Gaussian gates, including exact denominator and mean costs.

The field is F2. Coefficients are exact rational functions; optional pressure
optimization uses floating real mode heights. See coupled-gaussian-gates.md.
"""
from fractions import Fraction
from itertools import product
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from check_rectangular_defect_support import Rat, ZERO, ONE, Z, pmul, pdivmod, pgcd
from check_reflection_unrestricted_gaussian import closed_matrix
from reflection_dual_rail_gaussian import tensor
from reflection_unrestricted_gaussian import ReflectionGaussianPressure
from native_matrix_smith_f4 import MatrixOracle
from reflection_gaussian_physical_frames import norm, add, mul, coefficient


def local_order(value, infinity=False):
    assert value
    return (value.d.bit_length()-value.n.bit_length()) if infinity else value.order()


def selected_pairs(M, gate):
    m=M//2
    first,second=gate['legs']
    assert first!=second and 0<=first<4 and 0<=second<4
    matching=gate.get('matching','identity')
    if matching=='identity':
        return [(first*m+i,second*m+i) for i in range(m)]
    if matching=='reverse':
        return [(first*m+i,second*m+(m-1-i)) for i in range(m)]
    if matching=='all':
        return [(first*m+i,second*m+j) for i in range(m) for j in range(m)]
    if matching=='single':
        i,j=gate.get('i',0),gate.get('j',0)
        assert 0<=i<m and 0<=j<m
        return [(first*m+i,second*m+j)]
    raise ValueError(matching)


def kinds(gate):
    result=tuple(gate.get('kinds',('create','create')))
    assert len(result)==2 and all(kind in ('create','annihilate') for kind in result)
    return result


def quadratic_update(matrix, rows, columns, first, second, first_kind, second_kind, g):
    """Apply one physical quadratic factor; return B' and its determinant factor d."""
    M=len(rows)
    size=2*M
    flip=lambda index:(index+M)%size
    row_pos={physical:index for index,physical in enumerate(rows)}
    col_pos={physical:index for index,physical in enumerate(columns)}
    result=[list(line) for line in matrix]
    if first in row_pos and second in row_pos:
        i=row_pos[first]+M*(first_kind=='annihilate')
        j=row_pos[second]+M*(second_kind=='annihilate')
        result[i]=[a+g*b for a,b in zip(matrix[i],matrix[flip(j)])]
        result[j]=[a+g*b for a,b in zip(matrix[j],matrix[flip(i)])]
        return result,ONE
    if first in col_pos and second in col_pos:
        # Column physical action is transposed: C^T=A and A^T=C.
        i=col_pos[first]+M*(first_kind=='create')
        j=col_pos[second]+M*(second_kind=='create')
        for r in range(size):
            result[r][flip(j)]=matrix[r][flip(j)]+g*matrix[r][i]
            result[r][flip(i)]=matrix[r][flip(i)]+g*matrix[r][j]
        return result,ONE
    if first not in row_pos:
        first,second=second,first
        first_kind,second_kind=second_kind,first_kind
    i=row_pos[first]+M*(first_kind=='annihilate')
    j=col_pos[second]+M*(second_kind=='create')
    d=ONE+g*matrix[flip(i)][j]
    assert d and (d.n.bit_count()&1)==1 and (d.d.bit_count()&1)==1
    ratio=g/d
    for r in range(size):
        factor=ratio*matrix[r][j]
        for c in range(size):
            result[r][c]=matrix[r][c]+factor*matrix[flip(i)][c]
    result[i][flip(j)]=result[i][flip(j)]+g
    return result,d


def coupled_matrix(M, side, gates):
    matrix,rows,columns=closed_matrix(M,side)
    mean=Fraction(M*M if side==0 else 0)
    factors=[]
    for gate in gates:
        a,b=gate['powers']
        assert a!=b and gate.get('gamma',1)==1
        first_kind,second_kind=kinds(gate)
        g=Z**a+Z**b
        for first,second in selected_pairs(M,gate):
            matrix,d=quadratic_update(matrix,rows,columns,first,second,
                                      first_kind,second_kind,g)
            order=local_order(d,bool(side))
            mean+=Fraction(order,2)
            factors.append({'pair':[first,second],'d':list(d.pair()),'order':order})
    return matrix,rows,columns,mean,factors


class RationalOracle:
    """Clear a common polynomial denominator and restore it at both places."""
    def __init__(self,matrix):
        denominator=1
        for row in matrix:
            for value in row:
                quotient,remainder=pdivmod(value.d,pgcd(denominator,value.d))
                assert not remainder
                denominator=pmul(denominator,quotient)
        values=[]
        for row in matrix:
            line=[]
            for value in row:
                quotient,remainder=pdivmod(denominator,value.d)
                assert not remainder
                line.append((pmul(value.n,quotient),0))
            values.append(line)
        bits=max(max(value).bit_length() for row in values for value in row)
        self.precision=max(512,1<<(max(bits,8*len(matrix))-1).bit_length())
        self.values=values
        self.size=len(matrix)
        self.oracle=MatrixOracle(values,self.size,precision=self.precision)
        self.orders=((denominator&-denominator).bit_length()-1,1-denominator.bit_length())
        self.denominator=denominator
        self.bits=bits

    def profile(self,row_weights=None,column_weights=None,scale=1.,infinity=False):
        while True:
            try:
                result=self.oracle.profile(row_weights,column_weights,scale,infinity)
                break
            except ArithmeticError:
                if self.precision>=8192:
                    raise
                self.precision*=2
                self.oracle=MatrixOracle(self.values,self.size,precision=self.precision)
        shift=self.orders[int(bool(infinity))]
        result['orders']=[value-shift*scale for value in result['orders']]
        result['base_prefix']=[value-shift*(i+1) for i,value in enumerate(result['base_prefix'])]
        return result


class CoupledPressure(ReflectionGaussianPressure):
    def __init__(self,M,gates):
        super().__init__(M)
        self.gates=gates
        self.cuts=[]
        self.means=[]
        self.exact_means=[]
        self.factors=[]
        for side in (0,1):
            matrix,rows,cols,mean,factors=coupled_matrix(M,side,gates)
            self.cuts.append((RationalOracle(matrix),np.array(rows),np.array(cols)))
            self.exact_means.append(mean)
            self.means.append(float(mean))
            self.factors.append(factors)

    def evaluate(self,x,gradient=False):
        x=np.asarray(x,dtype=float)
        ell=x[-1]
        mean=sum(self.means)
        log_b=-mean*ell-.5*x[:-1].sum()
        deriv=np.r_[np.full(2*self.size,-.5),-mean]
        records=[]
        for side,(oracle,row,col) in enumerate(self.cuts):
            h=x[side*self.size:(side+1)*self.size]
            result=oracle.profile(np.r_[h[row],-h[row]],np.r_[h[col],-h[col]],ell,
                                  infinity=bool(side))
            orders=np.array(result['orders'])
            assert np.max(np.abs(orders+orders[::-1]))<1e-6
            log_b+=.5*np.logaddexp(orders/2,-orders/2).sum()
            coeff=.25*np.tanh(orders/2)
            base=np.diff([0]+result['base_prefix'])
            for c,r,j,b in zip(coeff,result['rows'],result['cols'],base):
                deriv[side*self.size+row[r%self.interval]]+=c*(1 if r<self.interval else -1)
                deriv[side*self.size+col[j%self.interval]]+=c*(1 if j<self.interval else -1)
                deriv[-1]+=c*b
            records.append(result)
        hh=x[:-1].reshape(2,4,self.modes)
        physical=np.array([hh[:,0,:].sum(axis=0)+hh[:,3,:].sum(axis=0),
                           hh[:,1,:].sum(axis=0)+hh[:,2,:].sum(axis=0)])
        def equation(r):
            return log_b-r*np.logaddexp(0,-physical/r).sum()
        if equation(1e-8)<=0:
            raise ArithmeticError('No positive stable root')
        high=4.
        while equation(high)>0:
            high*=2
        rate=brentq(equation,1e-8,high,xtol=1e-13)
        p=expit(-physical/rate)
        entropy=np.logaddexp(0,-physical/rate).sum()+(physical/rate*p).sum()
        if entropy<1e-8:
            raise ArithmeticError('Small physical entropy')
        for side in range(2):
            for leg in range(4):
                j=side*self.size+leg*self.modes
                deriv[j:j+self.modes]+=p[0 if leg in (0,3) else 1]
        self.queries+=1
        self.last={'rho':float(rate/2),'critical_pressure':float(equation(math.sqrt(2))),
                   'entropy':float(entropy),'parameters':x.tolist(),'profiles':records,
                   'means':self.means}
        return (rate/2,deriv/(2*entropy)) if gradient else rate/2


def literal_tensor(M,gates):
    """Independent original coefficient tensor; intended only for small controls."""
    values=[norm(a,b,0) for a,b in tensor(M)['values']]
    for gate in gates:
        assert gate.get('gamma',1)==1
        g=coefficient(gate)
        first_kind,second_kind=kinds(gate)
        for first,second in selected_pairs(M,gate):
            mask=(1<<first)|(1<<second)
            occupied=(int(first_kind=='annihilate')<<first)|(int(second_kind=='annihilate')<<second)
            for state in range(len(values)):
                if state&mask==occupied:
                    values[state^mask]=add(values[state^mask],mul(g,values[state]))
    clearing=max(0,-min(s for a,b,s in values if a|b))
    polynomial=[(a<<(s+clearing),b<<(s+clearing)) for a,b,s in values]
    return polynomial,clearing


def degree_bound(M,gates):
    """Each uniform two-leg block uses at most m quadratic factors per term."""
    m=M//2
    negative=sum(max(0,-min(gate['powers'])) for gate in gates)*m
    positive=sum(max(0,max(gate['powers'])) for gate in gates)*m
    return {'clearing_bound':negative,'cleared_degree_bound':3*M*M//2+negative+positive}
