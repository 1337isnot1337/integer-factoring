#!/usr/bin/env python3
"""Exact local Smith orders over binary extension fields by scalar restriction.

The supplied generic rank is a promise, as in the native matrix engine.
A degree-d extension becomes a d-times larger binary matrix. Every original
Smith order occurs d times. All field-extension costs are explicit.
"""
from functools import lru_cache
from collections import Counter
from native_matrix_smith_f4 import MatrixOracle


def remainder(a,b):
    while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
    return a


def gcd(a,b):
    while b:a,b=b,remainder(a,b)
    return a


class BinaryExtension:
    def __init__(self,polynomial):
        if type(polynomial)is not int or polynomial<2:raise ValueError('Invalid binary modulus.')
        self.polynomial=polynomial;self.degree=polynomial.bit_length()-1;self.size=1<<self.degree
        if self.degree>16:raise ValueError('This implementation supports extension degrees through16.')
        x=remainder(2,polynomial);v=x
        for i in range(1,self.degree+1):
            v=self.multiply(v,v)
            if i<=self.degree//2 and gcd(v^x,polynomial)!=1:raise ValueError('The binary modulus is reducible.')
        if v!=x:raise ValueError('The binary modulus is reducible.')

    def multiply(self,a,b):
        if type(a)is not int or type(b)is not int or not 0<=a<self.size or not 0<=b<self.size:raise ValueError("Invalid field elements.")
        result=0
        while b:
            if b&1:result^=a
            b>>=1;a<<=1
            if a&self.size:a^=self.polynomial
        return result

    def power(self,a,n):
        if type(n)is not int or n<0 or type(a)is not int or not 0<=a<self.size:raise ValueError("Use a field element and nonnegative exponent.")
        result=1
        while n:
            if n&1:result=self.multiply(result,a)
            a=self.multiply(a,a);n>>=1
        return result

    def inverse(self,a):
        if not a:raise ZeroDivisionError('Zero has no inverse.')
        return self.power(a,self.size-2)

    @lru_cache(None)
    def multiplication_columns(self,a):
        return tuple(self.multiply(a,1<<j) for j in range(self.degree))


class ExtensionMatrixOracle:
    def __init__(self,matrix,field_polynomial,rank,precision=128):
        self.field=BinaryExtension(field_polynomial);self.n=len(matrix);d=self.field.degree
        if not self.n or any(len(row)!=self.n for row in matrix):raise ValueError('Expected a nonempty square polynomial matrix.')
        if type(rank)is not int or not 1<=rank<=self.n:raise ValueError('Supply a valid generic rank.')
        if d*self.n>1024:raise ValueError('Scalar restriction exceeds the native order1024 limit.')
        if any(type(c)is not int or not 0<=c<self.field.size for row in matrix for p in row for c in p):raise ValueError('Invalid field coefficient.')
        restricted=[[(0,0) for _ in range(d*self.n)] for _ in range(d*self.n)]
        for i,row in enumerate(matrix):
            for j,p in enumerate(row):
                blocks=[[0]*d for _ in range(d)]
                for k,c in enumerate(p):
                    if c:
                        for b,column in enumerate(self.field.multiplication_columns(c)):
                            for a in range(d):blocks[a][b]|=((column>>a)&1)<<k
                for a in range(d):
                    for b in range(d):restricted[d*i+a][d*j+b]=(blocks[a][b],0)
        if type(precision)is not int or precision<128 or precision>8192 or precision&(precision-1):raise ValueError("Use a power-of-two precision from128 through8192.")
        maximum=max((len(p) for row in matrix for p in row),default=0)
        while precision<maximum:precision*=2
        if precision>8192:raise ValueError("Polynomial exceeds maximum precision8192.")
        self.rank=rank;self.restricted=restricted;self.precision=precision

    def profile(self,infinity=False):
        d=self.field.degree;precision=self.precision
        while precision<=8192:
            try:
                raw=MatrixOracle(self.restricted,d*self.rank,precision=precision).profile(infinity=infinity)['orders'];break
            except ArithmeticError:precision*=2
        else:raise ArithmeticError('The promised rank was not resolved within precision8192.')
        if any(v!=round(v) for v in raw):raise ArithmeticError('Unweighted local orders must be integers.')
        counts=Counter(map(int,raw))
        if any(c%d for c in counts.values()):raise ArithmeticError('Scalar-restriction multiplicity check failed.')
        orders=sorted(v for v,c in counts.items() for _ in range(c//d))
        assert len(orders)==self.rank
        return {'orders':orders,'extension_degree':d,'original_order':self.n,'restricted_order':d*self.n,'precision':precision,'generic_rank_supplied':self.rank,'local_place':'infinity' if infinity else 'zero'}
