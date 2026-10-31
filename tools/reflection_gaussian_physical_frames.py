#!/usr/bin/env python3
"""Exact growing physical Gaussian frames over F4. See the companion proof note."""
from itertools import product
import numpy as np
from check_rectangular_defect_support import pmul
from reflection_dual_rail_gaussian import tensor,logical_heights
from reflection_unrestricted_gaussian import bogoliubov_data,ReflectionGaussianPressure
from native_matrix_smith_f4 import MatrixOracle
from relabelled_gaussian_sums import flatten

ZERO=(0,0,0)
def norm(a,b,s):
    if not(a|b):return ZERO
    shift=((a|b)&-(a|b)).bit_length()-1
    return a>>shift,b>>shift,s+shift
def add(x,y):
    if not(x[0]|x[1]):return y
    if not(y[0]|y[1]):return x
    s=min(x[2],y[2])
    return norm((x[0]<<(x[2]-s))^(y[0]<<(y[2]-s)),(x[1]<<(x[2]-s))^(y[1]<<(y[2]-s)),s)
def mul(x,y):
    if not(x[0]|x[1]) or not(y[0]|y[1]):return ZERO
    high=pmul(x[1],y[1])
    return norm(pmul(x[0],y[0])^high,pmul(x[0],y[1])^pmul(x[1],y[0])^high,x[2]+y[2])
def coefficient(gate):
    a,b=sorted(gate['powers']);assert a!=b
    gamma=gate.get('gamma',1);assert gamma in(1,2,3)
    polynomial=1|(1<<(b-a))
    return (polynomial if gamma&1 else 0,polynomial if gamma&2 else 0,a)
def pairs(m,matching):
    if matching=='reverse':return [(i,m-1-i) for i in range(m//2)]
    if matching=='adjacent':return [(i,i+1) for i in range(0,m,2)]
    if matching=='all':return [(i,j) for i in range(m) for j in range(i+1,m)]
    raise ValueError(matching)

def frame_matrix(M,side,gates):
    data,shifts,rows,cols=bogoliubov_data(M,side)
    b=[[norm(x,y,s) for (x,y),s in zip(row,powers)] for row,powers in zip(data,shifts)]
    m=M//2
    for gate in gates:
        leg=gate['leg'];kind=gate['kind'];assert kind in('create','annihilate')
        g=coefficient(gate);indices=pairs(m,gate.get('matching','reverse'))
        left=leg*m in rows;physical=rows if left else cols
        locations={j:physical.index(leg*m+j) for j in range(m)}
        target=0 if kind=='create' else M;source=M-target
        if left:
            # Sources lie in the opposite half and are unchanged by the
            # whole alternating shear, including overlapping mode pairs.
            for i,j in indices:
                first,second=locations[i],locations[j]
                b[target+first]=[add(x,mul(g,y)) for x,y in zip(b[target+first],b[source+second])]
                b[target+second]=[add(x,mul(g,y)) for x,y in zip(b[target+second],b[source+first])]
        else:
            # Acting on a tensor's column leg multiplies its Fock matrix
            # by the transpose. Creation therefore becomes a right lower
            # Clifford shear, which modifies creation columns.
            for i,j in indices:
                first,second=locations[i],locations[j]
                for line in b:
                    line[target+first]=add(line[target+first],mul(g,line[source+second]))
                    line[target+second]=add(line[target+second],mul(g,line[source+first]))
    return b,rows,cols

class LaurentOracle:
    """Keep increasing the finite window until every promised pivot resolves."""
    def __init__(self, matrix):
        self.values=[[(a,b) for a,b,s in row] for row in matrix]
        self.shifts=[[s for a,b,s in row] for row in matrix]
        self.size=len(matrix)
        span=max((max(a.bit_length(),b.bit_length())
                  for row in matrix for a,b,s in row),default=1)
        self.precision=max(512,1<<(max(span,8*self.size)-1).bit_length())
        self.oracle=MatrixOracle(self.values,self.size,self.shifts,self.precision)

    def profile(self,*args,**kwargs):
        while True:
            try:
                return self.oracle.profile(*args,**kwargs)
            except ArithmeticError:
                if self.precision>=8192:
                    raise
                self.precision*=2
                self.oracle=MatrixOracle(self.values,self.size,self.shifts,self.precision)


def matrix_oracle(matrix):
    return LaurentOracle(matrix)

class FramePressure(ReflectionGaussianPressure):
    def __init__(self,M,gates):
        super().__init__(M);self.gates=gates;self.cuts=[]
        for side in(0,1):
            b,row,col=frame_matrix(M,side,gates)
            self.cuts.append((matrix_oracle(b),np.array(row),np.array(col)))
    def evaluate(self,parameters,gradient=False):
        answer=super().evaluate(parameters,gradient)
        self.last['base_kernel_internal_margin']=self.last.pop('internal_regularity_margin')
        return answer

def literal_tensor(M,gates):
    m=M//2;values=[norm(a,b,0) for a,b in tensor(M)['values']]
    for gate in gates:
        g=coefficient(gate);leg=gate['leg']
        for i,j in pairs(m,gate.get('matching','reverse')):
            mask=(1<<(leg*m+i))|(1<<(leg*m+j))
            for state in range(len(values)):
                occupied=state&mask
                if gate['kind']=='create' and not occupied:
                    values[state|mask]=add(values[state|mask],mul(g,values[state]))
                elif gate['kind']=='annihilate' and occupied==mask:
                    values[state^mask]=add(values[state^mask],mul(g,values[state]))
    clearing=max(0,-min(s for a,b,s in values if a|b))
    polynomial=[(a<<(s+clearing),b<<(s+clearing)) for a,b,s in values]
    for state,(a,b) in enumerate(polynomial):
        x,y,u,v=((state>>(leg*m))&((1<<m)-1) for leg in range(4))
        assert (a.bit_count()&1)+2*(b.bit_count()&1)==int(x==v and y==u)
    return polynomial,clearing

def literal_profile(M,side,heights,values):
    m=M//2;q=1<<m
    rowlegs=(0,1) if side==0 else (0,2);collegs=tuple(l for l in range(4) if l not in rowlegs)
    def weights(legs):return [sum(heights[leg*m+i] for leg,state in zip(legs,(a,b)) for i in range(m) if state>>i&1) for a,b in product(range(q),repeat=2)]
    bits=max(max(a.bit_length(),b.bit_length()) for a,b in values)
    oracle=MatrixOracle(flatten(values,m,bool(side)),1<<M,precision=max(512,1<<(bits-1).bit_length()))
    return oracle.profile(weights(rowlegs),weights(collegs),infinity=bool(side))['orders']

def predicted(M,side,heights,gates):
    b,row,col=frame_matrix(M,side,gates);h=np.asarray(heights)
    data=matrix_oracle(b).profile(np.r_[h[row],-h[row]],np.r_[h[col],-h[col]],infinity=bool(side))
    signed=np.asarray(data['orders']);assert np.array_equal(signed,-signed[::-1])
    positive=signed[M:];mean=(M*M if side==0 else 0)+sum(heights)/2
    result=[mean-positive.sum()/2]
    for order in positive:result += [a+order for a in result]
    return sorted(result)

