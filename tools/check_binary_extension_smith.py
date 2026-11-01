#!/usr/bin/env python3
"""Independent finite-field Toeplitz and determinant controls for scalar restriction."""
from itertools import permutations
from pathlib import Path
from hashlib import sha256
import json,random
from binary_extension_smith import BinaryExtension,ExtensionMatrixOracle


def reference_multiply(a,b,mod):
    v=0
    for i in range(b.bit_length()):
        if b>>i&1:v^=a<<i
    while v.bit_length()>=mod.bit_length():v^=mod<<(v.bit_length()-mod.bit_length())
    return v


def rank(A,M):
    A=[row[:] for row in A];r=0
    for j in range(len(A[0]) if A else 0):
        i=next((i for i in range(r,len(A)) if A[i][j]),None)
        if i is None:continue
        A[r],A[i]=A[i],A[r];inv=next(x for x in range(1,len(M)) if M[x][A[r][j]]==1);A[r]=[M[inv][x] for x in A[r]]
        for i in range(r+1,len(A)):
            c=A[i][j]
            if c:A[i]=[x^M[c][y] for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r


def pmul(a,b,M):
    if not a or not b:return []
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]^=M[x][y]
    while out and not out[-1]:out.pop()
    return out


def determinant(A,M):
    n=len(A);out=[]
    for perm in permutations(range(n)):
        value=[1]
        for i,j in enumerate(perm):value=pmul(value,A[i][j],M)
        if len(out)<len(value):out += [0]*(len(value)-len(out))
        for i,v in enumerate(value):out[i]^=v
    while out and not out[-1]:out.pop()
    return out


def main():
    rng=random.Random(202611010718);records=[];toeplitz=0;field_controls=0
    for mod in (3,7,11,19):
        f=BinaryExtension(mod);M=[[reference_multiply(a,b,mod) for b in range(f.size)] for a in range(f.size)]
        for a in range(f.size):
            for b in range(f.size):assert f.multiply(a,b)==M[a][b];field_controls+=1
        for n in range(1,5):
            for trial in range(4):
                A=[[[0]*7 for _ in range(n)] for _ in range(n)]
                for i in range(n):A[i][i][rng.randrange(4)]=1
                offset=trial%3
                for i in range(n):
                    for j in range(n):
                        for k in range(3):
                            c=rng.randrange(f.size);A[i][j][k+offset]^=c;A[i][j][k+offset+1]^=c
                        while A[i][j] and not A[i][j][-1]:A[i][j].pop()
                assert all((lambda p: __import__('functools').reduce(int.__xor__,p,0))(A[i][j])==int(i==j) for i in range(n) for j in range(n))
                D=max(len(p)-1 for row in A for p in row);o=ExtensionMatrixOracle(A,mod,n);zero=o.profile()['orders'];inf=[D+v for v in o.profile(True)['orders']]
                det=determinant(A,M);assert det;assert sum(zero)==next(i for i,v in enumerate(det) if v);assert sum(inf)==n*D-(len(det)-1)
                for values,reverse in ((zero,False),(inf,True)):
                    B=[[(p+[0]*(D+1-len(p)))[::-1] if reverse else p for p in row] for row in A]
                    for h in range(1,6):
                        T=[[0]*(n*h) for _ in range(n*h)]
                        for a in range(h):
                            for b in range(a+1):
                                for i in range(n):
                                    for j in range(n):
                                        if a-b<len(B[i][j]):T[a*n+i][b*n+j]=B[i][j][a-b]
                        assert rank(T,M)==sum(max(0,h-v) for v in values);toeplitz+=1
                records.append({'modulus':mod,'n':n,'trial':trial,'zero_orders':zero,'regular_infinity_orders':inf,'degree_bound':D})
    rank_deficient=[]
    for mod in (3,7,11,19):
        A=[[[] for _ in range(4)] for _ in range(4)];A[0][0]=[0,0,1];A[1][1]=[0]*5+[1]
        oracle=ExtensionMatrixOracle(A,mod,2)
        assert oracle.profile()['orders']==[2,5]
        assert oracle.profile(True)['orders']==[-5,-2]
        rank_deficient.append({'modulus':mod,'matrix_order':4,'generic_rank':2})
    rejected=0
    for mod in (0,1,5,15):
        try:BinaryExtension(mod)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid/reducible modulus accepted.')
    for action in (lambda:BinaryExtension(11).multiply(-1,1),lambda:BinaryExtension(11).power(1,-1),lambda:ExtensionMatrixOracle([[[8]]],11,1),lambda:ExtensionMatrixOracle([[[1]]],11,0)):
        try:action()
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid input accepted.')
    names=[Path(__file__).name,'binary_extension_smith.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp']
    report={'status':'passed','field_product_controls':field_controls,'matrix_cases':len(records),'rank_deficient_controls':rank_deficient,'toeplitz_rank_controls':toeplitz,'determinant_controls':2*len(records),'invalid_controls_rejected':rejected,'records':records,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Exact unweighted local orders at zero/infinity over supplied binary extensions; generic rank is charged as an external promise.'}
    Path('evidence/binary-extension-smith-controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256')}))

if __name__=='__main__':main()
