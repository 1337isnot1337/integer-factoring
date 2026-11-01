#!/usr/bin/env python3
"""Independent regular-concentration and literal-tensor controls for the full-tail core."""

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from coherent_full_tail_core import core, images, pmul


def dot(a,b):
    result=0
    for x,y in zip(a,b):
        result^=pmul(x,y)
    return result


def multiply(a,b):
    columns=list(zip(*b))
    return [[dot(row,column) for column in columns] for row in a]


def shift(value,power):
    if power>=0:
        return value<<power
    assert not value&((1<<-power)-1)
    return value>>-power


def concentration_numerators(M,k):
    """C and C^-1, both multiplied by the common unit S'."""
    N=M//k
    denominator=sum(1<<(2*k*j) for j in range(N-1))
    C=[[0]*M for _ in range(M)]
    U=[[0]*M for _ in range(M)]
    for residue in range(k):
        C[residue][residue]=U[residue][residue]=denominator
        first=residue+k
        for B in range(1,N):
            row=residue+k*B
            C[row][first]=1<<(2*k*(B-1))
            U[first][row]=denominator
            if B>=2:
                C[first][row]=C[row][row]=denominator
                for B2 in range(1,N):
                    U[row][residue+k*B2]=1<<(2*k*(B-1))
                U[row][row]^=denominator
    return C,U,denominator


def block_double(matrix):
    n=len(matrix)
    return [row+[0]*n for row in matrix]+[[0]*n+row for row in matrix]


def regular_concentration_controls():
    counts={"matrix_inverse_entries":0,"weighted_regular_entries":0,
            "transported_operator_entries":0,"canonical_pairings":0}
    for k,M in ((1,4),(1,8),(2,4),(2,8),(2,16),(2,32),(4,8),(4,16),(4,32)):
        C,U,denominator=concentration_numerators(M,k)
        assert denominator&1 and denominator.bit_count()%2==1
        squared=pmul(denominator,denominator)
        for actual in (multiply(C,U),multiply(U,C)):
            for i,j in product(range(M),repeat=2):
                assert actual[i][j]==(squared if i==j else 0)
                counts["matrix_inverse_entries"]+=1
        weighted=[]
        for matrix in (C,U):
            weighted.append([[shift(matrix[i][j],2*(j-i)) for j in range(M)] for i in range(M)])
            counts["weighted_regular_entries"]+=M*M
        fullC,fullU=block_double(C),block_double(U)
        fullV,fullVinv=map(block_double,weighted)
        expected,unused=images(M,k)
        active=list(range(2*k))+list(range(M,M+2*k))
        operators=[[],[],[],[]]
        for a in range(k):
            unit=lambda j:[int(i==j) for i in range(2*M)]
            rx=(unit(a),unit(a));ry=(unit(M+a),unit(M+a))
            cv=(unit(a),unit(a));cu=(unit(M+a),unit(M+a))
            for b in range(M):
                Fba=3<<(2*(b&~a));Fab=3<<(2*(a&~b))
                rx[0][M+b]^=Fba;ry[1][b]^=Fab
                cv[1][M+b]^=Fba;cu[0][b]^=Fab
            for side,op in enumerate((rx,ry,cv,cu)):
                left,right=(fullU,fullC) if side<2 else (fullVinv,fullV)
                creation=[dot(row,op[0]) for row in left]
                annihilation=[dot(op[1],column) for column in zip(*right)]
                for which,actual in enumerate((creation,annihilation)):
                    target=[0]*(2*M)
                    for index,value in zip(active,expected[side][a][which]):
                        target[index]=pmul(denominator,value)
                    assert actual==target
                    counts["transported_operator_entries"]+=2*M
            for side in (expected[:2],expected[2:]):
                pairs=side[0]+side[1]
                for i,(creator,_) in enumerate(pairs):
                    for j,(_,annihilator) in enumerate(pairs):
                        assert dot(creator,annihilator)==int(i==j)
                        counts["canonical_pairings"]+=1
    return counts


def gaussian_factorization_controls():
    from coherent_hot_compounds import one_particle
    checked=0
    for M in (4,8,16):
        B=one_particle(M)
        F=[[3<<(2*(a&~b)) for b in range(M)] for a in range(M)]
        G=[list(row) for row in zip(*F)]
        J=[1<<(2*a+1) for a in range(M)]
        for a,b in product(range(M),repeat=2):
            assert B[a][M+b]==(J[a] if a==b else 0)
            assert B[a][b]==pmul(J[a],G[a][b])==pmul(F[a][b],J[b])
            schur=0
            for c in range(M):
                schur^=pmul(pmul(F[a][c],J[c]),G[c][b])
            schur^=(J[a]<<(2*M)) if a==b else 0
            assert schur==(J[a] if a==b else 0)
            checked+=3
    return checked


def literal_control():
    from eight_letter_reflection_retraction import ambient_normalized_tensor
    from check_block_smith_crossings import smith_profile
    from coherent_hot_compounds import profiles
    words=[word for word in range(16) if word&3!=3]
    q=len(words)
    reverse=lambda word:sum(((word>>i)&1)<<(3-i) for i in range(4))
    selected=list(map(reverse,words))
    ambient=ambient_normalized_tensor()
    matrix=[[ambient[a|(b<<4)|(c<<8)|(d<<12)] for c in selected for d in selected]
            for a in selected for b in selected]
    endpoints=0
    for x,y,u,v in product(range(q),repeat=4):
        value=matrix[q*x+y][q*u+v]
        assert value.bit_count()%2==int(x==v and y==u)
        endpoints+=1
    literal=[]
    for parity in (0,1):
        indices=[i for i in range(q*q) if (words[i//q].bit_count()+words[i%q].bit_count())%2==parity]
        block=[[matrix[i][j] for j in indices] for i in indices]
        degree=max(value.bit_length()-1 for row in block for value in row)
        values,precision=smith_profile(block,degree)
        literal.extend(values)
    expected,blocks,degree=core(4)
    binary=[1<<word for word in words]
    old=profiles(binary,binary,4)
    assert sorted(literal)==expected==old['zero_orders_t']==old['infinity_orders_t']
    return {"literal_endpoint_entries":endpoints,"independent_python_smith_orders":len(literal),
            "original_compound_orders":2*q*q,"active_core_orders":len(expected)}


def replay(path):
    saved=json.loads(path.read_text())
    assert saved['complete'] and not saved['factoring_exponent_improved']
    for name,digest in saved['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()==digest
    checked=0
    for record in saved['records']:
        actual,blocks,degree=core(record['M'])
        assert actual==record['orders'] and blocks==record['blocks']
        assert degree==record['actual_coefficient_degree']
        M=record['M']
        assert record['physical_alphabet_dimension']==3*2**(M-2)
        assert record['spectator_orders']==[2*j+1 for j in range(4,M)]+[2*M+2*j+1 for j in range(4,M)]
        assert record['physical_tail_heights']==[2*j+1 for j in range(2,M)]
        checked+=len(actual)
    return {"intervals":len(saved['records']),"recomputed_active_orders":checked,
            "finite_pressure_values_are_numerical":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result={"status":"passed","factoring_exponent_improved":False,
            "concentration":regular_concentration_controls(),
            "one_particle_factorization_entries":gaussian_factorization_controls(),
            "literal_control":literal_control()}
    if args.evidence:
        result['complete_profile_replay']=replay(args.evidence)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
