#!/usr/bin/env python3
"""Exact original-tensor controls for general coupled Gaussian gates."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_rectangular_defect_support import ZERO, ONE, Z, identity, mm
from check_reflection_unrestricted_gaussian import (
    add, conjugation_matrix, from_polynomial, quadratic_control,
    rational_smith, subset_prediction, transpose,
)
from reflection_dual_rail_gaussian import laurent_kernel
from reflection_gaussian_coupled_gates import (
    CoupledPressure, coupled_matrix, degree_bound, literal_tensor,
    local_order, selected_pairs,
)
from reflection_gaussian_physical_frames import literal_profile


def inverse_det(matrix):
    size=len(matrix)
    work=[list(row)+unit for row,unit in zip(matrix,identity(size))]
    determinant=ONE
    for j in range(size):
        pivot=next(i for i in range(j,size) if work[i][j])
        work[j],work[pivot]=work[pivot],work[j]
        value=work[j][j]
        determinant=determinant*value
        work[j]=[x/value for x in work[j]]
        for i in range(size):
            if i!=j and work[i][j]:
                value=work[i][j]
                work[i]=[x+value*y for x,y in zip(work[i],work[j])]
    return [row[size:] for row in work],determinant


def independent_kernel(M,side,gates):
    raw=laurent_kernel(M)
    kernel=[[sum((Z**power for power in entry),ZERO) for entry in row] for row in raw]
    for gate in gates:
        assert tuple(gate.get('kinds',('create','create')))==('create','create')
        a,b=gate['powers']
        value=Z**a+Z**b
        for i,j in selected_pairs(M,gate):
            kernel[i][j]=kernel[i][j]+value
            kernel[j][i]=kernel[j][i]+value
    m=M//2
    legs=(0,1) if side==0 else (0,2)
    rows=[leg*m+i for leg in legs for i in range(m)]
    columns=[i for i in range(2*M) if i not in rows]
    def block(first,second):
        return [[kernel[i][j] for j in second] for i in first]
    p,k,q=block(rows,rows),block(rows,columns),block(columns,columns)
    ki,determinant=inverse_det(k)
    kt=transpose(ki)
    pk,kq=mm(p,kt),mm(kt,q)
    upper=add(k,mm(pk,q))
    matrix=[upper[i]+pk[i] for i in range(M)]+[kq[i]+kt[i] for i in range(M)]
    mean=Fraction(M*M if side==0 else -M*M,2)+Fraction(local_order(determinant,bool(side)),2)
    return matrix,mean


def exact_original_fock(values,rows,columns,infinity):
    M=len(rows)
    row_masks=[sum(((state>>k)&1)<<physical for k,physical in enumerate(rows))
               for state in range(1<<M)]
    col_masks=[sum(((state>>k)&1)<<physical for k,physical in enumerate(columns))
               for state in range(1<<M)]
    return [[from_polynomial(values[r|c][0],infinity=infinity) for c in col_masks]
            for r in row_masks]


def verify():
    rng=np.random.default_rng(339101)
    cases=[{'legs':list(legs),'matching':'reverse','kinds':list(kinds),'powers':(-3,2)}
           for legs in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
           for kinds in product(('create','annihilate'),repeat=2)]
    words=[[gate] for gate in cases]
    words += [[{'legs':[0,1],'matching':'all','powers':(-2,3)},
               {'legs':[0,3],'matching':'reverse','powers':(-1,4)}],
              [{'legs':[0,1],'matching':'identity','kinds':['create','annihilate'],'powers':(-2,3)},
               {'legs':[0,2],'matching':'reverse','kinds':['annihilate','create'],'powers':(-1,2)}],
              [{'legs':[1,3],'matching':'all','kinds':['annihilate','annihilate'],'powers':(-1,2)},
               {'legs':[0,3],'matching':'reverse','kinds':['create','create'],'powers':(-2,3)}]]
    records=[]
    kernel_controls=rational_controls=conjugations=denominator_controls=mean_controls=0
    for M in (4,8):
        selected=words if M==4 else [words[0],words[1],words[2],words[3],words[-2],words[-1]]
        for index,gates in enumerate(selected):
            model=CoupledPressure(M,gates)
            values,clearing=literal_tensor(M,gates)
            bound=degree_bound(M,gates)
            degree=max(max(value).bit_length()-1 for value in values)
            assert clearing<=bound['clearing_bound'] and degree<=bound['cleared_degree_bound']
            m=M//2
            for state,(low,high) in enumerate(values):
                assert not high
                x,y,u,v=((state>>(leg*m))&((1<<m)-1) for leg in range(4))
                assert low.bit_count()%2==int(x==v and y==u)
            profiles=[]
            for side,(oracle,rows,columns) in enumerate(model.cuts):
                matrix,rr,cc,mean,factors=coupled_matrix(M,side,gates)
                quadratic_control(matrix)
                heights=rng.integers(-2*M,2*M+1,2*M)
                weights=(np.r_[heights[rows],-heights[rows]],np.r_[heights[columns],-heights[columns]])
                signed=oracle.profile(*weights,infinity=bool(side))['orders']
                assert all(value==round(value) for value in signed)
                signed=list(map(round,signed))
                expected=subset_prediction(signed,mean+Fraction(int(sum(heights)),2))
                expected=[value+(clearing if side==0 else -clearing) for value in expected]
                actual=literal_profile(M,side,heights.tolist(),values)
                assert expected==actual,(M,index,side)
                if oracle.orders[side]:
                    uncorrected=oracle.oracle.profile(*weights,infinity=bool(side))['orders']
                    assert uncorrected!=signed
                    denominator_controls+=1
                if mean!=(M*M if side==0 else 0):
                    mean_controls+=1
                record={'side':side,'heights':heights.tolist(),'twice_mean':int(2*mean),
                        'signed_bogoliubov_orders':signed,'fock_orders':actual,
                        'common_denominator':oracle.denominator,
                        'denominator_local_order':oracle.orders[side],
                        'elementary_determinant_factors':factors}
                if all(tuple(gate.get('kinds',('create','create')))==('create','create') for gate in gates):
                    independent,independent_mean=independent_kernel(M,side,gates)
                    assert independent==matrix and independent_mean==mean
                    kernel_controls+=1
                    record['independent_kernel_identity']=True
                if M==4 and index<4:
                    fock=exact_original_fock(values,rr,cc,bool(side))
                    rw=[sum(int(heights[physical]) for k,physical in enumerate(rr) if state>>k&1)
                        for state in range(1<<M)]
                    cw=[sum(int(heights[physical]) for k,physical in enumerate(cc) if state>>k&1)
                        for state in range(1<<M)]
                    weighted=[[value*Z**(rw[i]+cw[j]) for j,value in enumerate(row)]
                              for i,row in enumerate(fock)]
                    assert rational_smith(weighted)==actual
                    rational_controls+=1
                    local_matrix=matrix if side==0 else [
                        [from_polynomial(value.n,infinity=True)/from_polynomial(value.d,infinity=True)
                         for value in row] for row in matrix]
                    assert conjugation_matrix(fock,M)==local_matrix
                    conjugations+=1
                    record['independent_rational_fock']=True
                    record['full_conjugation_identity']=True
                profiles.append(record)
            records.append({'M':M,'gates':gates,'clearing':clearing,'degree':degree,
                            'degree_bounds':bound,'profiles':profiles})
    assert denominator_controls and mean_controls
    return {'status':'exact original-tensor coupled-gate controls passed',
            'gate_words':len(records),'complete_profiles':2*len(records),
            'independent_kernel_identities':kernel_controls,
            'independent_rational_fock_profiles':rational_controls,
            'full_conjugation_identities':conjugations,
            'nonzero_denominator_corrections':denominator_controls,
            'nontrivial_mean_corrections':mean_controls,'records':records}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify()
    root=Path(__file__).resolve().parents[1]
    sources=['research/coupled-gaussian-gates.md',
             'tools/reflection_gaussian_coupled_gates.py',
             'tools/check_reflection_gaussian_coupled_gates.py',
             'tools/reflection_gaussian_physical_frames.py',
             'tools/reflection_unrestricted_gaussian.py',
             'tools/check_reflection_unrestricted_gaussian.py']
    result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in sources}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items()
                      if key not in ('records','source_sha256')},sort_keys=True),flush=True)


if __name__=='__main__':
    main()
