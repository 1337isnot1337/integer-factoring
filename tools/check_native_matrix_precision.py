#!/usr/bin/env python3
"""Exact wide-window controls and missing-pivot guards for the F4 engine."""
import argparse
import hashlib
from itertools import combinations,permutations
import json
from pathlib import Path
from check_rectangular_f4_repair import add,multiply
from native_matrix_smith_f4 import MatrixOracle

ZERO=(0,0);ONE=(1,0)
def matmul(a,b):
    out=[[ZERO for _ in b[0]] for _ in a]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(b)):out[i][j]=add(out[i][j],multiply(a[i][k],b[k][j]))
    return out
def determinant(matrix):
    total=ZERO
    for permutation in permutations(range(len(matrix))):
        term=ONE
        for i,j in enumerate(permutation):term=multiply(term,matrix[i][j])
        total=add(total,term)
    return total
def reference(matrix,infinity):
    minima=[0];count=0
    for k in range(1,len(matrix)+1):
        candidates=[]
        for rows in combinations(range(len(matrix)),k):
            for columns in combinations(range(len(matrix)),k):
                value=determinant([[matrix[i][j] for j in columns] for i in rows]);count+=1
                if value!=ZERO:
                    bits=value[0]|value[1];candidates.append(-(bits.bit_length()-1) if infinity else (bits&-bits).bit_length()-1)
        minima.append(min(candidates))
    return [minima[i+1]-minima[i] for i in range(len(matrix))],count
def check():
    lower=[[ONE,ZERO,ZERO],[(1,1<<600),ONE,ZERO],[(1<<37,1),(1^(1<<500),0),ONE]]
    upper=[[ONE,(0,1<<400),(1<<250,1)],[ZERO,ONE,(1,1<<300)],[ZERO,ZERO,ONE]]
    diagonal=[[ONE,ZERO,ZERO],[ZERO,(1<<700,0),ZERO],[ZERO,ZERO,(1<<1500,0)]]
    matrix=matmul(matmul(lower,diagonal),upper);cases=[];minors=0
    for infinity in (False,True):
        expected,count=reference(matrix,infinity);minors+=count
        for precision in (4096,8192):
            actual=MatrixOracle(matrix,3,precision=precision).profile(infinity=infinity)['orders']
            assert actual==expected
            cases.append({'infinity':infinity,'precision':precision,'orders':actual})
    shifted=[[ONE,ZERO],[ZERO,ONE]];shifts=[[0,0],[0,600]]
    try:MatrixOracle(shifted,2,shifts).profile()
    except ArithmeticError:guard=True
    else:raise AssertionError('An unresolved promised pivot must not be omitted.')
    assert MatrixOracle(shifted,2,shifts,precision=1024).profile()['orders']==[0.,600.]
    names=['native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp',Path(__file__).name]
    return {'status':'pass','exact_minors':minors,'wide_cases':cases,'narrow_window_missing_pivot_rejected':guard,'expanded_window_orders':[0,600],'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args();result=check();text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text)
if __name__=='__main__':main()
