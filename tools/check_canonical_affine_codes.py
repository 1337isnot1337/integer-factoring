#!/usr/bin/env python3
"""Exhaustive exact signed affine-code columns with fixed point rows.

Each projected zero is replayed against all original local coefficients.
Projection seeds affect runtime only. NumPy and a C++17 compiler are required.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from canonical_affine_code_jets import fixed_row_jets
from cauchy_packet_polynomials import add, cauchy_polynomials, multiply, scale, shift

ROOT=Path(__file__).resolve().parents[1]
SOURCES=['tools/canonical_affine_code_jets.py','tools/check_canonical_affine_codes.py',
         'tools/check_canonical_signed_affine_codes.cpp','tools/cauchy_packet_polynomials.py']

def all_spaces(length):
    for rank in range(length+1):
        for pivots in combinations(range(length),rank):
            free=[(i,j) for i,p in enumerate(pivots) for j in range(p+1,length) if j not in pivots]
            for mask in range(1<<len(free)):
                basis=[1<<p for p in pivots]
                for k,(i,j) in enumerate(free):
                    if mask>>k&1:basis[i]|=1<<j
                support=[0]
                for v in basis:support += [x^v for x in support]
                yield basis,support

def expected_counts(length):
    spaces=characters=0
    for k in range(length+1):
        numerator=denominator=1
        for i in range(k):
            numerator*=2**(length-i)-1;denominator*=2**(k-i)-1
        count=numerator//denominator
        spaces+=count;characters+=count*(1<<k)
    return spaces,characters

def verify_candidate(matrix,row,candidate,precision):
    support=[0]
    for basis in candidate['basis']:support += [x^basis for x in support]
    signs=[1-2*((i&candidate['character']).bit_count()%2) for i in range(len(support))]
    mask=(1<<precision)-1
    return all(sum(int(coefficients[row^x])*sign for x,sign in zip(support,signs))&mask==0
               for coefficients in matrix)

def native_row(binary,directory,length,row,precision,seed=987613):
    matrix=fixed_row_jets(length,row,precision)
    rng=np.random.default_rng(seed);projections=[]
    for _ in range(3):
        weights=rng.integers(0,1<<precision,size=matrix.shape[0],dtype=np.uint64)
        values=(weights@matrix)&np.uint64((1<<precision)-1)
        projections.append(values[np.arange(1<<length)^row])
    path=directory/'projection.bin'
    np.array(projections,dtype=np.uint64).tofile(path)
    result=subprocess.run([str(binary),str(length),str(precision),str(path)],capture_output=True,text=True,check=True)
    record=json.loads(result.stderr)
    assert (record['spaces'],record['characters'])==expected_counts(length)
    witnesses=[]
    for line in result.stdout.splitlines():
        candidate=json.loads(line)
        if verify_candidate(matrix,row,candidate,precision):witnesses.append(candidate)
    return {'length':length,'row':row,'precision':precision,**record,'exact_witnesses':witnesses,
            'matrix_sha256':sha256(matrix.tobytes()).hexdigest()}

def direct_controls():
    count=0
    rng=np.random.default_rng(11497)
    for length,samples in [(2,16),(4,16),(8,12)]:
        kernels=cauchy_polynomials(length)
        for _ in range(samples):
            row=int(rng.integers(1<<length));column=int(rng.integers(1<<length))
            center=-length+1+2*int(rng.integers(length))
            polynomial=[Fraction(1)]
            for x in range(length):
                value=[Fraction(0)]
                for y in range(length):
                    if column>>y&1:value=add(value,kernels[y-x])
                if not row>>x&1:value=add([Fraction(1)],scale(value,-1))
                polynomial=multiply(polynomial,value)
            exact=shift(polynomial,center,2*length)
            modulus=1<<32
            exact=[v.numerator*pow(v.denominator,-1,modulus)%modulus for v in exact]
            degree=length*(2*length-1)
            exact += [0]*(degree+1-len(exact))
            matrix=fixed_row_jets(length,row,32)
            block=(center+length-1)//2
            assert exact==list(map(int,matrix[block*(degree+1):(block+1)*(degree+1),column]))
            count+=1
    spaces=list(all_spaces(4))
    assert len(spaces)==67 and len({tuple(sorted(s)) for _,s in spaces})==67
    return {'independent_rational_local_polynomials':count,'unique_four_dimensional_subspaces':67}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true',help='Run all 2,013,762,304 eight-wire signed affine columns')
    parser.add_argument('--output',type=Path,default=Path('/tmp/canonical-affine-codes.json'))
    args=parser.parse_args()
    np.seterr(over='ignore')
    controls=direct_controls();rows=[]
    with tempfile.TemporaryDirectory(prefix='canonical-affine-') as temporary:
        directory=Path(temporary);binary=directory/'screen'
        subprocess.run(['c++','-O3','-std=c++17',str(ROOT/SOURCES[2]),'-o',str(binary)],check=True)
        positive=native_row(binary,directory,4,5,11)
        assert positive['exact_witnesses']
        for row in range(16):
            result=native_row(binary,directory,4,row,12)
            assert not result['exact_witnesses'];rows.append(result)
        eight_control=native_row(binary,directory,8,85,31)
        assert eight_control['exact_witnesses']
        if args.full:
            for row in range(256):
                result=native_row(binary,directory,8,row,32)
                assert not result['exact_witnesses'];rows.append(result)
                if row%32==31:print(json.dumps({'eight_wire_rows_completed':row+1}),flush=True)
    record={'scope':'Fixed point rows and all signed affine-code columns; no all-input decoder or factoring improvement.',
            'arithmetic':'Exact rational controls and unsigned integer residues modulo powers of two',
            'controls':controls,'full_eight_wire_search':args.full,'rows':rows,
            'positive_controls':[positive,eight_control],
            'total_signed_affine_pairs_checked':sum(x['characters'] for x in rows),
            'sources':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}}
    args.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'checks':controls,'total_pairs':record['total_signed_affine_pairs_checked'],'full':args.full}),flush=True)

if __name__=='__main__':main()
