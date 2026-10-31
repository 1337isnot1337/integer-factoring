#!/usr/bin/env python3
"""Common cofactor spans, Frobenius powers, and disjoint-span controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_joint_determinant_barrier import determinant,power,polynomial_gcd,rectangle_multiplicity
from check_multipoint_hermite_crossings import flatten,rank
from matching_curve import KEYS,matching_curve,cut_matrices,cofactor_matrices

ROOT = Path(__file__).resolve().parents[1]


def invert_with_rows(matrix,prime):
    work = [list(row) for row in matrix]
    labels = list(range(len(work)))
    selected = []
    n = len(work[0])
    for k in range(n):
        hit = next(i for i in range(k,len(work)) if work[i][k])
        work[k],work[hit] = work[hit],work[k]
        labels[k],labels[hit] = labels[hit],labels[k]
        selected.append(labels[k])
        inverse = pow(work[k][k],-1,prime)
        work[k] = [x*inverse % prime for x in work[k]]
        for i in range(k+1,len(work)):
            scalar = work[i][k]
            work[i] = [(a-scalar*b) % prime for a,b in zip(work[i],work[k])]
    square = [list(matrix[i]) for i in selected]
    augmented = [row+[int(i == j) for j in range(n)] for i,row in enumerate(square)]
    for k in range(n):
        hit = next(i for i in range(k,n) if augmented[i][k])
        augmented[k],augmented[hit] = augmented[hit],augmented[k]
        inverse = pow(augmented[k][k],-1,prime)
        augmented[k] = [x*inverse % prime for x in augmented[k]]
        for i in range(n):
            if i != k:
                scalar = augmented[i][k]
                augmented[i] = [(a-scalar*b) % prime for a,b in zip(augmented[i],augmented[k])]
    return selected,[row[n:] for row in augmented]


def degree_eleven_candidate():
    seed,values = matching_curve(11,(0,2,3,4,5,6))
    cofactors = cofactor_matrices(seed,11)
    matrix = [[cofactors[0][i][j][k] for i in range(4) for j in range(4)]
              +[-cofactors[1][i][j][k] % 11 for i in range(4) for j in range(4)]
              for k in range(34)]
    rows,inverse = invert_with_rows(matrix,11)
    return {'prime':11,'degree':11,'keys':KEYS,'seed':[seed[k] for k in KEYS],
            'values':values,'row_points':[0,2,3,4,5,6],
            'independent_cofactor_rows':rows,'cofactor_minor_inverse':inverse}


def check_degree_eleven(candidate):
    rebuilt = json.loads(json.dumps(degree_eleven_candidate()))
    assert candidate == rebuilt
    seed = dict(zip(KEYS,candidate['seed']))
    for a in range(11):
        assert all(sum(c*pow(a,k,11) for k,c in enumerate(seed[key])) % 11 == candidate['values'][a][i]
                   for i,key in enumerate(KEYS))
    assert [seed[key][-1] for key in KEYS] == candidate['values'][-1]
    pairs = []
    for a,value in enumerate(candidate['values']):
        table = {key:[x] for key,x in zip(KEYS,value)}
        pair = [rank(flatten(table,2,cut),11) for cut in (0,1)]
        expected = [4,4] if a == 1 else [1,4] if a in candidate['row_points'] else [4,1]
        assert pair == expected
        pairs.append(pair)
    assert rank(list(seed.values()),11) == 12
    ds = [determinant(m,11) for m in cut_matrices(seed)]
    assert len(polynomial_gcd(*ds,11)) == 1
    assert max(i for i,c in enumerate(ds[0]) if c) == 44
    cofactors = cofactor_matrices(seed,11)
    matrix = [[cofactors[0][i][j][k] for i in range(4) for j in range(4)]
              +[-cofactors[1][i][j][k] % 11 for i in range(4) for j in range(4)]
              for k in range(34)]
    square = [matrix[i] for i in candidate['independent_cofactor_rows']]
    inverse = candidate['cofactor_minor_inverse']
    assert len(square) == len(inverse) == 32
    for a,b in ((square,inverse),(inverse,square)):
        assert [[sum(a[i][k]*b[k][j] for k in range(32)) % 11
                 for j in range(32)] for i in range(32)] == [[int(i == j) for j in range(32)] for i in range(32)]
    return {'projective_point_cut_ranks':pairs,'entry_span_dimension':12,
            'cofactor_coefficient_matrix_shape':[34,32],'certified_cofactor_span_rank':32,
            'independent_inverse_product_entries':2048,
            'common_cofactor_span_dimension':0,'homogeneous_cut_determinant_common_degree':0}


def degree_five_span_control():
    candidate = json.loads((ROOT/'evidence/degree-five-matching-candidate.json').read_text())
    seed = dict(zip(KEYS,candidate['seed']))
    weights = [[row[0] for row in candidate['affine_coordinates'][offset:offset+16]]
               for offset in (0,16)]
    cofactors = cofactor_matrices(seed,5)
    spaces = [[c[i][j] for i in range(4) for j in range(4)] for c in cofactors]
    dimensions = [rank(space,5) for space in spaces]
    combined = rank(spaces[0]+spaces[1],5)
    assert dimensions == [10,12] and combined == 16
    polynomials = [[sum(w[4*i+j]*c[i][j][k] for i in range(4) for j in range(4)) % 5
                    for k in range(16)] for c,w in zip(cofactors,weights)]
    assert polynomials[0] == polynomials[1] and sum(polynomials[0]) % 5 == 1
    p = polynomials[0]
    for exponent in (1,5,25):
        actual = power(p,exponent,5)
        expected = [0]*(15*exponent+1)
        for i,c in enumerate(p):
            expected[i*exponent] = pow(c,exponent,5)
        assert actual == expected
    comparisons = []
    for r in range(5):
        t = 5**r
        h = rectangle_multiplicity(3,t)
        comparisons.append({'frobenius_exponent':t,'word_length':3*t,
                            'rank_floor_bits':h.bit_length(),
                            'rank_product_exceeds_rho_five_sevenths':h**14 > 2**(60*t)})
    return {'common_polynomial':p,'entry_covectors':weights,
            'cofactor_span_dimensions':dimensions,'common_cofactor_span_dimension':sum(dimensions)-combined,
            'cofactor_weight_matrix_ranks':[rank([w[4*i:4*i+4] for i in range(4)],5) for w in weights],
            'frobenius_identity_exponents':[1,5,25],'exact_rank_floor_comparisons':comparisons}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate-candidate',type=Path)
    parser.add_argument('--candidate',type=Path,default=ROOT/'evidence/degree-eleven-matching-candidate.json')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if args.generate_candidate:
        args.generate_candidate.write_text(json.dumps(degree_eleven_candidate(),separators=(',',':'))+'\n')
        args.candidate = args.generate_candidate
    candidate = json.loads(args.candidate.read_text())
    result = {'status':'pass','degree_five_common_span':degree_five_span_control(),
              'degree_eleven_candidate':check_degree_eleven(candidate),
              'scope':'Common-cofactor-span controls. Both curves are excluded by the broader binary site-symmetry proof; no factoring improvement.',
              'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                                for name in ('matching_curve.py','check_frobenius_cofactor_span.py')},
              'degree_eleven_candidate_sha256':sha256(args.candidate.read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
