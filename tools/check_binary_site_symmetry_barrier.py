#!/usr/bin/env python3
"""Integer identities and finite controls for the binary site-symmetry barrier."""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import permutations
import json
from math import comb, factorial, prod
from pathlib import Path
from random import Random

from check_multipoint_hermite_crossings import rank

ROOT = Path(__file__).resolve().parents[1]
MONOMIALS = ((0,6,9,15),(0,6,11,13),(2,4,9,15),(2,4,11,13),
             (1,7,8,14),(1,7,10,12),(3,5,8,14),(3,5,10,12))
SIGNS = (1,-1,-1,1,-1,1,1,-1)


def add(*polynomials):
    result = defaultdict(int)
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] += coefficient
    return {m:c for m,c in result.items() if c}


def scale(polynomial, scalar):
    return {m:scalar*c for m,c in polynomial.items() if scalar*c}


def multiply(a, b):
    result = defaultdict(int)
    for m,c in a.items():
        for n,d in b.items():
            result[tuple(sorted(m+n))] += c*d
    return {m:c for m,c in result.items() if c}


def determinant(matrix):
    n = len(matrix)
    result = {}
    for permutation in permutations(range(n)):
        term = {():1}
        for i,j in enumerate(permutation):
            term = multiply(term, matrix[i][j])
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(n) for j in range(i+1,n))
        result = add(result,scale(term,(-1)**inversions))
    return result


def universal_matrices():
    variables = [{(i,):1} for i in range(16)]
    return ([[variables[4*i+j] for j in range(4)] for i in range(4)],
            [[variables[8*(i//2)+4*(j//2)+2*(i%2)+j%2]
              for j in range(4)] for i in range(4)])


def adjugate(matrix):
    return [[scale(determinant([[matrix[r][s] for s in range(4) if s != i]
                                for r in range(4) if r != j]),(-1)**(i+j))
             for j in range(4)] for i in range(4)]


def trace_action(matrix, adj, left, right):
    # tr(diag(left) M diag(right) adj(M)).
    return add(*(scale(multiply(matrix[i][j],adj[j][i]),left[i]*right[j])
                 for i in range(4) for j in range(4)))


def identity_controls():
    expected = dict(zip(MONOMIALS,SIGNS))
    records = []
    for matrix in universal_matrices():
        adj = adjugate(matrix)
        det = determinant(matrix)
        formula = add(*(scale(multiply(matrix[i][j],adj[j][i]),sign)
                        for i,j,sign in ((0,1,1),(0,3,1),(1,0,-1),(1,2,-1))))
        assert formula == expected
        ordinary = trace_action(matrix,adj,(1,1,0,0),(0,1,0,1))
        assert add(ordinary,scale(det,-1)) == expected
        traceless = trace_action(matrix,adj,(1,1,-1,-1),(-1,1,-1,1))
        assert traceless == scale(expected,4)
        # In characteristic two, B and I both have trace zero. The two
        # output functionals need not be traceless. R(M)(I)=det(M) I.
        identity_entry = add(*(multiply(matrix[0][j],adj[j][0]) for j in range(4)))
        assert identity_entry == det
        assert sum((0,1,0,1)) % 2 == sum((1,1,1,1)) % 2 == 0
        records.append({'determinant_terms':len(det),'common_polynomial_terms':len(formula),
                        'odd_characteristic_multiplier':4,'characteristic_two_division':False})
    crossing = {0,6,9,15}
    assert sum(c for m,c in expected.items() if set(m) <= crossing) == 1
    return records


def tableau_multiplicity(shape):
    total = sum(shape)
    numerator = factorial(total)*prod(shape[i]-shape[j]+j-i
                                     for i in range(4) for j in range(i+1,4))
    denominator = prod(factorial(shape[i]+3-i) for i in range(4))
    assert numerator % denominator == 0
    return numerator//denominator


def h(t):
    numerator = comb(4*t,2*t)*comb(2*t,t)*(t+2)
    denominator = (2*t+1)*(2*t+2)
    assert numerator % denominator == 0
    return numerator//denominator


def dominated_shapes(t):
    for a in range(2*t+1):
        for b in range(a+1):
            if a+b > 3*t:
                continue
            for c in range(b+1):
                d = 4*t-a-b-c
                if 0 <= d <= c:
                    yield (a,b,c,d)


def multiplicity_controls():
    assert h(1) == 3 and h(2) == 56
    for t in range(1,65):
        assert h(t) == tableau_multiplicity((2*t,t,t,0))
        assert h(t)*2*(4*t+1)*(2*t+1)**2 >= 2**(6*t)
    shapes = 0
    for t in range(1,13):
        for shape in dominated_shapes(t):
            f = tableau_multiplicity(shape)
            assert f*(4*t+1)*(2*t+1)*(4*t+4)**6 >= 2**(6*t)
            shapes += 1
    digits = []
    for prime in (2,3,5,11):
        for t in (1,2,7,31,127,511):
            remainder,place,bound,loss = t,1,1,1
            while remainder:
                remainder,d = divmod(remainder,prime)
                bound *= h(place)**d
                loss *= (2*(4*place+1)*(2*place+1)**2)**d
                place *= prime
            assert bound*loss >= 2**(6*t)
            digits.append({'prime':prime,'t':t,'rank_floor_bits':bound.bit_length(),
                           'denominator_bits':loss.bit_length()})
    return {'hook_formula_t_range':[1,64],'dominated_shapes_checked':shapes,
            'base_prime_digit_checks':digits}


def orbit_records(monomial):
    return [tuple(sum(((word[j]>>(3-bit))&1)<<(3-j) for j in range(4))
                  for bit in range(4)) for word in permutations(monomial)]


def two_cut_ranks(data, prime):
    ranks = []
    for left,right in (((0,1),(2,3)),((0,2),(1,3))):
        rows = sorted({tuple(labels[i] for i in left) for labels in data})
        columns = sorted({tuple(labels[i] for i in right) for labels in data})
        matrix = []
        for r in rows:
            row = []
            for s in columns:
                entry = dict(zip(left,r))|dict(zip(right,s))
                row.append(data.get(tuple(entry[i] for i in range(4)),0))
            matrix.append(row)
        ranks.append(rank(matrix,prime))
    return ranks


def literal_controls():
    orbits = [orbit_records(monomial) for monomial in MONOMIALS]
    result = []
    for prime in (2,3,5):
        random = Random(1337)
        coefficients = ([[a>>i&1 for i in range(8)] for a in range(256)] if prime == 2
                        else [[random.randrange(prime) for _ in range(8)] for _ in range(64)])
        checked,minimum = 0,256
        for values in coefficients:
            if sum(a*b for a,b in zip(SIGNS,values)) % prime == 0:
                continue
            data = {labels:c for c,orbit in zip(values,orbits) if c for labels in orbit}
            ranks = two_cut_ranks(data,prime)
            assert min(ranks) >= h(1)
            checked += 1
            minimum = min(minimum,*ranks)
        result.append({'prime':prime,'nonzero_witness_controls':checked,
                       'minimum_observed_cut_rank':minimum,'proved_m4_floor':h(1)})
    # A nonsymmetric coordinate tensor shows the hypothesis cannot be dropped.
    point = orbit_records(MONOMIALS[0])[0]
    assert two_cut_ranks({point:1},2) == [1,1]
    return {'symmetric_controls':result,'nonsymmetric_negative_control_ranks':[1,1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = {'status':'pass','integer_cut_identities':identity_controls(),
              'multiplicities':multiplicity_controls(),'literal_ranks':literal_controls(),
              'scope':'Exact algebra and finite controls for the ordinary site-symmetry proof; not an independent proof of its representation-theoretic dependencies or a factoring improvement.',
              'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
