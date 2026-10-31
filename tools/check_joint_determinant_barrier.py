#!/usr/bin/env python3
"""Joint determinant certificates and positive-characteristic rank controls."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations,product
import json
from math import factorial
from pathlib import Path
from random import Random

from check_multipoint_hermite_crossings import flatten,multiply,rank
from check_multipoint_rank_one_seeds import orbit_control
from search_multipoint_rank_one_seeds import KEYS,SWAP,pure_tensors,representatives

ROOT = Path(__file__).resolve().parents[1]


def determinant(matrix,prime):
    n = len(matrix)
    out = [0]*(n*(len(matrix[0][0])-1)+1)
    for p in permutations(range(n)):
        sign = (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
        value = [1]
        for i,j in enumerate(p):
            value = multiply(value,matrix[i][j],prime)
        for i,a in enumerate(value):
            out[i] = (out[i]+sign*a) % prime
    return out


def power(value,exponent,prime):
    out = [1]
    for _ in range(exponent):
        out = multiply(out,value,prime)
    return out


def exact_division(a,b,prime):
    a = list(a)
    out = [0]*max(1,len(a)-len(b)+1)
    while a and not a[-1]:
        a.pop()
    while len(a) >= len(b):
        k = len(a)-len(b)
        x = a[-1]*pow(b[-1],-1,prime) % prime
        out[k] = x
        for i,v in enumerate(b):
            a[k+i] = (a[k+i]-x*v) % prime
        while a and not a[-1]:
            a.pop()
    assert not a
    while out and not out[-1]:
        out.pop()
    return out


def polynomial_gcd(a,b,prime):
    def trim(value):
        value = list(value)
        while value and not value[-1]:
            value.pop()
        return value
    a,b = trim(a),trim(b)
    while b:
        while len(a) >= len(b):
            k = len(a)-len(b)
            scalar = a[-1]*pow(b[-1],-1,prime) % prime
            for j,value in enumerate(b):
                a[k+j] = (a[k+j]-scalar*value) % prime
            a = trim(a)
        a,b = b,a
    return a


def span_witness(vectors,target,prime):
    pivots = {}
    for i,vector in enumerate(vectors):
        row = list(vector)
        coefficients = [int(j == i) for j in range(len(vectors))]
        for p,(basis,witness) in sorted(pivots.items()):
            scalar = row[p]
            row = [(a-scalar*b) % prime for a,b in zip(row,basis)]
            coefficients = [(a-scalar*b) % prime for a,b in zip(coefficients,witness)]
        p = next((j for j,v in enumerate(row) if v),None)
        if p is not None:
            inverse = pow(row[p],-1,prime)
            pivots[p] = ([v*inverse % prime for v in row],
                         [v*inverse % prime for v in coefficients])
    row = list(target)
    answer = [0]*len(vectors)
    for p,(basis,witness) in sorted(pivots.items()):
        scalar = row[p]
        row = [(a-scalar*b) % prime for a,b in zip(row,basis)]
        answer = [(a+scalar*b) % prime for a,b in zip(answer,witness)]
    if any(row):
        return None
    assert [sum(c*v[k] for c,v in zip(answer,vectors)) % prime
            for k in range(len(target))] == list(target)
    return answer


def rectangle_multiplicity(n,t):
    numerator = factorial(n*t)
    denominator = 1
    for j in range(n):
        numerator *= factorial(j)
        denominator *= factorial(t+j)
    assert numerator % denominator == 0
    return numerator//denominator


def tableau_count(n,t):
    states = {(0,)*n:1}
    for _ in range(n*t):
        next_states = Counter()
        for shape,count in states.items():
            for i in range(n):
                if shape[i] < t and (i == 0 or shape[i] < shape[i-1]):
                    updated = list(shape)
                    updated[i] += 1
                    next_states[tuple(updated)] += count
        states = next_states
    return states[(t,)*n]


def literal_representation_controls():
    random = Random(271828)
    checks = 0
    for prime in (2,3,5):
        for t in (1,2,3):
            for _ in range(4):
                matrix = []
                for i in range(2):
                    row = []
                    for j in range(2):
                        a,b = random.randrange(prime),random.randrange(prime)
                        row.append([a,b,(int(i == j)-a-b) % prime])
                    matrix.append(row)
                det = determinant(matrix,prime)
                witness = power(det,t,prime)
                count = 2*t
                work = [[[1]]]
                for _ in range(count):
                    n = len(work)
                    work = [[multiply(work[i//2][j//2],matrix[i%2][j%2],prime)
                             for j in range(2*n)] for i in range(2*n)]
                weights = [random.randrange(prime) for _ in range(2*count+1)]
                if sum(a*b for a,b in zip(weights,witness)) % prime == 0:
                    p = next(i for i,a in enumerate(witness) if a)
                    weights[p] = (weights[p]+1) % prime
                assert sum(a*b for a,b in zip(weights,witness)) % prime
                operator = [[sum(a*b for a,b in zip(weights,value)) % prime
                             for value in row] for row in work]
                assert rank(operator,prime) >= rectangle_multiplicity(2,t)
                checks += 1
    # The nonzero determinant-functional hypothesis cannot be discarded.
    assert rank([[0]*4 for _ in range(4)],2) < rectangle_multiplicity(2,1)
    return checks


def coverage_counts(tensors,reps):
    pure_set = set(tensors)
    counts = []
    for a in reps:
        if not any(a):
            counts.append(513**2)
            continue
        line = {tuple(s*x % 3 for x in a) for s in range(3)}
        count = 0
        for b in tensors:
            if b in line:
                count += 513
            else:
                span = {tuple((s*x+t*y) % 3 for x,y in zip(a,b))
                        for s,t in product(range(3),repeat=2)}
                assert len(span) == 9
                count += len(span.intersection(pure_set))
        counts.append(count)
    assert counts == [263169,4185,4185,4185,4185,4185]
    degenerate = sum(a*b for a,b in zip(counts,orbit_control()))
    assert degenerate == 2405889
    return counts,degenerate


def asymmetric_quadratic_controls():
    from search_multipoint_rank_one_seeds import seed_from_values
    records = json.loads((ROOT/'evidence/multipoint-rank-two-seeds.json').read_text())['catalogue']
    counts = Counter()
    examples = {}
    for number,record in enumerate(records):
        seed = seed_from_values(*record['values'])
        determinants = []
        for cut in (0,1):
            coefficients = [flatten(seed,2,cut,k) for k in range(3)]
            matrix = [[tuple(coefficients[k][i][j] for k in range(3))
                       for j in range(4)] for i in range(4)]
            determinants.append(determinant(matrix,3))
        f,g = determinants
        if f == g:
            label = 'equal determinants'
        else:
            assert rank(list(seed.values()),3) == 3
            degree_f = max(i for i,value in enumerate(f) if value)
            degree_g = max(i for i,value in enumerate(g) if value)
            common_degree = len(polynomial_gcd(f,g,3))-1+min(8-degree_f,8-degree_g)
            assert common_degree in (4,5)
            # The degree-four factor at 0 and 2 is common in every case.
            h = power([0,1,1],2,3)
            assert len(exact_division(f,h,3)) <= 5
            assert len(exact_division(g,h,3)) <= 5
            label = 'common determinant degree '+str(common_degree)
        counts[label] += 1
        examples.setdefault(label,number)
    assert counts == {'equal determinants':126,'common determinant degree 4':47,
                      'common determinant degree 5':148}
    return dict(sorted(counts.items())),examples


def check():
    tableau_checks = 0
    for n in range(2,6):
        for t in range(1,7):
            h = rectangle_multiplicity(n,t)
            assert h == tableau_count(n,t)
            assert h*(n*t+1)**n*(t+n)**(n*(n-1)//2) >= n**(n*t)
            tableau_checks += 1
    literal = literal_representation_controls()
    tensors,reps = pure_tensors(),representatives()
    counts,degenerate = coverage_counts(tensors,reps)
    h = power([0,1,1],3,3)
    records = json.loads((ROOT/'evidence/multipoint-pure-cubics.json').read_text())['catalogue']
    certificates = []
    dimensions = Counter()
    for number,record in enumerate(records):
        ai,bi,ci = record['indices']
        a,b,c = reps[ai],tensors[bi],tensors[ci]
        seed = {key:(x,(y-e-z)%3,(-x-y-e)%3,z)
                for key,x,y,z,e in zip(KEYS,a,b,c,SWAP)}
        determinants = []
        for cut in (0,1):
            coefficients = [flatten(seed,2,cut,k) for k in range(4)]
            matrix = [[tuple(coefficients[k][i][j] for k in range(4))
                       for j in range(4)] for i in range(4)]
            determinants.append(determinant(matrix,3))
        f,g = determinants
        qa,qb = exact_division(f,h,3),exact_division(g,h,3)
        assert len(qa) <= 4 and len(qb) <= 4
        qa += [0]*(4-len(qa))
        qb += [0]*(4-len(qb))
        entries = [seed[key] for key in KEYS]
        dimension = rank(entries,3)
        dimensions[dimension] += 1
        if dimension < 4:
            assert f == g
        certificate = {'signature_index':number,'indices':record['indices'],
                       'entry_span_dimension':dimension,'determinants':[f,g]}
        if f == g:
            certificate.update(kind='equal determinants',rho_floor='1')
        else:
            wa,wb = span_witness(entries,qa,3),span_witness(entries,qb,3)
            assert wa is not None and wb is not None
            certificate.update(kind='joint determinant quotient',rho_floor='4/5',
                               common_factor=h,quotients=[qa,qb],
                               entry_covectors=[wa,wb])
        certificates.append(certificate)
    quadratic_counts,quadratic_examples = asymmetric_quadratic_controls()
    return {'status':'pass','representation_source':'Stephen Doty, arXiv:math/0610591v3, Lemma 13 and formal-character argument before Proposition 14',
            'tableau_multiplicity_checks':tableau_checks,
            'literal_positive_characteristic_rank_checks':literal,
            'degenerate_representatives_by_first_value':counts,
            'degenerate_actual_coefficient_seeds':degenerate,
            'full_entry_span_actual_coefficient_seeds':513**3-degenerate,
            'total_actual_coefficient_seeds':513**3,
            'cubic_witness_entry_span_dimensions':dict(sorted(dimensions.items())),
            'equal_determinant_witnesses':sum(x['kind'] == 'equal determinants' for x in certificates),
            'joint_quotient_witnesses':sum(x['kind'] == 'joint determinant quotient' for x in certificates),
            'asymmetric_quadratic_determinant_classes':quadratic_counts,
            'asymmetric_quadratic_example_signature_indices':quadratic_examples,
            'certificates':certificates,
            'scope':'Ordinary all-field joint determinant barrier; finite controls check transcription and the representation-theoretic source supplies the uniform filtration.',
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'certificates'},indent=2))


if __name__ == '__main__':
    main()
