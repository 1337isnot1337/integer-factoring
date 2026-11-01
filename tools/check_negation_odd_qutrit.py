#!/usr/bin/env python3
"""Construct the F7 global negation-odd crossing and certify its early jets."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

import numpy as np

from check_all_field_pascal_crossings import smith_orders, toeplitz_rank
from check_bare_seven_to_twelve_letters import ln_bounds

ROOT = Path(__file__).resolve().parents[1]
PRIME = 7


def digest(matrix):
    return sha256(np.asarray(matrix, dtype=np.uint8).tobytes()).hexdigest()


def invert(matrix):
    n = len(matrix)
    work = np.concatenate((matrix.copy(), np.eye(n, dtype=np.int64)), axis=1)
    for j in range(n):
        i = j+int(np.flatnonzero(work[j:, j])[0])
        work[[i,j]] = work[[j,i]]
        work[j] = work[j]*pow(int(work[j,j]), -1, PRIME) % PRIME
        indices = np.flatnonzero(work[:,j])
        indices = indices[indices != j]
        work[indices] = (work[indices]-work[indices,j,None]*work[j]) % PRIME
    result = work[:,n:]
    assert np.array_equal(matrix@result % PRIME, np.eye(n, dtype=np.int64))
    return result


def rank_certificate(matrix):
    """Select a nonsingular minor and verify a full rank factorization."""
    matrix = np.asarray(matrix, dtype=np.int64) % PRIME
    work, row_ids = matrix.copy(), list(range(len(matrix)))
    columns, rank = [], 0
    for j in range(matrix.shape[1]):
        candidates = np.flatnonzero(work[rank:,j])
        if not len(candidates):
            continue
        i = rank+int(candidates[0])
        work[[i,rank]] = work[[rank,i]]
        row_ids[i], row_ids[rank] = row_ids[rank], row_ids[i]
        work[rank,j:] = work[rank,j:]*pow(int(work[rank,j]), -1, PRIME) % PRIME
        indices = rank+1+np.flatnonzero(work[rank+1:,j])
        work[indices,j:] = (work[indices,j:]-work[indices,j,None]*work[rank,j:]) % PRIME
        columns.append(j)
        rank += 1
        if rank == len(matrix):
            break
    rows = row_ids[:rank]
    minor_inverse = invert(matrix[np.ix_(rows,columns)])
    right = minor_inverse@matrix[rows,:] % PRIME
    assert np.array_equal(matrix[:,columns]@right % PRIME, matrix)
    return dict(rank=rank, shape=list(matrix.shape), pivot_rows=rows, pivot_columns=columns,
                matrix_sha256=digest(matrix), right_factor_sha256=digest(right),
                minor_inverse_sha256=digest(minor_inverse),
                factorization_entries_checked=int(matrix.size))


def cell():
    states = list(product(range(3), repeat=2))
    gates = [((1,1),(0,0),2), ((0,0),(1,2),4), ((1,1),(1,1),2)]
    coefficients = [np.eye(9, dtype=np.int64)]
    for a,b,phase in gates:
        U = np.zeros((9,9), dtype=np.int64)
        for j,s in enumerate(states):
            out = tuple((x+y) % 3 for x,y in zip(s,a))
            exponent = (2*sum(x*y for x,y in zip(a,b))+sum(x*y for x,y in zip(s,b))) % 3
            U[states.index(out),j] = pow(2, exponent, PRIME)
        assert np.array_equal(U@U@U % PRIME, np.eye(9,dtype=np.int64))
        P = (np.eye(9,dtype=np.int64)+U+U@U)*pow(3,-1,PRIME) % PRIME
        assert np.array_equal(P@P % PRIME,P)
        constant = (1-phase)*P % PRIME
        linear = (np.eye(9,dtype=np.int64)-constant) % PRIME
        result = [np.zeros((9,9),dtype=np.int64) for _ in range(len(coefficients)+1)]
        for k,A in enumerate(coefficients):
            result[k] = (result[k]+A@constant) % PRIME
            result[k+1] = (result[k+1]+A@linear) % PRIME
        coefficients = result
    tensor = np.array(coefficients).reshape(4,3,3,3,3)
    assert np.array_equal(sum(coefficients) % PRIME,np.eye(9,dtype=np.int64))
    first = rank_certificate(tensor[0].reshape(9,9))
    last = rank_certificate(tensor[-1].transpose(0,3,2,1).reshape(9,9))
    assert first['rank'] == last['rank'] == 1
    return tensor, dict(polynomial_sha256=digest(tensor), zero_rank=first, infinity_rank=last)


def representatives(copies):
    return np.array([a for a in product(range(3), repeat=copies)
                     if any(a) and next(x for x in a if x) == 1], dtype=np.int64)


def construct(copies, source):
    reps = representatives(copies)
    q, degree = len(reps), 3*copies
    L = np.zeros((3**copies,q),dtype=np.int64)
    powers = 3**np.arange(copies-1,-1,-1,dtype=np.int64)
    for j,a in enumerate(reps):
        L[int(a@powers),j], L[int((-a % 3)@powers),j] = 1, PRIME-1
    R = L.T*pow(2,-1,PRIME) % PRIME
    assert np.array_equal(R@L % PRIME, np.eye(q,dtype=np.int64))
    labels = np.indices((q,)*4).reshape(4,-1)
    result = np.zeros((degree+1,q**4),dtype=np.int64)
    for signs in product((1,-1),repeat=4):
        current = np.ones((1,q**4),dtype=np.int64)
        for site in range(copies):
            indices = [(signs[j]*reps[labels[j],site]) % 3 for j in range(4)]
            factor = source[:,*indices]
            updated = np.zeros((len(current)+3,q**4),dtype=np.int64)
            for i in range(len(current)):
                for j in range(4):
                    updated[i+j] += current[i]*factor[j]
            current = updated % PRIME
        result = (result+int(np.prod(signs))*pow(4,-1,PRIME)*current) % PRIME
    tensor = result.reshape((degree+1,)+(q,)*4)
    assert np.array_equal(tensor.sum(axis=0) % PRIME,np.eye(q*q,dtype=np.int64).reshape((q,)*4))
    return tensor, reps, dict(copies=copies, alphabet=q, degree=degree, split_entries_checked=q*q,
                             endpoint_entries_checked=q**4, polynomial_sha256=digest(tensor),
                             lift_sha256=digest(L), projection_sha256=digest(R))


def literal_entry(source, reps, labels):
    result = [0]*(3*len(reps[0])+1)
    for signs in product((1,-1),repeat=4):
        current = [1]
        for site in range(len(reps[0])):
            indices = tuple(int(signs[j]*reps[labels[j],site] % 3) for j in range(4))
            factor = source[(slice(None),)+indices].tolist()
            updated = [0]*(len(current)+3)
            for i,a in enumerate(current):
                for j,b in enumerate(factor):
                    updated[i+j] = (updated[i+j]+a*b) % PRIME
            current = updated
        multiplier = int(np.prod(signs))*pow(4,-1,PRIME)
        result = [(a+multiplier*b) % PRIME for a,b in zip(result,current)]
    return result


def toeplitz(coefficients, length):
    n = coefficients.shape[1]
    result = np.zeros((n*length,n*length),dtype=np.int64)
    for i in range(length):
        for j in range(i+1):
            if i-j < len(coefficients):
                result[n*i:n*(i+1),n*j:n*(j+1)] = coefficients[i-j]
    return result


def run():
    source, source_record = cell()
    controls = []
    large = None
    for copies in (2,3):
        tensor,reps,record = construct(copies,source)
        q, degree = len(reps), len(tensor)-1
        randomizer = random.Random(20261101)
        tuples = list(product(range(q),repeat=4)) if copies == 2 else [tuple(randomizer.randrange(q) for _ in range(4)) for _ in range(128)]
        for labels in tuples:
            assert literal_entry(source,reps,labels) == tensor[(slice(None),)+labels].tolist()
        record['literal_polynomial_entries_checked'] = len(tuples)
        matrices = (tensor.reshape(degree+1,q*q,q*q),
                    tensor[::-1].transpose(0,1,4,3,2).reshape(degree+1,q*q,q*q))
        cut_records = []
        for coefficients in matrices:
            ranks = {str(length):rank_certificate(toeplitz(coefficients,length)) for length in (4,5)}
            if copies == 2:
                literal = [[tuple(map(int,coefficients[:,i,j])) for j in range(q*q)] for i in range(q*q)]
                orders,precision = smith_orders(literal,PRIME)
                for length in (4,5):
                    assert ranks[str(length)]['rank'] == sum(max(0,length-order) for order in orders)
                    assert ranks[str(length)]['rank'] == toeplitz_rank(literal,PRIME,length)
                ranks['independent_complete_orders'] = list(orders)
                ranks['independent_precision'] = precision
            ranks['orders_at_most_four'] = ranks['5']['rank']-ranks['4']['rank']
            cut_records.append(ranks)
        record['cuts_zero_and_normalized_infinity'] = cut_records
        if copies == 2:
            controls.append(record)
        else:
            large = record
    gate = min(cut['orders_at_most_four'] for cut in large['cuts_zero_and_normalized_infinity']) >= 38
    lower, _ = ln_bounds(F(38))
    _, upper = ln_bounds(F(13))
    assert lower*lower > 2*upper*upper
    return dict(scope='Explicit F7 odd-sector crossing and finite exact jets; the stopping bound concerns bare balanced pressure only.',
                prime=PRIME, omega=2, source_cell=source_record, controls=controls, target=large,
                bare_all_scale_exclusion=gate, retained_channels=38, retained_order_upper_bound=4,
                exact_log_38_lower=str(lower), exact_log_13_upper=str(upper),
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                    ('tools/check_negation_odd_qutrit.py','tools/check_all_field_pascal_crossings.py',
                     'tools/check_bare_seven_to_twelve_letters.py')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/negation-odd-qutrit.json').read_text())
    print(json.dumps({'scope':result['scope'],'bare_all_scale_exclusion':result['bare_all_scale_exclusion'],
        'counts':[cut['orders_at_most_four'] for cut in result['target']['cuts_zero_and_normalized_infinity']]}))


if __name__ == '__main__':
    main()
