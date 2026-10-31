#!/usr/bin/env python3
"""Independent literal controls for the odd-cyclic matrix crossing interface."""
from __future__ import annotations
from fractions import Fraction
from itertools import product
import json
import random

from odd_cyclic_crossing import (flatten, matrix_rank, order_root, phase_point,
                                profiles, rational_pressure, smith_orders,
                                symbols, tensor_entry)


def toeplitz_orders(coefficients, prime):
    q = len(coefficients[0])
    degree = len(coefficients)-1
    prior_rank = prior_step = 0
    orders = []
    for length in range(1, q*degree+2):
        matrix = [[0]*(q*length) for _ in range(q*length)]
        for i in range(length):
            for j in range(i+1):
                if i-j > degree:
                    continue
                for a in range(q):
                    matrix[i*q+a][j*q:(j+1)*q] = coefficients[i-j][a][:]
        rank = matrix_rank(matrix, prime)
        step = rank-prior_rank
        assert prior_step <= step <= q
        orders += [length-1]*(step-prior_step)
        if step == q:
            return orders
        prior_rank, prior_step = rank, step
    raise AssertionError("failed to recover full polynomial matrix rank")


def verify():
    rng = random.Random(311026)
    endpoint_entries = block_entries = diagonal_entries = rank_controls = 0
    profile_controls = floor_controls = 0
    records = []
    for q, prime in ((3,7),(5,11),(7,29),(9,19)):
        root = order_root(q, prime)
        identity = [[int(i == j) for j in range(q)] for i in range(q)]
        for x,y,u,v in product(range(q), repeat=4):
            assert tensor_entry(identity,x,y,u,v) == int(x == v and y == u)
            endpoint_entries += 1
        for _ in range(8):
            matrix = [[rng.randrange(prime) for _ in range(q)] for _ in range(q)]
            sym = [value[0] for value in symbols([matrix], prime)]
            assert matrix_rank(flatten(matrix,0),prime) == sum(bool(x) for x in sym)
            assert matrix_rank(flatten(matrix,1),prime) == q*matrix_rank(matrix,prime)
            rank_controls += 2
            # Check the diagonalization entrywise in the normalized row block.
            for a,b,x in product(range(q), repeat=3):
                left = sum(tensor_entry(matrix,x,(a-x)%q,(a-v)%q,v)
                           *pow(root,b*v % q,prime) for v in range(q)) % prime
                right = sym[a*q+b]*pow(root,b*x % q,prime) % prime
                assert left == right
                diagonal_entries += 1
            # Check the column block permutation, not merely its rank.
            half = pow(2,-1,q)
            for c,x,v in product(range(q),repeat=3):
                assert tensor_entry(matrix,x,(v-c)%q,(x-c)%q,v) == (
                    matrix[(x-c*half)%q][(v-c*half)%q])
                block_entries += 1
        if q > 5:
            continue
        controls = []
        for degree in (1,2,3):
            for trial in range(5):
                coefficients = [[[rng.randrange(prime) for _ in range(q)]
                                 for _ in range(q)] for _ in range(degree)]
                coefficients.append([[(identity[i][j]-sum(a[i][j] for a in coefficients)) % prime
                                      for j in range(q)] for i in range(q)])
                controls.append(coefficients)
        # A sparse phase-space constant and a rank-one leading matrix.
        first = phase_point(q,prime)
        last = [[int(i==0 and j==0) for j in range(q)]for i in range(q)]
        middle = [[(identity[i][j]-first[i][j]-last[i][j]) % prime
                   for j in range(q)]for i in range(q)]
        controls.append([first,middle,last])
        for coefficients in controls:
            assert smith_orders(coefficients,prime) == toeplitz_orders(coefficients,prime)
            assert smith_orders(coefficients[::-1],prime) == toeplitz_orders(coefficients[::-1],prime)
            profile_controls += 2
            profile = profiles(coefficients,prime)
            for scale in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(4)):
                for orientation in (0,1):
                    assert rational_pressure(profile,q,scale,orientation) >= q**3
                    floor_controls += 1
        records.append({'alphabet':q,'prime':prime,'polynomial_seeds':len(controls)})
    # A pivot past the entry degree must be retained by the precision rule.
    delayed = [[[0,1],[0,0]],[[1,0],[0,1]]]
    assert smith_orders(delayed,7) == [0,2] == toeplitz_orders(delayed,7)
    return {'scope':'exact crossing normal form, complete local profiles, and finite controls for an ordinary all-degree bare floor',
            'endpoint_entries':endpoint_entries,'column_block_entries':block_entries,
            'row_eigenvector_entries':diagonal_entries,'literal_rank_controls':rank_controls,
            'toeplitz_profile_comparisons':profile_controls+1,
            'exact_rational_floor_comparisons':floor_controls,'seed_records':records,
            'delayed_pivot_negative_control':{'entry_degree':1,'orders':[0,2]},
            'bare_rho_lower_bound':'3/4','factoring_score_improvement':False}


if __name__ == '__main__':
    print(json.dumps(verify(),indent=2))
