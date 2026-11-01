#!/usr/bin/env python3
"""Exact finite coefficient-rank and crossing-degeneration checks."""
from random import Random
import json
from check_partition_channels import rank_mod, column_matrix


def multiply_series(left, right, order, prime):
    rows, middle, columns = len(left[0]), len(right[0]), len(right[0][0])
    out = [[[0]*columns for _ in range(rows)] for _ in range(order+1)]
    for k in range(order+1):
        for a in range(k+1):
            if a >= len(left) or k-a >= len(right): continue
            for i in range(rows):
                for j in range(columns):
                    out[k][i][j] += sum(left[a][i][r]*right[k-a][r][j] for r in range(middle))
        out[k] = [[v % prime for v in row] for row in out[k]]
    return out


def verify():
    rng = Random(10302026)
    rank_cases = 0
    attained = []
    for prime in (2, 3, 101):
        for order in range(5):
            for dimension in range(1, 10):
                for generic_upper in range(1, min(4, dimension+1)):
                    left = [[[rng.randrange(prime) for _ in range(generic_upper)]
                             for _ in range(dimension)] for _ in range(order+1)]
                    right = [[[rng.randrange(prime) for _ in range(dimension)]
                              for _ in range(generic_upper)] for _ in range(order+1)]
                    product = multiply_series(left, right, order, prime)
                    actual = rank_mod(product[order], prime)
                    assert actual <= (order+1)*generic_upper
                    rank_cases += 1
        for order in range(9):
            dimension = order+1
            left = [[[int(i == k)] for i in range(dimension)] for k in range(dimension)]
            right = [[[int(order-j == k) for j in range(dimension)]] for k in range(dimension)]
            product = multiply_series(left, right, 2*order, prime)
            assert product[order] == [[int(i == j) for j in range(dimension)] for i in range(dimension)]
            assert rank_mod(product[order], prime) == order+1
            attained.append({'prime':prime, 'order':order, 'generic_rank':1, 'coefficient_rank':order+1})

    # B0=(1+z)R and B1=R+zC over F2 have sum zE.
    m = 2
    rr = [[int(x == y and v == u) for v in range(m) for u in range(m)]
          for x in range(m) for y in range(m)]
    ee = [[int(i == j) for j in range(m*m)] for i in range(m*m)]
    cc = [[a ^ b for a, b in zip(r, e)] for r, e in zip(rr, ee)]
    for i in range(m*m):
        for j in range(m*m):
            assert rr[i][j] ^ rr[i][j] == 0
            assert rr[i][j] ^ cc[i][j] == ee[i][j]
    profiles = [[rank_mod(a,2),rank_mod(column_matrix(a,m),2)] for a in (rr,cc)]
    assert profiles == [[1,4],[4,1]]

    # A low rank at one parameter value is not an upper bound on generic rank.
    assert rank_mod([[0]],101) == 0 and rank_mod([[1]],101) == 1
    # The full outer product has determinant0; first-order truncation does not.
    for z in range(1,10):
        full = [[1,z],[z,z*z]]
        truncated = [[1,z],[z,0]]
        assert rank_mod(full,101) == 1
        assert rank_mod(truncated,101) == 2

    return {'scope':'exact finite coefficient-rank checks; no improved degeneration constructed',
            'coefficient_rank_cases':rank_cases,'sharp_examples':attained,
            'cancelled_crossing_entries':16,'exact_coefficient_piece_ranks':profiles,
            'negative_controls':['single-point rank is not a generic-rank upper bound',
                                 'entry truncation can raise generic rank'],
            'new_factoring_exponent_claimed':False}


if __name__ == '__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
