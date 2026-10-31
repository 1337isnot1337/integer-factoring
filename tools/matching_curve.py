#!/usr/bin/env python3
"""Projective matching-value curves on a binary physical alphabet."""
from __future__ import annotations

from itertools import product
from random import Random

from check_joint_determinant_barrier import power
from check_multipoint_hermite_crossings import flatten

KEYS = tuple(product(range(2),repeat=4))


def matching_curve(prime,row_points):
    random = Random(1337)
    group = [v for v in product(range(prime),repeat=4)
             if (v[0]*v[3]-v[1]*v[2]) % prime]
    def matching(cut,a,b):
        return tuple(a[2*x+(y if cut == 0 else u)]*b[2*(u if cut == 0 else y)+v] % prime
                     for x,y,u,v in KEYS)
    crossing = tuple(int(x == v and y == u) for x,y,u,v in KEYS)
    values = [matching(0,(1,0,0,1),(1,0,0,1)),crossing]
    for a in range(2,prime):
        values.append(matching(0 if a in row_points else 1,
                               random.choice(group),random.choice(group)))
    values.append(matching(1,random.choice(group),random.choice(group)))
    basis = []
    for a in range(prime):
        v = [-x % prime for x in power([-a % prime,1],prime-1,prime)]
        v[0] = (v[0]+1) % prime
        basis.append(v+[0])
    basis.append([0,prime-1]+[0]*(prime-2)+[1])
    seed = {key:[sum(values[a][i]*basis[a][k] for a in range(prime+1)) % prime
                 for k in range(prime+1)] for i,key in enumerate(KEYS)}
    return seed,values


def cut_matrices(seed):
    degree = len(next(iter(seed.values())))-1
    result = []
    for cut in (0,1):
        coefficients = [flatten(seed,2,cut,k) for k in range(degree+1)]
        result.append([[tuple(coefficients[k][i][j] for k in range(degree+1))
                        for j in range(4)] for i in range(4)])
    return result


def cofactor_matrices(seed,prime):
    from check_joint_determinant_barrier import determinant
    return [[[determinant([[m[a][b] for b in range(4) if b != j]
                          for a in range(4) if a != i],prime)
              for j in range(4)] for i in range(4)] for m in cut_matrices(seed)]
