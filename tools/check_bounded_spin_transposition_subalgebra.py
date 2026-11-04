#!/usr/bin/env python3
"""The replica swap and transposition class sum generate the eight-site row algebra."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from bounded_spin_crossing import paths,path_entry,path_norm
from native_gf5_flint import multiply,rank

PRIME = 5
ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/bounded-spin-transposition-subalgebra.json'


def inverse(value: int) -> int:
    return pow(int(value) % PRIME,-1,PRIME)


def frame(length: int):
    ps = paths(length)
    words = [w for w in range(1 << length) if w.bit_count() == length // 2]
    vectors = np.array([[path_entry(p,w) % PRIME for w in words] for p in ps],
                       dtype=np.int64)
    norms = np.array([path_norm(p) % PRIME for p in ps],dtype=np.int64)
    inverse_norms = np.array([inverse(x) for x in norms],dtype=np.int64)
    return ps,words,vectors,inverse_norms


def transposition_matrices(length: int,words,vectors,inverse_norms):
    index = {w:k for k,w in enumerate(words)}
    result = []
    for i in range(length):
        for j in range(i+1,length):
            def swap(word):
                left,right = (word >> i) & 1,(word >> j) & 1
                return word if left == right else word ^ (1 << i) ^ (1 << j)
            perm = [index[swap(w)] for w in words]
            matrix = ((vectors[:,perm] @ vectors.T) * inverse_norms[None,:]) % PRIME
            result.append(matrix.astype(np.uint8))
    return result


def replica_swap(length: int) -> np.ndarray:
    size = length * length
    matrix = np.zeros((size,size),dtype=np.uint8)
    for a in range(length):
        for b in range(length):
            matrix[a*length+b,b*length+a] = 1
    return matrix


def transposition_class_sum(matrices: list[np.ndarray]) -> np.ndarray:
    q = len(matrices[0])
    result = np.zeros((q*q,q*q),dtype=np.uint8)
    for matrix in matrices:
        result = (result + np.kron(matrix,matrix)) % PRIME
    return result


def nullity(matrix: np.ndarray) -> int:
    return matrix.shape[1] - rank(matrix)


def joint_dimensions(b: np.ndarray,r: np.ndarray) -> dict[tuple[int,int],int]:
    size = len(b)
    identity = np.eye(size,dtype=np.int64)
    powers = {}
    for value in range(PRIME):
        power = identity.copy()
        for _ in range(6):
            power = multiply(power,(b.astype(np.int64) - value * identity) % PRIME)
        powers[value] = power
    result = {}
    for value,power in powers.items():
        for sign in (1,PRIME-1):
            subspace = nullity(np.vstack((power,(r.astype(np.int64) - sign*identity) % PRIME)))
            if subspace:
                result[value,sign] = subspace
    assert sum(result.values()) == size
    return result


def algebra_dimension(generators: list[np.ndarray]) -> int:
    identity = np.eye(len(generators[0]),dtype=np.uint8)
    flat,basis = [],[]

    def add(matrix):
        vector = matrix.reshape(-1)
        if not flat:
            flat.append(vector)
            basis.append(matrix)
            return True
        trial = np.array(flat+[vector],dtype=np.uint8)
        if rank(trial) > len(flat):
            flat.append(vector)
            basis.append(matrix)
            return True
        return False

    for matrix in [identity]+generators:
        add(matrix)
    frontier = [identity]+generators
    for _ in range(12):
        new = []
        for left in frontier:
            for right in generators:
                product = multiply(left,right)
                if add(product):
                    new.append(product)
        frontier = new
        if not frontier:
            break
    return len(basis)


def check_length(length: int) -> dict:
    ps,words,vectors,inverse_norms = frame(length)
    q = len(ps)
    matrices = transposition_matrices(length,words,vectors,inverse_norms)
    b = transposition_class_sum(matrices)
    r = replica_swap(q)
    assert np.array_equal(multiply(b,r),multiply(r,b))
    dimensions = joint_dimensions(b,r)
    dimension = algebra_dimension([b,r])
    return {
        'length': length,
        'paths': q,
        'transpositions': len(matrices),
        'class_sum_rank': rank(b),
        'replica_swap_rank': rank(r),
        'commute': True,
        'generated_algebra_dimension': dimension,
        'joint_generalized_dimensions': [
            [value,sign,multiplicity]
            for (value,sign),multiplicity in sorted(dimensions.items())
        ],
        'primary_dimensions': sorted(dimensions.values(),reverse=True),
    }


def main() -> None:
    eight = check_length(8)
    ten = check_length(10)
    assert eight['generated_algebra_dimension'] == 6
    assert eight['primary_dimensions'] == [70,43,35,20,1]
    assert ten['generated_algebra_dimension'] == 6
    assert ten['primary_dimensions'] == [406,335,225,155,35]
    full_ten = sorted([34,1,225,300,35,406,155],reverse=True)
    assert sum(full_ten) == 1156
    assert 34+1+300 == 335
    receipt = {
        'passed': True,
        'prime': PRIME,
        'lengths': [eight,ten],
        'eight_site_full_algebra': True,
        'ten_site_full_primary_dimensions': full_ten,
        'ten_site_coarse_335_split': [34,1,300],
        'scope': (
            'The replica swap R and the site-transposition class sum T generate '
            'the full six-dimensional eight-site row algebra and reproduce its '
            'five primary dimensions. At ten sites they generate only a proper '
            'six-dimensional subalgebra with five joint generalized eigenspaces '
            '35,155,225,335,406. The published eight-dimensional row algebra '
            'splits the coarse 335 block as 34+1+300. Hence the second-coefficient '
            'identity is not enough at ten sites. No factor bound changes; the '
            'twelve-site coarse computation remains open.'
        ),
    }
    EVIDENCE.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'passed': True,'eight': eight['primary_dimensions'],
                      'ten': ten['primary_dimensions'],
                      'ten_algebra_dimension': ten['generated_algebra_dimension']}))


if __name__ == '__main__':
    main()
