#!/usr/bin/env python3
"""All-field reduction of square-zero scalar-polynomial energy dressings."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_joint_energy_dressings import certificate as sector_certificate
from joint_energy_conjugations import MULTIPLICATION, N, Space, binary_product, kernel, model, unit_characters

ROOT = Path(__file__).resolve().parents[1]


def supported_rows(equations, allowed, dimension):
    """Intersect a binary row space with a specified coordinate subspace.

    Put the allowed columns first (lowest pivots). Echelon rows whose highest
    pivot is allowed form a basis of the intersection. Return those columns
    in the order requested, also used for taking a pure-square Frobenius root.
    """
    selected = set(allowed)
    order = allowed+[i for i in range(dimension) if i not in selected]
    inverse = {value: i for i, value in enumerate(order)}
    rows = []
    for value in equations:
        row = 0
        while value:
            bit = value & -value
            row ^= 1 << inverse[bit.bit_length()-1]
            value ^= bit
        rows.append(row)
    space = Space(rows)
    return [row for pivot, (row, _) in space.pivots.items() if pivot < len(allowed)]


def equations(basis, degree, weighted):
    """Coefficients of P^2, P B_j P (j<d), and [P,B_0]+P B_d P.

    Here B_j means the coefficient at z^(-8+j). Linear and square monomials
    occupy DISTINCT columns; no finite-field polynomial reduction is used.
    """
    size = len(basis)
    monomials = {(i, j): size+k for k, (i, j) in enumerate(
        (i, j) for i in range(size) for j in range(i, size))}
    dimension = size+len(monomials)
    blocks = [[0]*(N*N) for _ in range(degree+2)]

    def add(block, matrix, column):
        while matrix:
            bit = matrix & -matrix
            blocks[block][bit.bit_length()-1] ^= 1 << column
            matrix ^= bit

    leading = weighted[-8]
    for i, a in enumerate(basis):
        add(degree+1, binary_product(a, leading) ^ binary_product(leading, a), i)
        for j in range(i, size):
            b = basis[j]
            value = binary_product(a, b)
            if i != j:
                value ^= binary_product(b, a)
            add(0, value, monomials[i, j])
        for offset in range(degree+1):
            coefficient = weighted.get(-8+offset, 0)
            if not coefficient:
                continue
            left = binary_product(a, coefficient)
            for j in range(i, size):
                value = binary_product(left, basis[j])
                if i != j:
                    value ^= binary_product(binary_product(basis[j], coefficient), a)
                add(offset+1, value, monomials[i, j])
    return blocks, monomials, dimension


def dense_product(a, b):
    # Independent literal F4 multiplication for the polynomial-assembly controls.
    multiply = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
    result = [0]*(N*N)
    for r in range(N):
        for k in range(N):
            if a[N*r+k]:
                for c in range(N):
                    result[N*r+c] ^= multiply[a[N*r+k]][b[N*k+c]]
    return result


def assembly_controls(initial, weighted):
    controls = 0
    for degree in (1, 2, 3, 17):
        blocks, monomials, _ = equations(initial, degree, weighted)
        for shift in (0, 1, 2):
            parameters = [(i+shift) % 4 for i in range(len(initial))]
            planes = tuple(sum(1 << i for i, v in enumerate(parameters) if v >> p & 1)
                           for p in range(2))
            matrix = Space(initial).matrix(planes)
            dense = [((matrix[0] >> i) & 1) | (((matrix[1] >> i) & 1) << 1)
                     for i in range(N*N)]
            expected = [dense_product(dense, dense)]
            for offset in range(degree+1):
                coefficient = [(weighted.get(-8+offset, 0) >> i) & 1 for i in range(N*N)]
                expected.append(dense_product(dense_product(dense, coefficient), dense))
            leading = [(weighted[-8] >> i) & 1 for i in range(N*N)]
            left, right = dense_product(dense, leading), dense_product(leading, dense)
            expected[-1] = [a ^ b ^ c for a, b, c in zip(expected[-1], left, right)]
            monomial_values = parameters[:]
            for i, j in monomials:
                monomial_values.append(MULTIPLICATION[parameters[i]][parameters[j]])
            low = sum((v & 1) << i for i, v in enumerate(monomial_values))
            high = sum((v >> 1) << i for i, v in enumerate(monomial_values))
            actual = [[((row & low).bit_count() % 2) | (((row & high).bit_count() % 2) << 1)
                       for row in block] for block in blocks]
            assert actual == expected
            controls += len(blocks)*N*N
    # alpha+alpha^2 has no valid pure-linear or pure-square consequence.
    assert not supported_rows([0b11], [0], 2)
    assert not supported_rows([0b11], [1], 2)
    # alpha^2+beta^2 has Frobenius root alpha+beta, leaving dimension one.
    assert kernel(supported_rows([0b1100], [2, 3], 4), 2) == [0b11]
    return controls


def certificate():
    data = model()
    assert min(data['weighted']) == -8 and max(data['weighted']) == 9
    sector = sector_certificate()
    characters, _ = unit_characters(data)
    radical = [data['leading'].matrix((v, 0))[0] for v in kernel(characters, 10)]
    initial = [1 << (N*r+c) for r, c in data['positions']]
    assert len(initial) == 69
    records = []
    for degree in range(1, 18):
        basis, stages = initial[:], []
        while True:
            size = len(basis)
            blocks, monomials, dimension = equations(basis, degree, data['weighted'])
            rows = [row for block in blocks for row in block if row]
            linear = supported_rows(rows, list(range(size)), dimension)
            square = supported_rows(rows, [monomials[i, i] for i in range(size)], dimension)
            remaining = kernel(linear+square, size)
            stages.append(dict(input_dimension=size, polynomial_row_rank=binary_rank(rows),
                               linear_consequences=linear, square_root_consequences=square,
                               remaining_dimension=len(remaining)))
            if len(remaining) == size:
                assert not rows and size == 4
                break
            space = Space(basis)
            basis = [space.matrix((v, 0))[0] for v in remaining]
        assert binary_rank(basis+radical) == 4
        records.append(dict(degree=degree, stages=stages, final_basis=basis))
    return dict(energy_block_dimension=69, nilpotent_sector_dimension=4,
                tested_degree_cases=records, dense_F4_entry_controls=assembly_controls(initial, data['weighted']),
                degree_at_least_18='B_(d-8)=0 forces commutation with B_-8; square zero then forces membership in the same four-dimensional sector.',
                sector_quadratic_identities=sector['quadratic_matrix_identities'],
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_energy_riccati.py', 'check_joint_energy_dressings.py',
                                            'joint_energy_conjugations.py', 'check_joint_leading_tangent.py',
                                            'search_f4_recentered_crossings.py')},
                scope='Every characteristic-two field, every rank, every nonconstant scalar polynomial f(t), for energy-block P with P^2=0 and U=I+f(t)P.',
                conclusion='No new infinity pole below -8 forces P into the four-dimensional sector. If the leading realignment also has rank at most one, P commutes with the reference and the seed is unchanged.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-energy-riccati.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'tested_degree_cases')}))
    print('Degree/dimension chains:', [(r['degree'], [s['input_dimension'] for s in r['stages']])
                                      for r in result['tested_degree_cases']])


if __name__ == '__main__':
    main()
