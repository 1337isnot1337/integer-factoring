#!/usr/bin/env python3
"""All-field classification of arbitrary-depth energy-projector lattices."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_joint_energy_riccati import dense_product, supported_rows
from check_joint_leading_tangent import realign
from joint_energy_conjugations import IDENTITY, MULTIPLICATION, N, Space, binary_product, kernel, model

ROOT = Path(__file__).resolve().parents[1]


def quadratic_system(basis, cold, c0, c1):
    size = len(basis)
    monomials = [(i,) for i in range(size)]+[(i, j) for i in range(size) for j in range(i, size)]
    index = {v: i for i, v in enumerate(monomials)}
    blocks = [[0]*(N*N) for _ in range(4)]

    def add(block, matrix, column):
        while matrix:
            bit = matrix & -matrix
            blocks[block][bit.bit_length()-1] ^= 1 << column
            matrix ^= bit

    for i, a in enumerate(basis):
        add(0, a, i)
        add(1, binary_product(cold, a), i)
        add(2, binary_product(a, c0), i)
        add(3, binary_product(a, c0) ^ binary_product(c0, a) ^ binary_product(a, c1), i)
        for j in range(i, size):
            b = basis[j]
            for block, coefficient in enumerate((None, cold, c0, c1)):
                value = binary_product(a, b) if coefficient is None else binary_product(binary_product(a, coefficient), b)
                if i != j:
                    value ^= binary_product(b, a) if coefficient is None else binary_product(binary_product(b, coefficient), a)
                add(block, value, index[i, j])
    return blocks, monomials


def unpack(value, monomials):
    result = set()
    while value:
        bit = value & -value
        result.add(monomials[bit.bit_length()-1])
        value ^= bit
    return result


def polynomial_product(first, second):
    result = set()
    for a in first:
        for b in second:
            monomial = tuple(sorted(a+b))
            if monomial in result:
                result.remove(monomial)
            else:
                result.add(monomial)
    return result


def reduced_polynomials(blocks, monomials, c0):
    base = [row for block in blocks[:3] for row in block if row]
    space = Space(base)

    def reduce(row):
        for pivot in sorted(space.pivots, reverse=True):
            if row >> pivot & 1:
                row ^= space.pivots[pivot][0]
        return row

    leading = [None]*(N*N)
    constant = realign(c0)
    for x, y, u, v in product(range(5), repeat=4):
        target, source = N*(x+5*u)+y+5*v, N*(x+5*y)+v+5*u
        leading[target] = unpack(reduce(blocks[3][source]), monomials)
        if constant >> target & 1:
            leading[target] ^= {()}
    # The singleton zero-energy state fixes the rank-one pivot over all fields.
    assert leading[0] == {()}
    equations = [unpack(row, monomials) for row in base]
    for r in range(1, N):
        for c in range(1, N):
            equations.append(leading[N*r+c] ^ polynomial_product(leading[N*r], leading[c]))
    return [p for p in equations if p]


def evaluate(polynomial, parameters):
    value = 0
    for monomial in polynomial:
        term = 1
        for i in monomial:
            term = MULTIPLICATION[term][parameters[i]]
        value ^= term
    return value


def dense_controls(initial, cold, c0, c1):
    blocks, monomials = quadratic_system(initial, cold, c0, c1)
    def numeric(matrix):
        return [(matrix >> i) & 1 for i in range(N*N)]
    cold, c0, c1 = map(numeric, (cold, c0, c1))
    def plus(*matrices):
        result = [0]*(N*N)
        for matrix in matrices:
            result = [a ^ b for a, b in zip(result, matrix)]
        return result
    for shift in range(4):
        parameters = [(i+shift) % 4 for i in range(len(initial))]
        planes = tuple(sum(1 << i for i, v in enumerate(parameters) if v >> p & 1) for p in range(2))
        matrix = Space(initial).matrix(planes)
        p = [((matrix[0] >> i) & 1) | (((matrix[1] >> i) & 1) << 1) for i in range(N*N)]
        right0, right1 = dense_product(p, c0), dense_product(p, c1)
        expected = [plus(dense_product(p, p), p),
                    plus(dense_product(cold, p), dense_product(dense_product(p, cold), p)),
                    plus(right0, dense_product(right0, p)),
                    plus(right0, dense_product(c0, p), right1, dense_product(right1, p))]
        actual = [[evaluate(unpack(row, monomials), parameters) for row in block] for block in blocks]
        assert actual == expected
    return 4*4*N*N


def certificate():
    data = model()
    cold = sum((data['operator'][r][c][0] & 1) << (N*r+c) for r in range(N) for c in range(N))
    c0, c1 = data['weighted'][-8], data['weighted'][-7]
    cold_next = sum(((data['operator'][r][c][0] >> 1) & 1) << (N*r+c) for r in range(N) for c in range(N))
    initial = [1 << (N*r+c) for r, c in data['positions']]
    cases = []
    for label in ('depth_one_with_leading_rank', 'depth_at_least_two_without_leading_rank',
                  'depth_one_infinity_only_with_leading_rank'):
        basis, stages = initial[:], []
        while True:
            size = len(basis)
            blocks, quadratic_monomials = quadratic_system(
                basis, 0 if label == 'depth_one_infinity_only_with_leading_rank' else cold, c0, c1)
            if label != 'depth_at_least_two_without_leading_rank':
                polynomials = reduced_polynomials(blocks, quadratic_monomials, c0)
            else:
                second, other_monomials = quadratic_system(basis, cold_next, c1, 0)
                assert quadratic_monomials == other_monomials
                polynomials = [unpack(row, quadratic_monomials)
                               for block in blocks[:3]+second[1:3] for row in block if row]
            monomials = set(m for p in polynomials for m in p) | {(i,) for i in range(size)} | {(i, i) for i in range(size)}
            ordered = sorted(monomials, key=lambda m: (len(m), m))
            index = {m: i for i, m in enumerate(ordered)}
            rows = [sum(1 << index[m] for m in polynomial) for polynomial in polynomials]
            linear = supported_rows(rows, [index[i,] for i in range(size)], len(ordered))
            square = supported_rows(rows, [index[i, i] for i in range(size)], len(ordered))
            remaining = kernel(linear+square, size)
            stages.append(dict(input_dimension=size, formal_monomials=len(ordered),
                               polynomial_row_rank=binary_rank(rows), linear_consequences=linear,
                               square_root_consequences=square, remaining_dimension=len(remaining)))
            if len(remaining) == size:
                assert size == 2
                break
            space = Space(basis)
            basis = [space.matrix((v, 0))[0] for v in remaining]
        cases.append(dict(label=label, stages=stages, terminal_basis=basis))
    assert [s['input_dimension'] for s in cases[0]['stages']] == [69, 13, 2]
    assert [s['input_dimension'] for s in cases[1]['stages']] == [69, 33, 10, 3, 2]
    assert [s['input_dimension'] for s in cases[2]['stages']] == [69, 14, 2]
    assert binary_rank(cases[0]['terminal_basis']+basis) == 2
    assert binary_rank(cases[1]['terminal_basis']+basis) == 2
    assert basis[0] ^ basis[1] == IDENTITY
    assert all(binary_product(a, a) == a for a in basis)
    assert binary_product(basis[0], basis[1]) == binary_product(basis[1], basis[0]) == 0
    assert all(binary_product(a, b) == binary_product(b, a) for a in basis for b in data['weighted'].values())
    supports = []
    for a in basis:
        support = [i for i in range(N) if a >> (N*i+i) & 1]
        assert a == sum(1 << (N*i+i) for i in support)
        supports.append(support)
    assert list(map(len, supports)) == [13, 12]
    return dict(cases=cases, neutral_component_supports=supports,
                terminal_idempotents=4,
                dense_F4_entry_controls=dense_controls(initial, cold, c0, c1)+dense_controls(initial, cold_next, c1, 0),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_projector_lattices.py', 'check_joint_energy_riccati.py',
                                            'joint_energy_conjugations.py', 'check_joint_leading_tangent.py',
                                            'search_f4_recentered_crossings.py')},
                scope='Every energy-block idempotent over every characteristic-two field. All positive depths preserve polynomiality and the pole bound; depth one requires leading rank one. The infinity-only depth-one theorem also covers every scalar rational f(t)=t+O(1), without a polynomiality assumption.',
                conclusion='Exactly four projectors: zero, identity, and the two neutral coordinate support-component projections. Every resulting tensor is the reference.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-projector-lattices.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in ('cases', 'source_sha256')}))
    print('All-field dimension chains:', {c['label']: [s['input_dimension'] for s in c['stages']]
                                          for c in result['cases']})


if __name__ == '__main__':
    main()
