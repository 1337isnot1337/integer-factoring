#!/usr/bin/env python3
"""All-field endpoint monomial gauges with arbitrary nilpotent generators."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_joint_endpoint_loops import hasse, operator_coefficients
from check_joint_energy_riccati import dense_product, supported_rows
from check_joint_projector_lattices import evaluate, polynomial_product
from joint_energy_conjugations import N, Space, binary_product, kernel, model, unit_characters

ROOT = Path(__file__).resolve().parents[1]


def linear_matrix(basis):
    return [{(i,) for i, value in enumerate(basis) if value >> position & 1}
            for position in range(N*N)]


def right_product(first, second, blocks):
    result = []
    for r in range(N):
        for c in range(N):
            value = set()
            for j in blocks[c]:
                value ^= polynomial_product(first[N*r+j], second[N*j+c])
            result.append(value)
    return result


def equation_blocks(basis, jets, blocks, depth):
    p = linear_matrix(basis)
    square = right_product(p, p, blocks)
    cube = right_product(square, p, blocks)
    fourth = right_product(cube, p, blocks)
    commutators = {k: linear_matrix([binary_product(a, jets[k]) ^ binary_product(jets[k], a)
                                    for a in basis]) for k in range(1, 19)}
    result = [fourth]
    for i in range(1, min(depth, 18)+1):
        value = [set() for _ in range(N*N)]
        for power in range(4):
            k = i+power*depth
            if k >= 19:
                continue
            term = commutators[k]
            for _ in range(power):
                term = right_product(term, p, blocks)
            value = [a ^ b for a, b in zip(value, term)]
        result.append(value)
    return result


def assembly_controls(initial, jets, blocks):
    dense_jets = [[(value >> i) & 1 for i in range(N*N)] for value in jets]
    identity = [int(r == c) for r in range(N) for c in range(N)]
    def add(first, second):
        return [a ^ b for a, b in zip(first, second)]
    count = 0
    for depth in (1, 2, 4, 9, 17):
        symbolic = equation_blocks(initial, jets, blocks, depth)
        for shift in (0, 1):
            parameters = [(i+shift) % 4 for i in range(len(initial))]
            coordinates = tuple(sum(1 << i for i, v in enumerate(parameters) if v >> p & 1) for p in range(2))
            planes = Space(initial).matrix(coordinates)
            p = [((planes[0] >> i) & 1) | (((planes[1] >> i) & 1) << 1) for i in range(N*N)]
            powers = [identity, p]
            for _ in range(3):
                powers.append(dense_product(powers[-1], p))
            commutators = [add(dense_product(p, value), dense_product(value, p)) for value in dense_jets]
            expected = [powers[4]]
            for i in range(1, depth+1):
                value = [0]*(N*N)
                for power in range(4):
                    k = i+power*depth
                    if k < 19:
                        value = add(value, dense_product(commutators[k], powers[power]))
                expected.append(value)
            actual = [[evaluate(polynomial, parameters) for polynomial in block] for block in symbolic]
            assert actual == expected
            count += len(symbolic)*N*N
    return count


def certificate():
    data = model()
    assert max(map(len, data['blocks'])) == 4
    jets = hasse(operator_coefficients(data))
    blocks = {i: block for block in data['blocks'] for i in block}
    initial = [1 << (N*r+c) for r, c in data['positions']]
    _, characters = unit_characters(data)
    trivial = [data['full'].matrix((v, 0))[0] for v in kernel(characters, 4)]
    generator = sum(1 << (N*r+c) for pair in ((4, 20), (7, 11), (8, 16)) for r, c in product(pair, repeat=2))
    radical = trivial+[generator]
    records = []
    for depth in range(1, 18):
        basis, stages = initial[:], []
        while True:
            size = len(basis)
            polynomials = [p for block in equation_blocks(basis, jets, blocks, depth) for p in block if p]
            monomials = (set(m for p in polynomials for m in p)
                         | {(i,) for i in range(size)} | {(i, i) for i in range(size)}
                         | {(i,)*4 for i in range(size)})
            ordered = sorted(monomials, key=lambda m: (len(m), m))
            index = {m: i for i, m in enumerate(ordered)}
            rows = [sum(1 << index[m] for m in p) for p in polynomials]
            linear = supported_rows(rows, [index[i,] for i in range(size)], len(ordered))
            square = supported_rows(rows, [index[i, i] for i in range(size)], len(ordered))
            fourth = supported_rows(rows, [index[(i,)*4] for i in range(size)], len(ordered))
            remaining = kernel(linear+square+fourth, size)
            stages.append(dict(input_dimension=size, formal_monomials=len(ordered),
                               linear_consequences=linear, square_root_consequences=square,
                               fourth_root_consequences=fourth, remaining_dimension=len(remaining)))
            if len(remaining) == size:
                assert not rows
                break
            space = Space(basis)
            basis = [space.matrix((v, 0))[0] for v in remaining]
        expected = radical if depth <= 3 else trivial
        assert len(basis) == len(expected) == len(Space(basis+expected).basis)
        assert all(binary_product(a, b) == 0 for a, b in product(basis, repeat=2))
        records.append(dict(depth=depth, stages=stages, terminal_basis=basis))
    return dict(maximum_energy_block_size=4, nilpotence_bound=4, exact_depth_cases=records,
                dense_F4_entry_controls=assembly_controls(initial, jets, blocks),
                depth_at_least_18='Every strictness equation reduces to [P,D_i]=0; the nilpotent generator is in the full stabilizer and the tensor is unchanged.',
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_nilpotent_depths.py', 'check_joint_endpoint_loops.py',
                                            'check_joint_energy_riccati.py', 'check_joint_projector_lattices.py',
                                            'joint_energy_conjugations.py', 'search_f4_recentered_crossings.py')},
                scope='Every characteristic-two field, every nilpotent constant energy-block generator, every positive integer depth d, for U=I+P/(t+1)^d yielding a strict polynomial crossing.',
                conclusion='Depths one through three force the three-dimensional square-zero endpoint sector; every depth at least four forces the two-dimensional full stabilizer. All nontrivial tensors therefore inherit the proven rho>=3/4 floor at natural heights.',
                limit='General rational multipliers of a higher-nilpotency generator are not classified.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-nilpotent-depths.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'exact_depth_cases')}))
    print('Depth/dimension chains:', [(r['depth'], [s['input_dimension'] for s in r['stages']])
                                     for r in result['exact_depth_cases']])


if __name__ == '__main__':
    main()
