#!/usr/bin/env python3
"""All-field obstruction to improving the first two zero orders by I+P/t."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_joint_endpoint_loops import operator_coefficients
from check_joint_energy_riccati import dense_product, supported_rows
from check_joint_projector_lattices import polynomial_product
from joint_energy_conjugations import IDENTITY, MULTIPLICATION, N, Space, binary_product, model

ROOT = Path(__file__).resolve().parents[1]


def monomials(size):
    return [()]+[(i,) for i in range(size)]+[(i, j) for i in range(size) for j in range(i, size)]


def matrix_equations(basis, first, second, constant=0):
    """Coefficients of constant+[P,first]+P second P; monomials stay formal."""
    size = len(basis)
    terms = monomials(size)
    columns = {term: i for i, term in enumerate(terms)}
    output = [int(bool(constant >> i & 1)) for i in range(N*N)]

    def add(value, column):
        while value:
            bit = value & -value
            output[bit.bit_length()-1] ^= 1 << column
            value ^= bit

    for i, a in enumerate(basis):
        add(binary_product(a, first) ^ binary_product(first, a), i+1)
        left = binary_product(a, second)
        for j in range(i, size):
            value = binary_product(left, basis[j])
            if i != j:
                value ^= binary_product(binary_product(basis[j], second), a)
            add(value, columns[i, j])
    return output


def affine_solution(rows, size):
    space = Space(rows)
    if 0 in space.pivots:
        return None
    pivots = {p: row for p, (row, _) in space.pivots.items()}
    free = [i for i in range(1, size+1) if i not in pivots]

    def solve(value):
        for p, row in sorted(pivots.items()):
            if (row & value).bit_count() % 2:
                value ^= 1 << p
        return value >> 1

    return solve(1), [solve(1 << i) for i in free]


def consequences(rows, size):
    terms = monomials(size)
    columns = {term: i for i, term in enumerate(terms)}
    linear = supported_rows(rows, list(range(size+1)), len(terms))
    squares = supported_rows(rows, [0]+[columns[i, i] for i in range(size)], len(terms))
    return linear, squares


def assembly_controls(initial, specifications):
    terms = monomials(len(initial))
    controls = 0
    for shift in range(3):
        parameters = [(i+shift) % 4 for i in range(len(initial))]
        values = [1]+parameters+[MULTIPLICATION[parameters[i]][parameters[j]]
                                 for i in range(len(initial)) for j in range(i, len(initial))]
        planes = tuple(sum(((v >> p) & 1) << i for i, v in enumerate(parameters)) for p in range(2))
        packed = Space(initial).matrix(planes)
        matrix = [sum(((plane >> i) & 1) << p for p, plane in enumerate(packed)) for i in range(N*N)]
        for first, second, constant in specifications:
            numeric = lambda x: [(x >> i) & 1 for i in range(N*N)]
            left, right = dense_product(matrix, numeric(first)), dense_product(numeric(first), matrix)
            sandwich = dense_product(dense_product(matrix, numeric(second)), matrix)
            expected = [a ^ b ^ c ^ d for a, b, c, d in zip(left, right, sandwich, numeric(constant))]
            rows = matrix_equations(initial, first, second, constant)
            for row, target in zip(rows, expected):
                value = 0
                while row:
                    bit = row & -row
                    value ^= values[bit.bit_length()-1]
                    row ^= bit
                assert value == target
                controls += 1
    assert len(terms) == len(values)
    return controls


def two_place_pencil(data, coefficients):
    generator = sum(1 << (N*r+c) for pair in ((4, 20), (7, 11), (8, 16))
                    for r in pair for c in pair)
    allowed = {N*r+c for r, c in data['positions']}
    assert all(not (generator >> i & 1) or i in allowed for i in range(N*N))
    assert binary_product(generator, generator) == 0
    changes = [binary_product(generator, a) ^ binary_product(a, generator) for a in coefficients]
    assert all(not value for i, value in enumerate(changes) if i not in (8, 12, 16))
    assert changes[8] and changes[16] and changes[12] == changes[8] ^ changes[16]
    assert all(binary_product(binary_product(generator, a), generator) == 0 for a in coefficients)
    return dict(generator=generator, first_matrix=changes[8], second_matrix=changes[16],
                commutator_factor='t^8*(t+1)^4*(K+t^4*L)', sandwich_identities=19,
                strict_polynomial_family='A+q(t)*(K+t^4*L), degree(q)<=11, q(1)=0',
                parameter_dimension=11,
                normalized_gauge='I+q(t)*N/(t^8*(t+1)^4)')


def certificate():
    data = model()
    coefficients = operator_coefficients(data)
    regular = [0]*29
    for j, value in enumerate(coefficients):
        while value:
            bit = value & -value
            position = bit.bit_length()-1
            regular[j+data['energies'][position//N]] ^= bit
            value ^= bit
    assert not any(regular[:8]) and regular[8] & 1
    assert data['blocks'][0] == [0]
    # Strict polynomiality, then no zero-cut entries below degree eight.
    specifications = ((0, IDENTITY), (0, coefficients[0]), (coefficients[0], coefficients[1]),
                      (0, regular[8]), (regular[8], regular[9]))
    basis = [1 << (N*r+c) for r, c in data['positions']]
    controls = assembly_controls(basis, [(a, b, 0) for a, b in specifications]
                                 +[(regular[9], regular[10], regular[8]),
                                   (regular[10], regular[11], regular[9])])
    assert controls == 13125
    base_stages = []
    while True:
        size = len(basis)
        rows = [row for first, second in specifications
                for row in matrix_equations(basis, first, second) if row]
        linear, squares = consequences(rows, size)
        solution = affine_solution(linear+squares, size)
        assert solution is not None
        particular, coordinates = solution
        assert particular == 0
        base_stages.append(dict(input_dimension=size, remaining_dimension=len(coordinates),
                                linear_consequences=linear, Frobenius_consequences=squares))
        if len(coordinates) == size:
            break
        space = Space(basis)
        basis = [space.matrix((v, 0))[0] for v in coordinates]
    assert [stage['remaining_dimension'] for stage in base_stages] == [41, 36, 36]
    old_terms = monomials(size)

    def decode(row):
        result = set()
        while row:
            bit = row & -row
            result.add(old_terms[bit.bit_length()-1])
            row ^= bit
        return result

    polynomials = [decode(row) for row in rows]
    variable = lambda i: {(i,)}
    # A rank-one leading matrix u v^T, and its tangent u b^T+a v^T.
    # The unchanged singleton (0,0) entry permits u_0=v_0=1, a_0=0.
    u = [{()}]+[variable(size+i) for i in range(24)]
    v = [{()}]+[variable(size+24+i) for i in range(24)]
    a = [set()]+[variable(size+48+i) for i in range(24)]
    b = [variable(size+72+i) for i in range(25)]
    leading = matrix_equations(basis, regular[9], regular[10], regular[8])
    next_coefficient = matrix_equations(basis, regular[10], regular[11], regular[9])
    for r in range(N):
        for c in range(N):
            position = N*r+c
            polynomials.append(decode(leading[position]) ^ polynomial_product(u[r], v[c]))
            polynomials.append(decode(next_coefficient[position])
                               ^ polynomial_product(u[r], b[c]) ^ polynomial_product(a[r], v[c]))
    size += 97
    stages = []
    while True:
        terms = monomials(size)
        columns = {term: i for i, term in enumerate(terms)}
        rows = [sum(1 << columns[m] for m in poly) for poly in polynomials if poly]
        linear, squares = consequences(rows, size)
        solution = affine_solution(linear+squares, size)
        stages.append(dict(input_dimension=size, equations=len(rows),
                           linear_consequences=linear, Frobenius_consequences=squares,
                           remaining_dimension=None if solution is None else len(solution[1])))
        if solution is None:
            break
        particular, coordinates = solution
        assert len(coordinates) < size
        forms = []
        for i in range(size):
            form = {()} if particular >> i & 1 else set()
            form |= {(j,) for j, value in enumerate(coordinates) if value >> i & 1}
            forms.append(form)
        substitutions = {(): {()}}
        for term in {m for poly in polynomials for m in poly}:
            if len(term) == 1:
                substitutions[term] = forms[term[0]]
            elif len(term) == 2:
                substitutions[term] = polynomial_product(forms[term[0]], forms[term[1]])
        updated = []
        for poly in polynomials:
            value = set()
            for term in poly:
                value ^= substitutions[term]
            if value:
                updated.append(value)
        polynomials, size = updated, len(coordinates)
    assert [(stage['input_dimension'], stage['remaining_dimension']) for stage in stages] == [(133, 60), (60, None)]
    names = ('check_joint_zero_prefix.py', 'check_joint_energy_riccati.py',
             'check_joint_projector_lattices.py', 'joint_energy_conjugations.py',
             'search_f4_recentered_crossings.py')
    return dict(energy_parameters=69, necessary_linear_reduction=base_stages,
                terminal_generator_basis=basis, rank_one_tangent_reduction=stages,
                all_field_contradiction=True, independent_dense_F4_entry_controls=controls,
                two_place_pencil=two_place_pencil(data, coefficients),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in names},
                scope='Every characteristic-two field and every constant energy-block P with P^2=0, for a strict polynomial gauge I+P/t.',
                conclusion='If the smallest zero Smith order stays eight, the second order is at most nine. In particular orders (8,10,...) cannot be reached by this move, even with a changed rank-one leading matrix.',
                limit='Higher poles, products of moves, other starting tensors and pressure improvements trading away the first order are not excluded.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-zero-prefix.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_sha256', 'terminal_generator_basis',
                                   'necessary_linear_reduction', 'rank_one_tangent_reduction',
                                   'two_place_pencil')}))


if __name__ == '__main__':
    main()
