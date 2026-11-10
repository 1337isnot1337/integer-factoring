#!/usr/bin/env python3
"""Exact abstract structure of the twelve-site R,T orbit algebra.

The generic orbit model in the R,T orbit certificate gives two explicit
8 x 8 matrices.  This builder closes them into an eight-dimensional matrix
algebra over F5, records its multiplication table and primitive central
idempotents, and verifies that it is F5^4 x dual x dual.  The coefficient
algebra is the larger F5^6 x dual x dual.  This is an abstract structure
statement; the embedding of the orbit model into the coefficient algebra and
all cut-module ranks remain open.
"""
from __future__ import annotations

import argparse
import json
import time
from hashlib import sha256
from pathlib import Path

import build_bounded_spin_twelve_coefficient_idempotents as decomposition

ROOT = Path(__file__).resolve().parents[1]
PRIME = 5
ORBIT = ROOT / 'evidence/bounded-spin-twelve-RT-orbit.json'
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-RT-orbit-algebra.json'


def flat(matrix):
    return [entry for row in matrix for entry in row]


def unflat(vector, size):
    return [vector[i * size:(i + 1) * size] for i in range(size)]


def multiply(left, right):
    size = len(left)
    return [[sum(left[i][k] * right[k][j] for k in range(size)) % PRIME
             for j in range(size)] for i in range(size)]


def solve(columns, target):
    """Solve sum_j columns[j] * x_j = target over F5."""
    width = len(columns)
    height = len(target)
    augmented = [[columns[j][i] % PRIME for j in range(width)]
                 + [target[i] % PRIME] for i in range(height)]
    row = 0
    pivots = []
    for column in range(width):
        pivot = next((i for i in range(row, height)
                      if augmented[i][column] % PRIME), None)
        if pivot is None:
            continue
        augmented[row], augmented[pivot] = augmented[pivot], augmented[row]
        inverse = pow(augmented[row][column], -1, PRIME)
        augmented[row] = [value * inverse % PRIME
                          for value in augmented[row]]
        for i in range(height):
            if i != row and augmented[i][column] % PRIME:
                factor = augmented[i][column]
                augmented[i] = [(left - factor * right) % PRIME
                                for left, right in zip(augmented[i],
                                                       augmented[row])]
        pivots.append(column)
        row += 1
    if row != width:
        raise ValueError('the proposed basis is not independent')
    solution = [0] * width
    for i, column in enumerate(pivots):
        solution[column] = augmented[i][-1] % PRIME
    return solution


def matrix_algebra(class_sum, replica):
    size = len(class_sum)
    identity = [[int(i == j) for j in range(size)] for i in range(size)]
    basis = [identity]
    frontier = [identity]
    while frontier:
        new_frontier = []
        for current in frontier:
            for generator in (class_sum, replica):
                product = multiply(current, generator)
                candidate = [flat(item) for item in basis] + [flat(product)]
                if decomposition.rank(candidate) > len(basis):
                    basis.append(product)
                    new_frontier.append(product)
        frontier = new_frontier
    columns = [[item[i][k] for i in range(size) for k in range(size)]
               for item in basis]
    table = [[solve(columns, flat(multiply(left, right)))
              for right in basis] for left in basis]
    unit = solve(columns, flat(identity))
    return basis, table, unit


def primitive_blocks(algebra, idempotents):
    blocks = []
    for index, idempotent in enumerate(idempotents):
        block_basis = algebra.block_basis(idempotent)
        dimension = len(block_basis)
        form = [[algebra.trace_left(algebra.mul(block_basis[a],
                                                block_basis[b]))
                 for b in range(dimension)] for a in range(dimension)]
        form_rank = decomposition.rank(form)
        if dimension == 1:
            assert form == [[1]], form
            kind = 'scalar'
            nilpotent = None
        elif dimension == 2:
            assert form_rank == 1, form
            kind = 'dual'
            nilpotent = None
            for a in range(PRIME):
                for b in range(PRIME):
                    if not a and not b:
                        continue
                    candidate = [(a * block_basis[0][k]
                                  + b * block_basis[1][k]) % PRIME
                                 for k in range(algebra.size)]
                    if (algebra.mul(candidate, candidate) == [0] * algebra.size
                            and decomposition.rank([block_basis[0], candidate]) == 2):
                        nilpotent = candidate
                        break
                if nilpotent is not None:
                    break
            assert nilpotent is not None
            assert algebra.mul(block_basis[0], nilpotent) == nilpotent
            assert algebra.mul(nilpotent, nilpotent) == [0] * algebra.size
        else:
            raise AssertionError(f'unexpected block dimension {dimension}')
        blocks.append({
            'index': index,
            'idempotent': idempotent,
            'regular_rank': dimension,
            'block_basis': block_basis,
            'trace_form': form,
            'trace_form_rank': form_rank,
            'kind': kind,
            'nilpotent': nilpotent,
        })
    return blocks


def compute():
    started = time.time()
    raw = ORBIT.read_text()
    orbit = json.loads(raw)
    assert orbit['passed'] is True
    assert orbit['prime'] == PRIME
    model = orbit['twelve_site']['generic_orbit_model']
    assert model['generated_algebra_dimension'] == 8
    assert orbit['twelve_site']['full_coefficient_dimension'] == 10
    assert orbit['twelve_site']['remaining_dimensions'] == 2
    class_sum = model['class_sum_action']
    replica = model['replica_action']
    size = len(class_sum)
    assert len(replica) == len(class_sum) == size

    basis, table, unit = matrix_algebra(class_sum, replica)
    assert len(basis) == 8
    algebra = decomposition.Algebra(table, unit)
    n = algebra.size
    for i in range(n):
        e_i = [1 if k == i else 0 for k in range(n)]
        assert algebra.mul(unit, e_i) == e_i
        assert algebra.mul(e_i, unit) == e_i
    assert all(table[i][j] == table[j][i]
               for i in range(n) for j in range(n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                assert (algebra.mul(algebra.mul(
                    [1 if x == i else 0 for x in range(n)],
                    [1 if x == j else 0 for x in range(n)]),
                    [1 if x == k else 0 for x in range(n)])
                    == algebra.mul([1 if x == i else 0 for x in range(n)],
                                   algebra.mul(
                                       [1 if x == j else 0 for x in range(n)],
                                       [1 if x == k else 0 for x in range(n)])))

    idempotents = decomposition.decompose(algebra)
    total = [0] * n
    for idempotent in idempotents:
        for k in range(n):
            total[k] = (total[k] + idempotent[k]) % PRIME
    assert total == unit
    blocks = primitive_blocks(algebra, idempotents)
    assert sorted(block['regular_rank'] for block in blocks) == [1] * 4 + [2] * 2
    assert sum(1 for block in blocks if block['kind'] == 'scalar') == 4
    assert sum(1 for block in blocks if block['kind'] == 'dual') == 2
    assert not any(block['kind'] == 'field' for block in blocks)
    radicals = [block['nilpotent'] for block in blocks
                if block['kind'] == 'dual']
    assert decomposition.rank(radicals) == 2
    for radical in radicals:
        assert algebra.mul(radical, radical) == [0] * n
    assert algebra.mul(radicals[0], radicals[1]) == [0] * n

    full_form = [[algebra.trace_left(algebra.mul(
        [1 if k == i else 0 for k in range(n)],
        [1 if k == j else 0 for k in range(n)])) for j in range(n)]
        for i in range(n)]
    full_rank = decomposition.rank(full_form)
    kernel = decomposition.nullspace(full_form)
    assert full_rank == 6
    assert len(kernel) == 2
    assert decomposition.rank(kernel + radicals) == 2

    result = {
        'status': 'pass',
        'prime': PRIME,
        'length': 12,
        'orbit_dimension': 8,
        'full_coefficient_dimension': 10,
        'codimension': 2,
        'structure': 'F5^4 x F5[e]/(e^2) x F5[e]/(e^2)',
        'full_coefficient_structure': 'F5^6 x F5[e]/(e^2) x F5[e]/(e^2)',
        'regular_rank_partition': [block['regular_rank'] for block in blocks],
        'primitive_central_idempotents': blocks,
        'trace_form': {
            'rank': full_rank,
            'radical_dimension': len(kernel),
            'radical_basis': radicals,
            'square_zero': True,
        },
        'orbit_basis': [flat(matrix) for matrix in basis],
        'multiplication_table': table,
        'unit': unit,
        'class_sum_action': class_sum,
        'replica_action': replica,
        'joint_generalized_dimensions': model['joint_generalized_dimensions'],
        'orbit_evidence_sha256': sha256(raw.encode()).hexdigest(),
        'checks': {
            'generated_dimension': len(basis),
            'commutative': True,
            'associative': True,
            'unital': True,
            'idempotent_count': len(blocks),
            'four_scalar_blocks': True,
            'two_dual_blocks': True,
            'regular_rank_sum': sum(block['regular_rank'] for block in blocks),
            'trace_form_rank_plus_radical': full_rank + len(kernel),
        },
        'seconds': round(time.time() - started, 3),
        'scope': ('Exact abstract structure of the algebra generated by the '
                  'twelve-site R,T orbit actions: it is F5^4 x dual x dual, '
                  'dimension 8. The coefficient algebra is F5^6 x dual x dual, '
                  'dimension 10, so the orbit model has the abstract type of a '
                  'quotient by two scalar blocks. The actual embedding or '
                  'quotient map, which scalar blocks are involved, cut-module '
                  'ranks, and the two-cut pairing remain open; no score change.'),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = compute()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({
        'status': result['status'],
        'structure': result['structure'],
        'orbit_dimension': result['orbit_dimension'],
        'full_coefficient_dimension': result['full_coefficient_dimension'],
        'idempotents': len(result['primitive_central_idempotents']),
        'regular_rank_partition': result['regular_rank_partition'],
        'trace_form': result['trace_form']['rank'],
        'radical_dimension': result['trace_form']['radical_dimension'],
        'seconds': result['seconds'],
    }, indent=2))


if __name__ == '__main__':
    main()
