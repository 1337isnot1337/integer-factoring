#!/usr/bin/env python3
"""Exact primitive central idempotents of the twelve-site coefficient algebra."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT/'evidence/bounded-spin-twelve-coefficient-table.json'
PRIME = 5
N = 10
GENERIC = [1, 0, 0, 2, 1, 2, 2, 4, 1, 4]


def trim(poly):
    poly = [int(value) % PRIME for value in poly]
    while poly and poly[-1] == 0:
        poly.pop()
    return poly or [0]


def poly_add(left, right):
    result = [0]*max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] = (result[i]+value) % PRIME
    for i, value in enumerate(right):
        result[i] = (result[i]+value) % PRIME
    return trim(result)


def poly_sub(left, right):
    return poly_add(left, [(-value) % PRIME for value in right])


def poly_mul(left, right):
    left, right = trim(left), trim(right)
    result = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                result[i+j] = (result[i+j]+a*b) % PRIME
    return trim(result)


def poly_divmod(left, right):
    left, right = trim(left), trim(right)
    if len(left) < len(right):
        return [0], left
    quotient = [0]*(len(left)-len(right)+1)
    inverse = pow(right[-1], -1, PRIME)
    while len(left) >= len(right) and any(left):
        shift = len(left)-len(right)
        coefficient = left[-1]*inverse % PRIME
        quotient[shift] = coefficient
        for i, value in enumerate(right):
            left[i+shift] = (left[i+shift]-coefficient*value) % PRIME
        left = trim(left)
    return trim(quotient), left


def poly_egcd(left, right):
    if not any(right):
        return trim(left), [1], [0]
    quotient, remainder = poly_divmod(left, right)
    divisor, s, t = poly_egcd(right, remainder)
    return divisor, t, poly_sub(s, poly_mul(quotient, t))


def poly_inverse(value, modulus):
    divisor, s, _ = poly_egcd(value, modulus)
    assert len(divisor) == 1 and divisor[0], (divisor, value, modulus)
    return trim(poly_divmod(
        [entry*pow(divisor[0], -1, PRIME) % PRIME for entry in s],
        modulus)[1])


def matrix_mul(left, right):
    result = [[0]*N for _ in range(N)]
    for i in range(N):
        for k in range(N):
            coefficient = left[i][k]
            if coefficient:
                for j in range(N):
                    result[i][j] = (result[i][j]+coefficient*right[k][j]) % PRIME
    return result


def matrix_eval(matrix, poly):
    result = [[0]*N for _ in range(N)]
    power = [[int(i == j) for j in range(N)] for i in range(N)]
    for coefficient in poly:
        if coefficient:
            for i in range(N):
                for j in range(N):
                    result[i][j] = (result[i][j]+coefficient*power[i][j]) % PRIME
        power = matrix_mul(power, matrix)
    return result


def regular_matrix(table, vector):
    matrix = [[0]*N for _ in range(N)]
    for i, coefficient in enumerate(vector):
        if coefficient:
            for j in range(N):
                for k in range(N):
                    matrix[k][j] = (matrix[k][j]
                                    +coefficient*table[i][j][k]) % PRIME
    return matrix


def vectorized(matrix):
    return [matrix[i][j] for i in range(N) for j in range(N)]


def solve_consistent(columns, rhs):
    rows, count = len(rhs), len(columns)
    augmented = [[columns[j][i] for j in range(count)]+[rhs[i]]
                 for i in range(rows)]
    rank = 0
    pivots = []
    for column in range(count):
        pivot = next((i for i in range(rank, rows) if augmented[i][column]), None)
        if pivot is None:
            continue
        augmented[rank], augmented[pivot] = augmented[pivot], augmented[rank]
        inverse = pow(augmented[rank][column], -1, PRIME)
        augmented[rank] = [value*inverse % PRIME for value in augmented[rank]]
        for i in range(rows):
            if i != rank and augmented[i][column]:
                factor = augmented[i][column]
                augmented[i] = [(x-factor*y) % PRIME
                                for x, y in zip(augmented[i], augmented[rank])]
        pivots.append(column)
        rank += 1
        if rank == rows:
            break
    for i in range(rank, rows):
        if augmented[i][-1] and not any(augmented[i][:count]):
            return None
    solution = [0]*count
    for i, column in enumerate(pivots):
        solution[column] = augmented[i][-1]
    return solution


def minimal_polynomial(matrix):
    identity = [[int(i == j) for j in range(N)] for i in range(N)]
    columns, power = [vectorized(identity)], identity
    for degree in range(1, N+1):
        power = matrix_mul(power, matrix)
        solution = solve_consistent(columns, vectorized(power))
        if solution is not None:
            return [(-solution[i]) % PRIME for i in range(degree)]+[1]
        columns.append(vectorized(power))
    raise ValueError('no minimal polynomial of degree at most ten')


def root_multiplicities(poly):
    result = []
    for root in range(PRIME):
        factor = [(-root) % PRIME, 1]
        multiplicity = 0
        remaining = poly[:]
        while True:
            quotient, remainder = poly_divmod(remaining, factor)
            if any(remainder):
                break
            multiplicity += 1
            remaining = quotient
        if multiplicity:
            result.append((root, multiplicity))
    return result


def spectral_idempotent(poly, root, multiplicity, matrix):
    factor = [(-root) % PRIME, 1]
    power = [1]
    for _ in range(multiplicity):
        power = poly_mul(power, factor)
    cofactor, remainder = poly_divmod(poly, power)
    assert not any(remainder)
    value = poly_divmod(poly_mul(cofactor, poly_inverse(cofactor, power)),
                        poly)[1]
    evaluated = matrix_eval(matrix, value)
    return [evaluated[i][0] for i in range(N)]


def multiply(table, left, right):
    result = [0]*N
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    coefficient = a*b % PRIME
                    for k, value in enumerate(table[i][j]):
                        result[k] = (result[k]+coefficient*value) % PRIME
    return result


def rank(matrix):
    matrix = [row[:] for row in matrix]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(rank, len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, PRIME)
        matrix[rank] = [value*inverse % PRIME for value in matrix[rank]]
        for i in range(len(matrix)):
            if i != rank and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [(x-factor*y) % PRIME
                             for x, y in zip(matrix[i], matrix[rank])]
        rank += 1
    return rank


def image_basis(table, vector):
    columns = [multiply(table, vector, [int(i == j) for j in range(N)])
               for i in range(N)]
    basis, pivots = [], []
    for column in columns:
        reduced = column[:]
        for row, pivot in zip(basis, pivots):
            factor = reduced[pivot]*pow(row[pivot], -1, PRIME) % PRIME
            reduced = [(x-factor*y) % PRIME for x, y in zip(reduced, row)]
        pivot = next((i for i, value in enumerate(reduced) if value), None)
        if pivot is not None:
            inverse = pow(reduced[pivot], -1, PRIME)
            basis.append([value*inverse % PRIME for value in reduced])
            pivots.append(pivot)
    return basis


def block_idempotents(table, vector):
    basis = image_basis(table, vector)
    result = []
    for coefficients in product(range(PRIME), repeat=len(basis)):
        element = [sum(coefficient*basis[j][i]
                       for j, coefficient in enumerate(coefficients)) % PRIME
                   for i in range(N)]
        if multiply(table, element, element) == element:
            result.append(element)
    return basis, result


def analyze(table):
    matrix = regular_matrix(table, GENERIC)
    minimal = minimal_polynomial(matrix)
    assert minimal == [0, 1, 0, 4, 0, 4, 0, 1], minimal
    roots = root_multiplicities(minimal)
    assert roots == [(0, 1), (1, 2), (2, 1), (3, 1), (4, 2)], roots
    spectral = [spectral_idempotent(minimal, root, multiplicity, matrix)
                for root, multiplicity in roots]
    spectral_dimensions = [rank(regular_matrix(table, value))
                           for value in spectral]
    atoms = []
    block_sizes = []
    for value in spectral:
        basis, idempotents = block_idempotents(table, value)
        nonzero = [element for element in idempotents if any(element)]
        block_atoms = [element for element in nonzero
                       if len([other for other in nonzero
                               if multiply(table, element, other) == other]) == 1]
        atoms.extend(block_atoms)
        block_sizes.append(len(basis))
    return {
        'generic_element': GENERIC,
        'minimal_polynomial': minimal,
        'root_multiplicities': roots,
        'spectral_idempotents': spectral,
        'spectral_dimensions': spectral_dimensions,
        'idempotents': atoms,
        'block_dimensions': [rank(regular_matrix(table, value)) for value in atoms],
        'block_basis_dimensions': block_sizes,
        'minimal_polynomial_sha256': sha256(
            json.dumps(minimal).encode()).hexdigest(),
    }


def check(table_path=TABLE_PATH):
    receipt = json.loads(table_path.read_text())
    assert receipt['passed'] and receipt['prime'] == PRIME
    assert receipt['basis_indices'] == [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    table, unit = receipt['table'], receipt['unit']
    analysis = analyze(table)
    idempotents = analysis['idempotents']
    assert len(idempotents) == 8
    assert [sum(element[i] for element in idempotents) % PRIME
            for i in range(N)] == unit
    for i, left in enumerate(idempotents):
        assert multiply(table, left, left) == left
        for j, right in enumerate(idempotents):
            assert multiply(table, left, right) == (
                left if i == j else [0]*N)
        basis, block = block_idempotents(table, left)
        nonzero = [element for element in block if any(element)]
        assert nonzero == [left]
        assert len(basis) == analysis['block_dimensions'][i]
    assert sum(analysis['block_dimensions']) == N
    assert analysis['spectral_dimensions'] == [1, 4, 1, 1, 3]
    result = {
        'status': 'pass',
        'prime': PRIME,
        'length': 12,
        'basis_indices': receipt['basis_indices'],
        'table_sha256': sha256(table_path.read_bytes()).hexdigest(),
        'source_table': str(table_path.relative_to(ROOT)),
        'method': ('spectral projectors of one generic element, followed by '
                   'complete idempotent enumeration inside each generalized '
                   'eigenspace'),
        **analysis,
        'scope': ('Primitive central idempotents and regular-block dimensions '
                  'of the twelve-site coefficient algebra over F5 only. No '
                  'cut-module rank, integer lift, or factoring score change.'),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, default=TABLE_PATH)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check(args.table)
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
