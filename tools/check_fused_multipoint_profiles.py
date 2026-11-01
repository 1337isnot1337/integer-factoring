#!/usr/bin/env python3
"""Exact all-F16-place profiles of literal fused six-vertex crossings.

Every profile is obtained from all nonzero minors of the genuine support
blocks. Independent local block-Toeplitz ranks check the complete orders.
The finite catalogue is input to mixed-place searches, not a score claim.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

from fused_six_vertex_crossings import construct, power
from f16_identity_operator_profiles import MUL, INV, add, multiply, evaluate_one


def determinant(matrix):
    result = {}
    for permutation in permutations(range(len(matrix))):
        term = {0: 1}
        for i, j in enumerate(permutation):
            term = multiply(term, matrix[i][j])
            if not term:
                break
        # Signs coincide in characteristic two.
        result = add(result, term)
    return result


def local_polynomial(polynomial, point, degree):
    if point is None:
        return {degree - k: value for k, value in polynomial.items()}
    result = {}
    for k in range(max(polynomial, default=-1) + 1):
        coefficient = 0
        for j, value in polynomial.items():
            if k & ~j == 0:
                coefficient ^= MUL[value][power(point, j - k)]
        if coefficient:
            result[k] = coefficient
    return result


def valuation(polynomial, point, homogeneous_degree):
    assert polynomial
    if point is None:
        return homogeneous_degree - max(polynomial)
    if point == 0:
        return min(polynomial)
    return min(local_polynomial(polynomial, point, homogeneous_degree))


def field_rank(matrix):
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    r = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(r, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        inverse = INV[work[r][column]]
        work[r] = [MUL[inverse][value] for value in work[r]]
        for i in range(r + 1, len(work)):
            scalar = work[i][column]
            if scalar:
                work[i] = [a ^ MUL[scalar][b] for a, b in zip(work[i], work[r])]
        r += 1
        if r == len(work):
            break
    return r


def toeplitz_rank(matrix, length):
    size = len(matrix)
    result = [[0] * (size * length) for _ in range(size * length)]
    for i in range(length):
        for j in range(i + 1):
            for a in range(size):
                for b in range(size):
                    result[i * size + a][j * size + b] = matrix[a][b].get(i - j, 0)
    return field_rank(result)


def cut_blocks(matrix, alphabet, cut):
    size = alphabet * alphabet
    flattened = [[{} for _ in range(size)] for _ in range(size)]
    for x, y, u, v in product(range(alphabet), repeat=4):
        row, column = ((alphabet * x + y, alphabet * u + v) if cut == 0 else
                       (alphabet * x + u, alphabet * y + v))
        flattened[row][column] = matrix[alphabet * y + x][alphabet * u + v]
    unseen = set(range(size))
    blocks = []
    while unseen:
        rows, columns = {min(unseen)}, set()
        while True:
            new_columns = {j for i in rows for j in range(size) if flattened[i][j]}
            new_rows = {i for j in new_columns for i in range(size) if flattened[i][j]}
            if new_rows == rows and new_columns == columns:
                break
            rows, columns = new_rows, new_columns
        assert rows and len(rows) == len(columns)
        assert rows <= unseen
        unseen -= rows
        blocks.append([[flattened[i][j] for j in sorted(columns)] for i in sorted(rows)])
    assert sum(map(len, blocks)) == size
    return blocks


def complete_profiles(matrix, alphabet, degree):
    points = list(range(16)) + [None]
    profiles = {point: [[], []] for point in points}
    controls = []
    for cut in (0, 1):
        blocks = cut_blocks(matrix, alphabet, cut)
        minor_count = rank_checks = 0
        for block in blocks:
            size = len(block)
            prefixes = {point: 0 for point in points}
            orders = {point: [] for point in points}
            for k in range(1, size + 1):
                minors = []
                for rows in combinations(range(size), k):
                    for columns in combinations(range(size), k):
                        value = determinant([[block[i][j] for j in columns] for i in rows])
                        if value:
                            minors.append(value)
                assert minors
                minor_count += len(minors)
                for point in points:
                    current = min(valuation(value, point, degree * k) for value in minors)
                    orders[point].append(current - prefixes[point])
                    prefixes[point] = current
            for point in points:
                assert orders[point] == sorted(orders[point])
                local = [[local_polynomial(value, point, degree) for value in row]
                         for row in block]
                for length in range(1, max(orders[point]) + 2):
                    actual = toeplitz_rank(local, length)
                    expected = sum(max(0, length - order) for order in orders[point])
                    assert actual == expected
                    rank_checks += 1
                profiles[point][cut].extend(orders[point])
        for point in points:
            profiles[point][cut].sort()
            assert len(profiles[point][cut]) == alphabet * alphabet
            assert min(profiles[point][cut]) >= 0
        covered = sum(sum(profiles[point][cut]) for point in points)
        total_degree = alphabet * alphabet * degree
        assert covered <= total_degree
        controls.append({'cut': cut, 'block_sizes': list(map(len, blocks)),
                         'nonzero_minors': minor_count, 'toeplitz_ranks': rank_checks,
                         'rational_place_determinant_order': covered,
                         'homogeneous_determinant_degree': total_degree,
                         'uncovered_closed_point_degree': total_degree - covered})
    assert profiles[1] == [[0] * (alphabet * alphabet)] * 2
    return profiles, controls


def check_instance(m, theta):
    seed = construct(m, theta)
    assert seed['invariant'] and seed['endpoint_identity']
    alphabet, degree = m + 1, seed['degree']
    for i, row in enumerate(seed['matrix']):
        for j, value in enumerate(row):
            assert evaluate_one(value) == int(i == j)
    profiles, controls = complete_profiles(seed['matrix'], alphabet, degree)
    encoded = json.dumps(seed['matrix'], sort_keys=True, separators=(',', ':')).encode()
    return {'m': m, 'theta': theta, 'alphabet': alphabet, 'degree': degree,
            'cancelled_common_degree': seed['common_degree'],
            'literal_invariant_subspace_check': True, 'literal_endpoint_check': True,
            'endpoint_entries': alphabet ** 4, 'tensor_sha256': sha256(encoded).hexdigest(),
            'profiles': [{'point': point, 'cuts': profiles[point]} for point in profiles],
            'controls': controls}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [check_instance(m, theta) for m, theta in ((2, 2), (3, 4), (4, 4), (5, 4))]
    directory = Path(__file__).parent
    sources = ['check_fused_multipoint_profiles.py', 'fused_six_vertex_crossings.py',
               'f16_identity_operator_profiles.py', 'f4_lattice_linear.py']
    totals = {'instances': len(cases), 'complete_place_profiles': 34 * len(cases),
              'endpoint_entries': sum(case['endpoint_entries'] for case in cases),
              'nonzero_minors': sum(control['nonzero_minors'] for case in cases
                                    for control in case['controls']),
              'independent_toeplitz_ranks': sum(control['toeplitz_ranks'] for case in cases
                                               for control in case['controls'])}
    result = {'schema_version': 1,
              'claim': 'Exact finite all-F16-place profiles; no stronger factoring exponent.',
              'coefficient_field': 'F16 = F4[u]/(u^2+u+omega), omega^2=omega+1',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest()
                                for name in sources}, 'cases': cases, 'totals': totals}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(totals, sort_keys=True))


if __name__ == '__main__':
    main()
