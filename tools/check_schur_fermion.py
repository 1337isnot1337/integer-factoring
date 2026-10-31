#!/usr/bin/env python3
"""Independent exact audit of the charged-fermion Schur query algorithm.

The reference uses Jacobi--Trudi power sums and the original bosonic recursion.
It is adapted, with attribution, from PR #30 at 28105165 (Apache-2.0).
Neither the reference nor its characteristic-zero basis change is used by the
new query algorithm in schur_fermion.py.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from math import factorial
from pathlib import Path
import json

import check_partition_channels as source
import schur_fermion as new


def clean(poly):
    return tuple(sorted((part, c) for part, c in poly.items() if c))


def multiply(left, right):
    answer = defaultdict(Fraction)
    for p, a in left:
        for q, b in right:
            answer[tuple(sorted(p + q))] += a * b
    return clean(answer)


@lru_cache(None)
def complete(weight):
    return tuple((p, Fraction(1, source.partition_denominator(p)))
                 for p in source.partitions(weight))


@lru_cache(None)
def schur(part):
    lam = part[::-1]
    length = len(lam)
    answer = defaultdict(Fraction)
    for sigma in permutations(range(length)):
        inversions = sum(sigma[i] > sigma[j]
                         for i in range(length) for j in range(i + 1, length))
        term = (((), Fraction((-1) ** inversions)),)
        for i, j in enumerate(sigma):
            term = multiply(term, complete(lam[i] - i + j))
            if not term:
                break
        for p, c in term:
            answer[p] += c
    return clean(answer)


def coefficient(poly, schur_label):
    coefficients = dict(poly)
    return sum(c * source.partition_denominator(p) * coefficients.get(p, 0)
               for p, c in schur(schur_label))


def mode(parameter, input_, output_weight):
    answer = {}
    for out in source.partitions(output_weight):
        value = sum(a * b * source.operator(p, out, q)
                    for p, a in parameter for q, b in input_)
        if value:
            answer[out] = value
    return clean(answer)


def matrix_multiply(left, right):
    return [[sum(a * b for a, b in zip(row, column))
             for column in zip(*right)] for row in left]


def tensor_square(matrix):
    d = len(matrix)
    return [[matrix[a][x] * matrix[b][y] for x in range(d) for y in range(d)]
            for a in range(d) for b in range(d)]


def source_schur_piece(h, shift):
    basis = source.partitions(h)
    change = [[dict(schur(label)).get(part, Fraction(0)) for part in basis]
              for label in basis]
    characters = [[a * source.partition_denominator(part)
                   for a, part in zip(row, basis)] for row in change]
    old = source.piece_matrix(h, shift)
    transformed = matrix_multiply(matrix_multiply(tensor_square(change), old),
                                  list(map(list, zip(*tensor_square(characters)))))
    answer = [[entry / factorial(h) for entry in row] for row in transformed]
    assert all(entry.denominator == 1 for row in answer for entry in row)
    return [[int(entry) for entry in row] for row in answer]


def verify(max_weight=4, max_crossing=4):
    integer_modes = field_modes = binary_coefficients = rational_values = 0
    maximum_integer_mode_bits = 0
    for h in range(1, max_weight + 1):
        for parameter in source.partitions(h):
            for ein in range(max_weight + 1):
                for input_ in source.partitions(ein):
                    for eout in range(max_weight + 3):
                        old = mode(schur(parameter), schur(input_), eout)
                        for output in source.partitions(eout):
                            expected = coefficient(old, output)
                            actual = new.schur_mode(parameter, output, input_)
                            assert actual == expected
                            integer_modes += 1
                            maximum_integer_mode_bits = max(maximum_integer_mode_bits,
                                                            abs(actual).bit_length())
                            for prime in (2, 3, 5):
                                assert new.schur_mode(parameter, output, input_, prime) == actual % prime
                                field_modes += 1

    crossing_records = []
    for h in range(1, max_crossing + 1):
        basis = source.partitions(h)
        d = len(basis)
        pieces = [source_schur_piece(h, s) for s in range(-2*h, 2*h+1)]
        fields = new.binary_interpolation_data(h)[0]
        support_union = set()
        for x, px in enumerate(basis):
            for y, py in enumerate(basis):
                for u, pu in enumerate(basis):
                    for v, pv in enumerate(basis):
                        expected = tuple(piece[x*d+y][v*d+u] for piece in pieces)
                        actual = new.crossing_coefficients_binary(px, py, pu, pv)
                        assert actual == tuple(value % 2 for value in expected)
                        assert sum(actual) % 2 == (x == v and y == u)
                        support_union.update(k-2*h for k, value in enumerate(actual) if value)
                        binary_coefficients += len(actual)
                        if h <= 3:
                            for z in (2, 3):
                                value = new.crossing_value(px, py, pu, pv, z)
                                assert value == sum(c * z**k for k, c in enumerate(expected))
                                rational_values += 1
        crossing_records.append({'h': h, 'dimension': d, 'pieces': 4*h+1,
                                 'extension_degree': fields.degree,
                                 'irreducible_polynomial_bits': fields.modulus,
                                 'nonzero_shift_support': sorted(support_union)})

    large_queries = []
    for x, y in (((4,4,4,4), (1,15)),
                 ((4,4,4,4,4,4,4,4), (1,1,2,3,5,8,12))):
        h = sum(x)
        assert h == sum(y)
        coefficients = new.crossing_coefficients_binary(x, y, y, x)
        assert sum(coefficients) % 2 == 1
        field = new.binary_interpolation_data(h)[0]
        # Evaluate away from the interpolation nodes as an independent check.
        extra = field.order - 1
        assert extra not in new.binary_interpolation_data(h)[1]
        value = 0
        for bit in reversed(coefficients):
            value = field.multiply(value, extra) ^ bit
        assert value == new.crossing_value_binary(x, y, y, x, extra, field)
        large_queries.append({'h': h, 'x': x, 'y': y,
                              'nonzero_shifts': [k-2*h for k, c in enumerate(coefficients) if c],
                              'sum': 1, 'extra_evaluation_point': extra})

    assert new.schur_mode((1,1), (), (1,1)) == -1
    correct = new.crossing_coefficients_binary((2,), (2,), (2,), (2,))
    wrong_conjugation = new.crossing_coefficients_binary((1,1), (2,), (2,), (2,))
    assert sum(correct) % 2 == 1 and sum(wrong_conjugation) % 2 == 0
    return {'status':'pass',
            'scope':'Exact finite identities and polynomial query construction; no stronger factoring exponent.',
            'source_revision':'16a5ce41304fcdaea4e0be7c132e1c361dc26da5',
            'integral_schur_predecessor':'28105165edb8e9a0427468dc833f3c09bc0a7fb7',
            'integer_mode_entries':integer_modes,
            'prime_field_mode_entries':field_modes,
            'maximum_integer_mode_bits':maximum_integer_mode_bits,
            'exact_rational_crossing_values':rational_values,
            'binary_crossing_coefficients':binary_coefficients,
            'crossings':crossing_records,'larger_entry_queries':large_queries,
            'negative_controls':{'unsigned_homogeneous_self_pairing':'actual integer value is -1',
                                 'wrong_label_conjugation':'endpoint changes from one to zero'},
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'query_sha256':sha256(Path(new.__file__).read_bytes()).hexdigest(),
            'boson_reference_sha256':sha256(Path(source.__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-weight',type=int,default=4)
    parser.add_argument('--max-crossing',type=int,default=4)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    report=verify(args.max_weight,args.max_crossing)
    rendered=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
