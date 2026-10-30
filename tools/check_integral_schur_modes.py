#!/usr/bin/env python3
"""Independent small integral-Schur checks for the pinned Appendix B modes.

Schur functions are built by Jacobi--Trudi, not a character recursion.
This finite audit does not prove an asymptotic factoring improvement.
The source recursion is imported from check_partition_channels.py.
"""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from math import comb, factorial, sqrt
from pathlib import Path
import argparse
import json

import check_partition_channels as source


def clean(poly):
    return tuple(sorted((p, c) for p, c in poly.items() if c))


def multiply(left, right):
    result = defaultdict(Fraction)
    for p, a in left:
        for q, b in right:
            result[tuple(sorted(p + q))] += a * b
    return clean(result)


@lru_cache(None)
def complete(weight):
    return tuple((p, Fraction(1, source.partition_denominator(p)))
                 for p in source.partitions(weight))


@lru_cache(None)
def schur(part):
    if not part:
        return (((), Fraction(1)),)
    lam = part[::-1]
    length = len(lam)
    result = defaultdict(Fraction)
    for sigma in permutations(range(length)):
        inversion = sum(sigma[i] > sigma[j]
                        for i in range(length) for j in range(i + 1, length))
        term = (((), Fraction((-1) ** inversion)),)
        for i, j in enumerate(sigma):
            term = multiply(term, complete(lam[i] - i + j))
            if not term:
                break
        for p, coefficient in term:
            result[p] += coefficient
    return clean(result)


def coefficient(poly, schur_label):
    coefficients = dict(poly)
    return sum(c * source.partition_denominator(p) * coefficients.get(p, 0)
               for p, c in schur(schur_label))


def transform(poly):
    result = defaultdict(Fraction)
    for p, c in poly:
        for q, d in source.dual_polynomial(p):
            result[q] += c * d
    return clean(result)


def mode(parameter, input_, output_weight):
    result = {}
    for out in source.partitions(output_weight):
        value = sum(a * b * source.operator(p, out, q)
                    for p, a in parameter for q, b in input_)
        if value:
            result[out] = value
    return clean(result)


def matrix_multiply(left, right):
    return [[sum(a * b for a, b in zip(row, column))
             for column in zip(*right)] for row in left]


def tensor_square(matrix):
    dimension = len(matrix)
    return [[matrix[a][x] * matrix[b][y]
             for x in range(dimension) for y in range(dimension)]
            for a in range(dimension) for b in range(dimension)]


def schur_piece(weight, shift, exchanged=True):
    basis = source.partitions(weight)
    change = [[dict(schur(label)).get(part, Fraction(0)) for part in basis]
              for label in basis]
    characters = [[entry * source.partition_denominator(part)
                   for entry, part in zip(row, basis)] for row in change]
    assert all(v.denominator == 1 for row in characters for v in row)
    old = source.piece_matrix(weight, shift, exchanged=exchanged)
    left = tensor_square(change)
    right = list(map(list, zip(*tensor_square(characters))))
    transformed = matrix_multiply(matrix_multiply(left, old), right)
    transformed = [[v / factorial(weight) for v in row] for row in transformed]
    assert all(v.denominator == 1 for row in transformed for v in row)
    return [[int(v) for v in row] for row in transformed]


def verify(max_weight=4, max_crossing=3):
    if not 1 <= max_weight <= 5:
        raise ValueError("The independent Jacobi--Trudi audit supports weights 1 through 5.")
    if not 1 <= max_crossing <= 4:
        raise ValueError("Use crossing weights 1 through 4 for this direct audit.")
    transform_checks = vacuum_checks = mode_checks = pairing_checks = 0
    max_mode_bit_length = 0
    for weight in range(max_weight + 1):
        basis = source.partitions(weight)
        for a in basis:
            for b in basis:
                assert coefficient(schur(a), b) == (a == b)
                pairing_checks += 1
        if not weight:
            continue
        expected = defaultdict(Fraction)
        for j in range(1, weight + 1):
            scalar = (-1) ** (weight - j) * comb(weight - 1, j - 1)
            for p, c in complete(j):
                expected[p] += scalar * (-1) ** (j - len(p)) * c
        assert transform(complete(weight)) == clean(expected)
        for label in basis:
            dual = transform(schur(label))
            assert transform(dual) == schur(label)
            for output_weight in range(weight + 1):
                for output in source.partitions(output_weight):
                    value = coefficient(dual, output)
                    assert value.denominator == 1
                    transform_checks += 1
            for input_label in basis:
                vacuum = dict(mode(dual, schur(input_label), 0)).get((), 0)
                assert vacuum == (label == input_label)
                vacuum_checks += 1
            for input_weight in range(max_weight + 1):
                for input_label in source.partitions(input_weight):
                    for output_weight in range(max_weight + 3):
                        result = mode(schur(label), schur(input_label), output_weight)
                        for output in source.partitions(output_weight):
                            value = coefficient(result, output)
                            assert value.denominator == 1
                            max_mode_bit_length = max(max_mode_bit_length,
                                                      abs(value.numerator).bit_length())
                            mode_checks += 1

    crossing_records = []
    exchange_entries = summed_entries = 0
    for weight in range(1, max_crossing + 1):
        dimension = len(source.partitions(weight))
        size = dimension * dimension
        summed = [[0] * size for _ in range(size)]
        ranks = []
        for shift in range(-2 * weight, 2 * weight + 1):
            piece = schur_piece(weight, shift)
            assert piece == schur_piece(weight, shift, exchanged=False)
            exchange_entries += size * size
            for a, row in enumerate(piece):
                for b, value in enumerate(row):
                    summed[a][b] += value
            other = source.column_matrix(piece, dimension)
            ranks.append({"shift": shift, "ranks_mod_2": [source.rank_mod(piece, 2),
                                                          source.rank_mod(other, 2)],
                          "ranks_mod_3": [source.rank_mod(piece, 3),
                                          source.rank_mod(other, 3)]})
        assert summed == [[int(a == b) for b in range(size)] for a in range(size)]
        summed_entries += size * size
        crossing_records.append({"weight": weight, "dimension": dimension,
                                 "pieces": ranks,
                                 "rank_sum_mod_2": sum(sqrt(r["ranks_mod_2"][0]
                                                             * r["ranks_mod_2"][1])
                                                       for r in ranks)})

    old_sum = [[sum(source.piece_matrix(2, s)[a][b] for s in range(-4, 5))
                for b in range(4)] for a in range(4)]
    assert all(value % 2 == 0 for row in old_sum for value in row)
    new_sum = [[sum(schur_piece(2, s)[a][b] for s in range(-4, 5))
                for b in range(4)] for a in range(4)]
    assert new_sum == [[int(a == b) for b in range(4)] for a in range(4)]
    central = schur_piece(2, 0)
    assert central == [[13, 4, 5, 18], [9, 6, 3, 12],
                       [12, 3, 6, 9], [18, 5, 4, 13]]
    # Forgetting the Hall denominator destroys Schur orthonormality.
    bad_pairing = sum(c * c for _, c in schur((2,)))
    assert bad_pairing == Fraction(1, 2)
    return {"status": "pass", "scope": "Finite exact checks; no improved asymptotic exponent.",
            "method": "Jacobi--Trudi over Q, integral output checked before reduction",
            "maximum_parameter_and_input_weight": max_weight,
            "maximum_output_weight": max_weight + 2,
            "schur_pairing_checks": pairing_checks,
            "integral_transform_coefficients": transform_checks,
            "unit_vacuum_pairings": vacuum_checks,
            "integral_mode_entries": mode_checks,
            "maximum_mode_coefficient_bit_length": max_mode_bit_length,
            "exchange_entries": exchange_entries, "summed_crossing_entries": summed_entries,
            "crossings": crossing_records, "weight_two_central_integer_matrix": central,
            "negative_controls": {"old_weight_two_endpoint_mod_2": "zero",
                                  "new_weight_two_endpoint_mod_2": "identity",
                                  "omitted_Hall_denominator_pairing": "1/2"},
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_checker_sha256": sha256(Path(source.__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-weight", type=int, default=4)
    parser.add_argument("--max-crossing", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = verify(args.max_weight, args.max_crossing)
    rendered = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
