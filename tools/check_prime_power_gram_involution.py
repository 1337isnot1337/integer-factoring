#!/usr/bin/env python3
"""Exact Gram involutions and streamed symmetric-function coefficient oracles.

The character and mode-entry routines below retain only polynomial-size
integer/rational working data. The verification harness also materializes
small comparison matrices; those matrices are not part of the oracle.
No floating point or inversion of a characteristic-divisible denominator
is used.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import argparse
import json

from check_partition_channels import operator, partition_denominator
from check_schur_fusion_gram import (
    character, divided_shift, homogeneous_pairing_polynomial,
)
from check_integral_schur_crossing import gram_bits


def stream_partitions(weight, least=1):
    """No cache and no retained list of all partitions."""
    if weight == 0:
        yield ()
    for first in range(least, weight + 1):
        for tail in stream_partitions(weight - first, first):
            yield (first,) + tail


def determinant_bareiss(matrix):
    """Fraction-free integer determinant, with exact divisions."""
    size = len(matrix)
    if not size:
        return 1
    values = [row[:] for row in matrix]
    sign, previous = 1, 1
    for column in range(size - 1):
        pivot_row = next((r for r in range(column, size) if values[r][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            values[pivot_row], values[column] = values[column], values[pivot_row]
            sign = -sign
        pivot = values[column][column]
        for row in range(column + 1, size):
            for j in range(column + 1, size):
                numerator = pivot * values[row][j] - values[row][column] * values[column][j]
                assert numerator % previous == 0
                values[row][j] = numerator // previous
            values[row][column] = 0
        previous = pivot
    return sign * values[-1][-1]


def character_finite_difference(shape, cycle_type):
    """One irreducible character in exp(O(sqrt(n))) time/poly(n) space."""
    weight = sum(shape)
    assert weight == sum(cycle_type)
    multiplicities = sorted(Counter(cycle_type).items())
    denominator = factorial(weight)
    descending = shape[::-1]
    length = len(shape)
    accumulated = 0
    ranges = [range(multiplicity + 1) for _, multiplicity in multiplicities]
    for point in product(*ranges):
        power_sums = dict(zip((part for part, _ in multiplicities), point))
        complete = [Fraction(1)]
        for degree in range(1, weight + 1):
            complete.append(sum((Fraction(value) * complete[degree - part]
                                 for part, value in power_sums.items() if part <= degree),
                                Fraction(0)) / degree)
        entries = []
        for row, part in enumerate(descending):
            values = []
            for column in range(length):
                degree = part - row + column
                value = denominator * complete[degree] if degree >= 0 else Fraction(0)
                assert value.denominator == 1
                values.append(value.numerator)
            entries.append(values)
        sign_weight = prod((-1) ** (multiplicity - value) * comb(multiplicity, value)
                           for (_, multiplicity), value in zip(multiplicities, point))
        accumulated += sign_weight * determinant_bareiss(entries)
    numerator = accumulated * prod(part ** multiplicity for part, multiplicity in multiplicities)
    common = denominator ** length
    assert numerator % common == 0
    return numerator // common


def power_sum_mode(parameter, output, input_):
    """One literal normal-order coefficient, with no recursion-state cache."""
    input_counts = Counter(input_)
    output_counts = Counter(output)
    shared = sorted(set(input_counts) & set(output_counts))
    tail_ranges = [range(min(input_counts[k], output_counts[k]) + 1) for k in shared]
    parameter_counts = Counter(parameter)
    answer = 0
    for tail_values in product(*tail_ranges):
        tail = dict(zip(shared, tail_values))
        created = [(k, count - tail.get(k, 0)) for k, count in sorted(output_counts.items())
                   if count > tail.get(k, 0)]
        destroyed = [(k, count - tail.get(k, 0)) for k, count in sorted(input_counts.items())
                     if count > tail.get(k, 0)]
        variables = [(False, k, count) for k, count in created]
        variables += [(True, k, count) for k, count in destroyed]
        if sum(count for _, _, count in variables) != len(parameter):
            continue
        accumulated = 0
        for point in product(*(range(count + 1) for _, _, count in variables)):
            value = 1
            for mode, multiplicity in parameter_counts.items():
                linear = 0
                for (annihilate, target, _), coordinate in zip(variables, point):
                    if annihilate:
                        coefficient = (-1) ** (mode - 1) * target * comb(target + mode - 1, mode - 1)
                    else:
                        coefficient = comb(target - 1, mode - 1) if target >= mode else 0
                    linear += coefficient * coordinate
                value *= linear ** multiplicity
            sign_weight = prod((-1) ** (count - coordinate) * comb(count, coordinate)
                               for (_, _, count), coordinate in zip(variables, point))
            accumulated += sign_weight * value
        denominator = prod(factorial(count) for _, _, count in variables)
        assert accumulated % denominator == 0
        coefficient = accumulated // denominator
        coefficient *= prod(factorial(count) // factorial(tail.get(k, 0))
                            for k, count in input_counts.items())
        answer += coefficient
    return answer


def matrix_multiply(left, right):
    columns = list(zip(*right))
    return [[sum(a * b for a, b in zip(row, column)) for column in columns] for row in left]


def schur_gram(weight):
    """Small verification matrix; each scalar entry uses the new oracles."""
    basis = tuple(stream_partitions(weight))
    denominator = factorial(weight)
    characters = [[character_finite_difference(shape, mu) for mu in basis] for shape in basis]
    scaled = [[denominator * value // partition_denominator(mu)
               for mu, value in zip(basis, row)] for row in characters]
    source = [[power_sum_mode(a, (), b) for b in basis] for a in basis]
    raw = matrix_multiply(matrix_multiply(scaled, source), list(zip(*scaled)))
    scale = denominator ** 2
    assert all(value % scale == 0 for row in raw for value in row)
    return [[value // scale for value in row] for row in raw]


def check_integer_square():
    operator_cases = coefficient_cases = primitive_cases = 0
    for weight in range(1, 10):
        for mu in stream_partitions(weight):
            left, right = defaultdict(int), defaultdict(int)
            for nu, coefficient in homogeneous_pairing_polynomial(mu):
                for rho, other in homogeneous_pairing_polynomial(nu):
                    left[rho] += coefficient * other
            for k in range(weight):
                shifted = divided_shift(divided_shift({mu: 1}, -1, k), 1, k)
                for rho, coefficient in shifted.items():
                    right[rho] += comb(2 * weight + k - 1, k) * coefficient
            assert {mu: value for mu, value in left.items() if value} == {
                mu: value for mu, value in right.items() if value}
            operator_cases += 1
    for weight in range(1, 41):
        for k in range(weight):
            collected = sum(comb(k, a) * comb(k, b) * comb(2 * weight - k - 1, a + b - k)
                            for a in range(k + 1) for b in range(k + 1) if a + b >= k)
            assert collected == comb(2 * weight + k - 1, k)
            coefficient_cases += 1
        for t in range(weight):
            assert comb(2 * weight - 1, t) ** 2 == sum(
                comb(2 * weight + k - 1, k) * comb(t, k) * comb(2 * weight - t - 1, k)
                for k in range(t + 1))
            primitive_cases += 1
    return {"integer_operator_square_cases": operator_cases,
            "divided_power_collection_cases": coefficient_cases,
            "primitive_chain_scalar_cases": primitive_cases}


def verify():
    character_cases = source_cases = 0
    maximum_character_grid = 0
    for weight in range(9):
        for shape in stream_partitions(weight):
            for cycles in stream_partitions(weight):
                actual = character_finite_difference(shape, cycles)
                assert actual == character(shape, cycles)
                character_cases += 1
                maximum_character_grid = max(maximum_character_grid,
                                             prod(value + 1 for value in Counter(cycles).values()))
    parameters = [mu for weight in range(7) for mu in stream_partitions(weight)]
    states = [mu for weight in range(6) for mu in stream_partitions(weight)]
    for mu in parameters:
        for output in states:
            for input_ in states:
                assert power_sum_mode(mu, output, input_) == operator(mu, output, input_)
                source_cases += 1

    matrices = {}
    records = []
    for prime, weight in ((2, 1), (2, 2), (2, 4), (2, 8),
                          (3, 1), (3, 3), (3, 9), (5, 1), (5, 5)):
        if weight not in matrices:
            matrices[weight] = schur_gram(weight)
        gram = matrices[weight]
        reduced = [[value % prime for value in row] for row in gram]
        square = matrix_multiply(reduced, reduced)
        dimension = len(gram)
        assert all(value % prime == int(i == j) for i, row in enumerate(square)
                   for j, value in enumerate(row))
        assert all(reduced[i][j] == reduced[j][i] for i in range(dimension) for j in range(dimension))
        if prime == 2:
            assert [sum(value << j for j, value in enumerate(row)) for row in reduced] == gram_bits(weight)
        records.append({"prime": prime, "weight": weight, "dimension": dimension,
                        "entrywise_identity": reduced == [[int(i == j) for j in range(dimension)]
                                                          for i in range(dimension)],
                        "square_identity_entries": dimension ** 2,
                        "maximum_integral_gram_entry_bits": max(abs(value).bit_length()
                                                               for row in gram for value in row)})

    # A power of the coefficient prime matters. Degree three is not a
    # power of two, and its full-Schur Gram square is not the identity.
    non_prime_power = [[value % 2 for value in row] for row in matrices[3]]
    negative_square = matrix_multiply(non_prime_power, non_prime_power)
    assert any(value % 2 != int(i == j) for i, row in enumerate(negative_square)
               for j, value in enumerate(row))
    assert not next(record for record in records if record["prime"] == 2
                    and record["weight"] == 8)["entrywise_identity"]
    assert power_sum_mode((1,), (1,), (1, 1)) == 2
    assert power_sum_mode((2,), (), (2,)) == -6

    result = {"scope": "integral Gram involution and streamed coefficient algorithms; no improved joint width exponent",
              "independent_character_entries": character_cases,
              "largest_character_grid_in_small_suite": maximum_character_grid,
              "independent_source_mode_entries": source_cases,
              "prime_power_gram_cases": records,
              "negative_controls": {"non_prime_power_square_fails": True,
                                    "involution_need_not_be_identity": True,
                                    "annihilation_multiplicities_and_signs_retained": True}}
    result.update(check_integer_square())
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = verify()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(rendered, end="")
