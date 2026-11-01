#!/usr/bin/env python3
"""Exact checks for the SoS crossing translation and transition profile."""

from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path

from check_partition_channels import rank_mod


def digits(value: int, length: int) -> tuple[int, ...]:
    result = [0] * length
    for position in range(length - 1, -1, -1):
        value, result[position] = divmod(value, 4)
    assert value == 0
    return tuple(result)


def encode(values: tuple[int, ...] | list[int]) -> int:
    result = 0
    for value in values:
        result = 4 * result + value
    return result


def partial_transpose(matrix: list[list[int]], length: int, mask: int) -> list[list[int]]:
    size = 4 ** length
    words = [digits(value, length) for value in range(size)]
    result = [[0] * size for _ in range(size)]
    for row, word in enumerate(words):
        for column, other in enumerate(words):
            left, right = list(word), list(other)
            for position in range(length):
                if mask >> position & 1:
                    left[position], right[position] = right[position], left[position]
            result[row][column] = matrix[encode(left)][encode(right)]
    return result


def determinant_gram(length: int, prime: int) -> list[list[int]]:
    """Literal coefficients of the determinant of a product of matrices."""
    size = 4 ** length
    result = [[0] * size for _ in range(size)]
    half = pow(2, -1, prime)
    for first_path in range(1 << (length + 1)):
        first = [first_path >> position & 1 for position in range(length + 1)]
        for inner_path in range(1 << (length - 1)):
            second = ([1 - first[0]]
                      + [inner_path >> position & 1 for position in range(length - 1)]
                      + [1 - first[-1]])
            row = encode([2 * first[j] + first[j + 1] for j in range(length)])
            column = encode([2 * second[j] + second[j + 1] for j in range(length)])
            result[row][column] = half * (-1) ** (first[0] + first[-1]) % prime
    return result


def regroup(word: tuple[int, ...]) -> tuple[int, ...]:
    first = [value // 2 for value in word]
    second = [value % 2 for value in word]
    return (2 * first[0] + second[-1],) + tuple(
        2 * second[j] + first[j + 1] for j in range(len(word) - 1))


def coefficient_rank(length: int, degree: int, prime: int) -> int:
    size = length - 1
    return 4 * sum(comb(size, k) * 3 ** k for k in range(degree + 1)
                   if comb(size - k, degree - k) % prime)


def gamma_two(matrix: list[list[int]], dimension: int) -> list[list[int]]:
    return [[matrix[x * dimension + u][v * dimension + y]
             for v in range(dimension) for u in range(dimension)]
            for x in range(dimension) for y in range(dimension)]


def verify_pairs() -> int:
    first, second = (1, 0, 0, -1), (0, 1, 1, 0)
    gram = [[first[i] * first[j] + second[i] * second[j] for j in range(4)]
            for i in range(4)]
    perturbation = [[0] * 4 for _ in range(4)]
    perturbation[0][3] = 1
    transposed = gamma_two(perturbation, 2)
    skew = [[perturbation[i][j] - transposed[i][j] for j in range(4)]
            for i in range(4)]
    cases = 0
    for prime in (3, 5, 7):
        for order in (0, 1, 5):
            for parameter in range(1, prime):
                matrix = [[(pow(parameter, order, prime) * gram[i][j]
                            + pow(parameter, order + 1, prime) * skew[i][j]) % prime
                           for j in range(4)] for i in range(4)]
                other = gamma_two(matrix, 2)
                assert all((matrix[i][j] + other[i][j]) % prime
                           == 2 * pow(parameter, order, prime) * (i == j) % prime
                           for i in range(4) for j in range(4))
                row, column = rank_mod(matrix, prime), rank_mod(other, prime)
                assert min(row, column) >= 2 and row + column >= 4
                assert row * column >= 2 * (4 - 2)
                cases += 1
    # In characteristic two the displayed identity would give no rank
    # information: the zero matrix satisfies A+Gamma(A)=2I.
    zero = [[0] * 4 for _ in range(4)]
    assert all((zero[i][j] + gamma_two(zero, 2)[i][j]) % 2 == 2 * (i == j) % 2
               for i in range(4) for j in range(4))
    assert rank_mod(zero, 2) ** 2 < 2 * (4 - 2)
    return cases


def verify_matrices() -> tuple[list[dict], dict[str, int]]:
    reports = []
    counts = {"orbit_factorization_entries": 0,
              "crossing_sum_entries": 0,
              "transpose_entries": 0,
              "literal_rank_checks": 0}
    for prime in (3, 5, 7):
        for length in range(1, 5):
            size = 4 ** length
            number = length - 1
            base = determinant_gram(length, prime)
            orbit = [partial_transpose(base, length, mask)
                     for mask in range(1 << number)]
            transitions = [(mask ^ (mask >> 1)).bit_count()
                           for mask in range(1 << number)]
            complete_mask = transitions.index(number)
            complete = orbit[complete_mask]
            permutation, inverses = [], []
            for row in complete:
                nonzero = [j for j, value in enumerate(row) if value]
                assert len(nonzero) == 1
                column = nonzero[0]
                permutation.append(column)
                inverses.append(pow(row[column], -1, prime))
            assert len(set(permutation)) == size
            coordinates = [regroup(digits(index, length)) for index in range(size)]

            for mask, matrix in enumerate(orbit):
                transition_bits = mask ^ (mask >> 1)
                for row, left in enumerate(coordinates):
                    for column, right in enumerate(coordinates):
                        expected = int(left[0] == right[0])
                        for edge in range(number):
                            a, b = left[edge + 1], right[edge + 1]
                            if transition_bits >> edge & 1:
                                expected *= int(a == b)
                            else:
                                expected *= -int(a in (0, 3) and b in (0, 3))
                        actual = matrix[row][permutation[column]] * inverses[column] % prime
                        assert actual == expected % prime
                        counts["orbit_factorization_entries"] += 1

            pieces = [[[sum(orbit[mask][row][column] for mask in range(1 << number)
                            if transitions[mask] == degree) % prime
                        for column in range(size)] for row in range(size)]
                      for degree in range(length)]
            half = pow(2, -1, prime)
            for row in range(size):
                for column in range(size):
                    # 2^(d-1) Q^(tensor d), directly in the determinant basis.
                    expected = (half * (-1) ** row.bit_count() % prime
                                if column == size - 1 - row else 0)
                    assert sum(piece[row][column] for piece in pieces) % prime == expected
                    counts["crossing_sum_entries"] += 1

            alternating = sum(1 << position for position in range(0, length, 2))
            ranks = []
            for degree, piece in enumerate(pieces):
                other = partial_transpose(piece, length, alternating)
                assert other == pieces[number - degree]
                counts["transpose_entries"] += size * size
                rank = rank_mod(piece, prime)
                assert rank == coefficient_rank(length, degree, prime)
                ranks.append(rank)
                counts["literal_rank_checks"] += 1
            reports.append({"prime": prime, "physical_positions": length,
                            "row_ranks": ranks,
                            "maximum_rank_product": max(a * b for a, b in zip(ranks, reversed(ranks)))})
    return reports, counts


def verify_integer_envelopes() -> int:
    cases = 0
    for prime in (3, 5, 7):
        for length in range(1, 65):
            number = length - 1
            ranks = [coefficient_rank(length, degree, prime) for degree in range(length)]
            for degree in range(length):
                product = ranks[degree] * ranks[number - degree]
                assert product >= 16 * comb(number, degree) ** 2 * 3 ** number
                assert product <= 16 * (number + 1) * 12 ** number
                cases += 1
            middle = number // 2
            assert (ranks[middle] * ranks[number - middle] * (number + 1) ** 2
                    >= 16 * 12 ** number)
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("evidence/sos-partial-transpose.json"))
    args = parser.parse_args()
    pair_cases = verify_pairs()
    reports, counts = verify_matrices()
    result = {"status": "passed",
              "scope": "Exact paired-degeneration cases and finite determinant-orbit profiles; the all-size statements have separate proofs.",
              "paired_identity_and_rank_cases": pair_cases,
              **counts,
              "integer_rank_envelope_cases": verify_integer_envelopes(),
              "negative_control": "Characteristic-two zero matrix satisfies the paired identity but violates the rank conclusion.",
              "finite_profiles": reports}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "finite_profiles"}, indent=2))


if __name__ == "__main__":
    main()
