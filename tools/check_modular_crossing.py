#!/usr/bin/env python3
"""Exact checks of the normalized modular partition crossing.

Ordinary research checks, not a Lean proof or an asymptotic rank certificate.
The general family retains full polynomial operator labels and quotients
only input/output states. Source operators: anthropics/formal-math,
16a5ce41304fcdaea4e0be7c132e1c361dc26da5, Appendix B (Apache-2.0).
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial
import argparse
import json

from check_partition_channels import (
    column_matrix, dual_polynomial, operator, partition_denominator,
    partitions, rank_mod,
)


def regular(parts: tuple[int, ...], prime: int) -> bool:
    return all(k % prime for k in parts) and all(
        count < prime for count in Counter(parts).values())


@lru_cache(None)
def regular_partitions(prime: int, weight: int,
                       least: int = 1) -> tuple[tuple[int, ...], ...]:
    if weight == 0:
        return ((),)
    if weight < 0:
        return ()
    return tuple((part,) + tail for part in range(least, weight + 1)
                 if part % prime
                 for tail in regular_partitions(prime, weight - part, part)
                 if tail.count(part) < prime - 1)


def erase(parts: tuple[int, ...], part: int) -> tuple[int, ...]:
    index = parts.index(part)
    return parts[:index] + parts[index + 1:]


@lru_cache(None)
def modular_operator(prime: int, parts: tuple[int, ...],
                     output: tuple[int, ...], input_: tuple[int, ...]) -> int:
    """Source recursion; parts is a full polynomial monomial parameter."""
    if not parts:
        return int(output == input_)
    first, tail = parts[0], parts[1:]
    value = sum(comb(k - 1, first - 1)
                * modular_operator(prime, tail, erase(output, k), input_)
                for k in set(output) if k >= first)
    value += sum((-1) ** (first - 1) * comb(k + first - 1, first - 1)
                 * k * input_.count(k)
                 * modular_operator(prime, tail, output, erase(input_, k))
                 for k in set(input_))
    return value % prime


@lru_cache(None)
def modular_angle(prime: int, parts: tuple[int, ...], index: int,
                  output: tuple[int, ...], input_: tuple[int, ...]) -> int:
    weight = sum(parts)
    # Deliberately retain every operator-label monomial, including ones
    # excluded from the quotient's state basis.
    return sum(coefficient * comb(weight + sum(term), index)
               * modular_operator(prime, term, output, input_)
               for term, coefficient in dual_polynomial(parts)
               if weight + sum(term) >= index) % prime


def normalized_piece(prime: int, weight: int, shift: int,
                     labels: tuple[tuple[int, ...], ...],
                     exchanged: bool = True) -> list[list[int]]:
    dimension = len(labels)
    size = dimension * dimension
    result = [[0] * size for _ in range(size)]
    scales = [pow(partition_denominator(v) % prime, -1, prime) for v in labels]
    for index in range(2 * weight + 1):
        energy = shift + index if exchanged else 2 * weight - shift - index
        for intermediate in regular_partitions(prime, energy):
            if exchanged:
                left = [(-1) ** index * modular_angle(prime, x, index, intermediate, y)
                        % prime for x in labels for y in labels]
                right = [scale * modular_operator(prime, v, u, intermediate) % prime
                         for v, scale in zip(labels, scales) for u in labels]
                nonzero = [(j, value) for j, value in enumerate(right) if value]
                for row, value in zip(result, left):
                    if value:
                        for j, other in nonzero:
                            row[j] = (row[j] + value * other) % prime
            else:
                for xi, x in enumerate(labels):
                    for ui, u in enumerate(labels):
                        left = (-1) ** index * modular_angle(
                            prime, x, index, u, intermediate) % prime
                        if left:
                            for vi, (v, scale) in enumerate(zip(labels, scales)):
                                for yi, y in enumerate(labels):
                                    right = scale * modular_operator(prime, v, intermediate, y)
                                    row = result[xi * dimension + yi]
                                    col = vi * dimension + ui
                                    row[col] = (row[col] + left * right) % prime
    return result


def determinant_two(rows: list[int]) -> int:
    rows = rows[:]
    for column in range(len(rows)):
        pivot = next((i for i in range(column, len(rows))
                      if (rows[i] >> column) & 1), None)
        if pivot is None:
            return 0
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for i in range(column + 1, len(rows)):
            if (rows[i] >> column) & 1:
                rows[i] ^= rows[column]
    return 1


@lru_cache(None)
def determinant_operator_two(parts: tuple[int, ...], output: tuple[int, ...],
                             input_: tuple[int, ...]) -> int:
    """Independent permanent=determinant formula for regular p=2 labels."""
    assert regular(parts, 2) and regular(output, 2) and regular(input_, 2)
    remaining = len(output) + len(input_) - len(parts)
    if remaining < 0 or remaining % 2:
        return 0
    remaining //= 2
    common = tuple(sorted(set(output) & set(input_)))
    if remaining > len(common):
        return 0
    answer = 0
    for unchanged in combinations(common, remaining):
        unchanged = set(unchanged)
        columns = [(True, k) for k in output if k not in unchanged]
        columns += [(False, k) for k in input_ if k not in unchanged]
        matrix = []
        for m in parts:
            row = 0
            for j, (creation, k) in enumerate(columns):
                bit = ((m - 1) & ~(k - 1) == 0) if creation else ((m - 1) & k == 0)
                row |= int(bit) << j
            matrix.append(row)
        answer ^= determinant_two(matrix)
    return answer


def inverse_two(rows: list[int], dimension: int) -> list[int]:
    rows = [row | (1 << (dimension + i)) for i, row in enumerate(rows)]
    for column in range(dimension):
        pivot = next((i for i in range(column, dimension)
                      if (rows[i] >> column) & 1), None)
        assert pivot is not None, "homogeneous pairing is singular"
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for i in range(dimension):
            if i != column and (rows[i] >> column) & 1:
                rows[i] ^= rows[column]
    assert all(row & ((1 << dimension) - 1) == 1 << i for i, row in enumerate(rows))
    return [row >> dimension for row in rows]


def rank_two(rows: list[int]) -> int:
    pivots = {}
    for row in rows:
        while row:
            index = row.bit_length() - 1
            if index not in pivots:
                pivots[index] = row
                break
            row ^= pivots[index]
    return len(pivots)


def homogeneous_family_two(weight: int):
    labels = regular_partitions(2, weight)
    dimension = len(labels)
    gram = [sum(determinant_operator_two(x, (), y) << j
                for j, y in enumerate(labels)) for x in labels]
    inverse = inverse_two(gram, dimension)
    duals = [tuple(labels[j] for j in range(dimension) if (row >> j) & 1)
             for row in inverse]
    for i, terms in enumerate(duals):
        for j, y in enumerate(labels):
            assert sum(determinant_operator_two(x, (), y) for x in terms) % 2 == int(i == j)
    return labels, duals


def homogeneous_piece_two(weight: int, shift: int, labels, duals,
                          exchanged: bool = True) -> list[int]:
    dimension = len(labels)
    result = [0] * (dimension * dimension)

    def dual_operator(xi, output, input_):
        return sum(determinant_operator_two(x, output, input_) for x in duals[xi]) % 2

    for index in (0, 2 * weight):
        energy = shift + index if exchanged else 2 * weight - shift - index
        for intermediate in regular_partitions(2, energy):
            if exchanged:
                right = sum(determinant_operator_two(v, u, intermediate) << (vi * dimension + ui)
                            for vi, v in enumerate(labels) for ui, u in enumerate(labels))
                if right:
                    for xi in range(dimension):
                        for yi, y in enumerate(labels):
                            if dual_operator(xi, intermediate, y):
                                result[xi * dimension + yi] ^= right
            else:
                for xi in range(dimension):
                    for ui, u in enumerate(labels):
                        if dual_operator(xi, u, intermediate):
                            for vi, v in enumerate(labels):
                                for yi, y in enumerate(labels):
                                    if determinant_operator_two(v, intermediate, y):
                                        result[xi * dimension + yi] ^= 1 << (vi * dimension + ui)
    return result


def column_bits(rows: list[int], dimension: int) -> list[int]:
    return [sum(((rows[x * dimension + y] >> (v * dimension + u)) & 1)
                << (v * dimension + y) for v in range(dimension) for y in range(dimension))
            for x in range(dimension) for u in range(dimension)]


def verify() -> dict:
    counts = {}
    cancellation_entries = 0
    for h in range(1, 7):
        for mu in partitions(h):
            for nu in partitions(h):
                value = sum(c * operator(a, (), nu) for a, c in dual_polynomial(mu))
                assert value == partition_denominator(nu) * (mu == nu)
                cancellation_entries += 1
    counts["integral_cancellation_entries"] = cancellation_entries

    quotient_checks = 0
    count_checks = 0
    for prime in (2, 3, 5):
        for weight in range(17):
            actual = regular_partitions(prime, weight)
            assert set(actual) == {x for x in partitions(weight) if regular(x, prime)}
            count_checks += 1
        for weight in range(1, 8):
            for input_ in partitions(weight):
                if regular(input_, prime):
                    continue
                for energy in range(8):
                    for output in regular_partitions(prime, energy):
                        for m in range(1, 6):
                            assert modular_operator(prime, (m,), output, input_) == 0
                            quotient_checks += 1
    counts["regular_partition_counts"] = count_checks
    counts["state_ideal_operator_entries"] = quotient_checks

    label_checks = 0
    for prime in (2, 3, 5, 7):
        for width in range(1, 6):
            blocks = []
            for j in range(width):
                if prime == 2:
                    blocks.append(((8 * j + 1, 8 * j + 7), (8 * j + 3, 8 * j + 5)))
                else:
                    pairs = [(12 * j + a, 12 * j + 13 - a) for a in range(1, 7)
                             if (12 * j + a) % prime and (12 * j + 13 - a) % prime]
                    assert len(pairs) >= 2
                    blocks.append(tuple(pairs[:2]))
            labels = {tuple(sorted(k for pair in choices for k in pair))
                      for choices in product(*blocks)}
            expected = 8 * width ** 2 if prime == 2 else 12 * width ** 2 + width
            assert len(labels) == 2 ** width
            assert all(regular(label, prime) and sum(label) == expected for label in labels)
            label_checks += len(labels)
    counts["quadratic_energy_labels"] = label_checks

    general_records = []
    exchange_entries = crossing_entries = normalized_duality_entries = 0
    for prime, width in ((2, 1), (3, 1), (5, 2)):
        h = next(h for h in range(1, 12 * width ** 2 + width + 1)
                 if len(regular_partitions(prime, h)) >= 2 ** width)
        labels = tuple(sorted(regular_partitions(prime, h),
                              key=lambda x: tuple(reversed(x))))[:2 ** width]
        dimension = len(labels)
        size = dimension * dimension
        for xi, x in enumerate(labels):
            for vi, v in enumerate(labels):
                value = sum(c * modular_operator(prime, a, (), v)
                            for a, c in dual_polynomial(x))
                value = value * pow(partition_denominator(v) % prime, -1, prime) % prime
                assert value == int(xi == vi)
                normalized_duality_entries += 1
        total = [[0] * size for _ in range(size)]
        profile = []
        for shift in range(-2 * h, 2 * h + 1):
            matrix = normalized_piece(prime, h, shift, labels)
            assert matrix == normalized_piece(prime, h, shift, labels, exchanged=False)
            exchange_entries += size * size
            for i in range(size):
                for j in range(size):
                    total[i][j] = (total[i][j] + matrix[i][j]) % prime
            profile.append([shift, rank_mod(matrix, prime),
                            rank_mod(column_matrix(matrix, dimension), prime)])
        assert total == [[int(i == j) for j in range(size)] for i in range(size)]
        crossing_entries += size * size
        if prime == 2:
            assert factorial(h) % prime == 0 and total[0][0] == 1
        general_records.append({"prime": prime, "label_bits": width, "weight": h,
                                "row_column_rank_profile": profile})
    counts["normalized_duality_entries"] = normalized_duality_entries
    counts["general_exchange_entries"] = exchange_entries
    counts["general_crossing_entries"] = crossing_entries

    determinant_checks = 0
    regular_states = tuple(x for h in range(13) for x in regular_partitions(2, h))
    for h in range(1, 13):
        for a in regular_partitions(2, h):
            for output in regular_states:
                for input_ in regular_states:
                    assert determinant_operator_two(a, output, input_) == operator(a, output, input_) % 2
                    determinant_checks += 1
    counts["determinant_operator_entries"] = determinant_checks

    homogeneous_records = []
    homogeneous_exchange = homogeneous_crossing = 0
    for h in (8, 16, 32):
        labels, duals = homogeneous_family_two(h)
        dimension = len(labels)
        size = dimension * dimension
        assert [i for i in range(2 * h + 1) if comb(2 * h, i) % 2] == [0, 2 * h]
        central = homogeneous_piece_two(h, 0, labels, duals)
        row_rank = rank_two(central)
        col_rank = rank_two(column_bits(central, dimension))
        assert row_rank == col_rank == {8: 2, 16: 12, 32: 261}[h]
        record = {"prime": 2, "weight": h, "label_count": dimension,
                  "central_row_rank": row_rank, "central_column_rank": col_rank,
                  "top_channel_ambient_dimension": len(regular_partitions(2, 2 * h))}
        if h <= 16:
            total = [0] * size
            for shift in range(-2 * h, 2 * h + 1):
                matrix = homogeneous_piece_two(h, shift, labels, duals)
                assert matrix == homogeneous_piece_two(h, shift, labels, duals, exchanged=False)
                homogeneous_exchange += size * size
                total = [a ^ b for a, b in zip(total, matrix)]
            assert total == [1 << i for i in range(size)]
            homogeneous_crossing += size * size
            record["all_shifts_exchange_and_sum_checked"] = True
        homogeneous_records.append(record)
    counts["homogeneous_exchange_entries"] = homogeneous_exchange
    counts["homogeneous_crossing_entries"] = homogeneous_crossing

    odd_prime_grams = []
    for prime, h in ((3, 3), (3, 9), (5, 5)):
        labels = regular_partitions(prime, h)
        gram = [[modular_operator(prime, a, (), b) for b in labels] for a in labels]
        assert rank_mod(gram, prime) == len(labels)
        assert all(comb(2 * h - 1, h - j) % prime for j in range(1, h + 1))
        assert [i for i in range(2 * h + 1) if comb(2 * h, i) % prime] == [0, h, 2 * h]
        odd_prime_grams.append({"prime": prime, "weight": h, "dimension": len(labels),
                                "pairing_rank": len(labels)})

    return {"status": "pass", "scope": "finite exact modular identities and ranks; no asymptotic rank improvement",
            "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
            "checks": counts, "normalized_original_family": general_records,
            "homogeneous_characteristic_two": homogeneous_records,
            "odd_prime_power_pairings": odd_prime_grams}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="write the same JSON report to this file")
    args = parser.parse_args()
    report = verify()
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
