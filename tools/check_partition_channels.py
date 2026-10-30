#!/usr/bin/env python3
"""Exact finite checks of the Appendix B partition pieces.

This implements the pinned source's integer definitions, not its Lean proof.
Full rank modulo a prime proves full rank over Q for the same integer matrix.
Small weights do not prove an asymptotic rank bound or a factoring rate.

Source: anthropics/formal-math, 16a5ce41304fcdaea4e0be7c132e1c361dc26da5,
eth/proof/rates/ClauseBound/Boson/{Defs,Defs2}.lean and
eth/proof/statements/AppHalfPiecesStatement.lean. The prior construction is
copyright 2026 Anthropic, PBC, Apache-2.0. This checker used OpenAI assistance.
"""

from collections import Counter, defaultdict
from functools import lru_cache
import argparse
import json
from math import comb, factorial


@lru_cache(None)
def partitions(weight: int, least: int = 1) -> tuple[tuple[int, ...], ...]:
    if weight == 0:
        return ((),)
    if weight < 0:
        return ()
    return tuple((part,) + tail for part in range(least, weight + 1)
                 for tail in partitions(weight - part, part))


def erase(parts: tuple[int, ...], part: int) -> tuple[int, ...]:
    index = parts.index(part)
    return parts[:index] + parts[index + 1:]


@lru_cache(None)
def operator(parts: tuple[int, ...], output: tuple[int, ...],
             input_: tuple[int, ...]) -> int:
    """Coefficient of output in (t^parts)[r] input, r=|output|-|input|.

    The recursion is (O2), equivalently Boson.nop. Creation removes one
    output part without a multiplicity. Annihilation includes its input
    multiplicity. The empty operator is the identity at r=0.
    """
    if not parts:
        return int(output == input_)
    first, tail = parts[0], parts[1:]
    value = sum(comb(part - 1, first - 1)
                * operator(tail, erase(output, part), input_)
                for part in set(output) if part >= first)
    value += sum((-1) ** (first - 1) * comb(part + first - 1, first - 1)
                 * part * input_.count(part)
                 * operator(tail, output, erase(input_, part))
                 for part in set(input_))
    return value


@lru_cache(None)
def dual_polynomial(parts: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], int], ...]:
    """a^mu = product (-1)^(mu_i-1) sum_j binom(mu_i,j) t_j."""
    if not parts:
        return (((), 1),)
    terms: dict[tuple[int, ...], int] = defaultdict(int)
    first = parts[0]
    for tail, coefficient in dual_polynomial(parts[1:]):
        for part in range(1, first + 1):
            terms[tuple(sorted((part,) + tail))] += (
                (-1) ** (first - 1) * comb(first, part) * coefficient)
    return tuple(sorted(terms.items()))


@lru_cache(None)
def angle_operator(parts: tuple[int, ...], index: int, output: tuple[int, ...],
                   input_: tuple[int, ...]) -> int:
    weight = sum(parts)
    return sum(coefficient * comb(sum(term) + weight, index)
               * operator(term, output, input_)
               for term, coefficient in dual_polynomial(parts)
               if sum(term) + weight >= index)


def partition_denominator(parts: tuple[int, ...]) -> int:
    value = 1
    for part, multiplicity in Counter(parts).items():
        value *= part ** multiplicity * factorial(multiplicity)
    return value


def piece_matrix(weight: int, shift: int, modulus: int | None = None,
                 exchanged: bool = True, alternating: bool = True) -> list[list[int]]:
    """Rows are (x,y), columns are (v,u), each index a partition of weight.

    exchanged=True implements Q'_s, using intermediate weight s+i.
    exchanged=False implements Q_s, using intermediate weight 2h-s-i.
    Comparing both forms checks the operator exchange identity independently.
    """
    basis = partitions(weight)
    dimension = len(basis)
    result = [[0] * (dimension * dimension) for _ in range(dimension * dimension)]
    scales = [factorial(weight) // partition_denominator(part) for part in basis]
    for index in range(2 * weight + 1):
        sign = (-1) ** index if alternating else 1
        intermediate_weight = shift + index if exchanged else 2 * weight - shift - index
        for intermediate in partitions(intermediate_weight):
            if exchanged:
                left = [sign * angle_operator(x, index, intermediate, y)
                        for x in basis for y in basis]
                right = [scale * operator(v, u, intermediate)
                         for v, scale in zip(basis, scales) for u in basis]
                if modulus is not None:
                    left = [value % modulus for value in left]
                    right = [value % modulus for value in right]
                nonzero_right = [(column, value) for column, value in enumerate(right) if value]
                for row, value in enumerate(left):
                    if value:
                        target = result[row]
                        for column, other in nonzero_right:
                            target[column] += value * other
            else:
                for x_index, x in enumerate(basis):
                    for u_index, u in enumerate(basis):
                        left = sign * angle_operator(x, index, u, intermediate)
                        if left:
                            for v_index, (v, scale) in enumerate(zip(basis, scales)):
                                for y_index, y in enumerate(basis):
                                    right = scale * operator(v, intermediate, y)
                                    result[x_index * dimension + y_index][
                                        v_index * dimension + u_index] += left * right
            if modulus is not None:
                result = [[value % modulus for value in row] for row in result]
    return result


def column_matrix(row_matrix: list[list[int]], dimension: int) -> list[list[int]]:
    return [[row_matrix[x * dimension + y][v * dimension + u]
             for v in range(dimension) for y in range(dimension)]
            for x in range(dimension) for u in range(dimension)]


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    values = [[value % prime for value in row] for row in matrix]
    rows, columns = len(values), len(values[0])
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if values[row][column]), None)
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        inverse = pow(values[rank][column], -1, prime)
        values[rank][column:] = [value * inverse % prime for value in values[rank][column:]]
        for row in range(rank + 1, rows):
            multiplier = values[row][column]
            if multiplier:
                values[row][column:] = [
                    (value - multiplier * pivot_value) % prime
                    for value, pivot_value in zip(values[row][column:], values[rank][column:])]
        rank += 1
        if rank == rows:
            break
    return rank


def verify(max_weight: int = 6) -> dict:
    if not 1 <= max_weight <= 8:
        raise ValueError("use weights 1 through 8; larger direct matrices need a separate budget")
    equality_entries = exchange_entries = 0
    for weight in range(1, min(max_weight, 3) + 1):
        dimension = len(partitions(weight))
        size = dimension * dimension
        summed = [[0] * size for _ in range(size)]
        for shift in range(-2 * weight, 2 * weight + 1):
            exchanged = piece_matrix(weight, shift)
            original = piece_matrix(weight, shift, exchanged=False)
            assert exchanged == original
            exchange_entries += size * size
            for row in range(size):
                for column in range(size):
                    summed[row][column] += original[row][column]
        for row in range(size):
            for column in range(size):
                assert summed[row][column] == factorial(weight) * (row == column)
                equality_entries += 1

    # A missing alternating sign invalidates exchange or the summed crossing.
    bad = sum(piece_matrix(1, shift, alternating=False)[0][0] for shift in range(-2, 3))
    assert bad != factorial(1)
    # A missing annihilation multiplicity also changes a literal operator entry.
    assert operator((1,), (1,), (1, 1)) == 2

    # The natural central matrix is not a positive Gram matrix.
    gram_counterexample = piece_matrix(2, 0)
    assert gram_counterexample == [[70, 0, 0, 34], [0, 2, -6, 0],
                                   [0, -10, -6, 0], [10, 0, 0, 10]]
    middle_minor = (gram_counterexample[1][1] * gram_counterexample[2][2]
                    - gram_counterexample[1][2] * gram_counterexample[2][1])
    assert middle_minor == -72
    assert gram_counterexample[0][3] != gram_counterexample[3][0]

    partition_counts = [1] + [0] * 30
    for part in range(1, 31):
        for weight in range(part, 31):
            partition_counts[weight] += partition_counts[weight - part]
    first_saving = next(weight for weight in range(1, 16)
                        if sum(partition_counts[:2 * weight + 1]) < partition_counts[weight] ** 2)
    assert first_saving == 15
    assert (sum(partition_counts[:31]), partition_counts[15] ** 2) == (28629, 30976)

    records = []
    binary_records = []
    for weight in range(1, max_weight + 1):
        dimension = len(partitions(weight))
        ranks = []
        for prime in (65521, 1000003):
            matrix = piece_matrix(weight, 0, prime)
            row_rank = rank_mod(matrix, prime)
            col_rank = rank_mod(column_matrix(matrix, dimension), prime)
            assert row_rank == col_rank == dimension * dimension
            ranks.append({"prime": prime, "row_rank": row_rank, "column_rank": col_rank})
            for width in range(1, dimension.bit_length()):
                if weight < 2 or not len(partitions(weight - 1)) < 2 ** width <= dimension:
                    continue
                # The source orders descending partition lists lexicographically.
                labels = sorted(range(dimension), key=lambda i: tuple(reversed(partitions(weight)[i])))[:2 ** width]
                indices = [x * dimension + y for x in labels for y in labels]
                restricted = [[matrix[row][column] for column in indices] for row in indices]
                binary_row_rank = rank_mod(restricted, prime)
                binary_col_rank = rank_mod(column_matrix(restricted, 2 ** width), prime)
                assert binary_row_rank == binary_col_rank == 2 ** (2 * width)
                binary_records.append({"label_bits": width, "weight": weight, "prime": prime,
                                       "row_rank": binary_row_rank, "column_rank": binary_col_rank})
        top_basis = partitions(2 * weight)
        split_count = 0
        for partition in top_basis:
            sums = {0}
            for part in partition:
                sums |= {value + part for value in sums if value + part <= weight}
            split_count += weight in sums
        top_map = [[operator(x, output, y) % 1000003 for output in top_basis]
                   for x in partitions(weight) for y in partitions(weight)]
        top_rank = rank_mod(top_map, 1000003)
        assert split_count <= top_rank <= len(top_basis)
        records.append({"weight": weight, "partition_count": dimension,
                        "matrix_dimension": dimension * dimension,
                        "declared_channels": sum(len(partitions(i)) for i in range(2 * weight + 1)),
                        "top_channel_dimension": len(top_basis),
                        "top_channel_rank_mod_1000003": top_rank,
                        "splittable_partitions": split_count,
                        "modular_ranks": ranks})
        operator.cache_clear()
        angle_operator.cache_clear()
    return {
        "source_commit": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "scope": "finite integer identities and exact finite-field ranks; no asymptotic rank claim",
        "crossing_identity_entries": equality_entries,
        "operator_exchange_entries": exchange_entries,
        "negative_controls": ["omitted alternating sign", "omitted annihilation multiplicity"],
        "natural_gram_counterexample": {"weight": 2, "matrix": gram_counterexample,
                                        "middle_principal_minor": middle_minor,
                                        "determinant": (70 * 10 - 34 * 10) * middle_minor},
        "first_declared_channel_saving": {"weight": first_saving, "channels": 28629,
                                          "matrix_dimension": 30976},
        "central_pieces": records,
        "binary_label_restrictions": binary_records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-weight", type=int, default=6)
    args = parser.parse_args()
    print(json.dumps(verify(args.max_weight), indent=2, sort_keys=True))
