#!/usr/bin/env python3
"""Exact checks for the bosonic moment and cumulative-profile identities.

The all-weight proofs are in research/bosonic-hermite-profile.md.  These
finite checks use the independent integral source-operator recursion from
check_partition_channels.py.  They do not establish a new factoring rate.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import argparse
import hashlib
import json
from math import comb, factorial
from pathlib import Path

from check_partition_channels import (
    column_matrix,
    dual_polynomial,
    operator,
    partition_denominator,
    partitions,
    piece_matrix,
    rank_mod,
)


Polynomial = dict[tuple[int, ...], Fraction]


def add_scaled(target: Polynomial, source: Polynomial, scale) -> None:
    for monomial, value in source.items():
        target[monomial] = target.get(monomial, 0) + scale * value
        if not target[monomial]:
            del target[monomial]


def lower(poly: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for monomial, value in poly.items():
        for part, multiplicity in Counter(monomial).items():
            if part > 1:
                changed = list(monomial)
                changed.remove(part)
                changed.append(part - 1)
                add_scaled(result, {tuple(sorted(changed)): value}, part * multiplicity)
    return result


def virasoro_minus(index: int, poly: Polynomial) -> Polynomial:
    """L_(-index), with L_0 equal to the energy operator."""
    if index < 0:
        raise ValueError("nonnegative index required")
    result: Polynomial = {}
    for monomial, value in poly.items():
        for part, multiplicity in Counter(monomial).items():
            changed = list(monomial)
            changed.remove(part)
            changed.append(part + index)
            add_scaled(result, {tuple(sorted(changed)): value}, part * multiplicity)
        for part in range(1, index):
            changed = tuple(sorted(monomial + (part, index - part)))
            add_scaled(result, {changed: value}, Fraction(1, 2))
    return result


@lru_cache(None)
def k_basis(monomial: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], Fraction], ...]:
    derivative: Polynomial = {monomial: Fraction(1)}
    result: Polynomial = {}
    index = 0
    while derivative:
        add_scaled(result, virasoro_minus(index, derivative), Fraction((-1) ** index, factorial(index)))
        derivative = lower(derivative)
        index += 1
    return tuple(sorted(result.items()))


def k_operator(poly: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for monomial, value in poly.items():
        add_scaled(result, dict(k_basis(monomial)), value)
    return result


@lru_cache(None)
def gamma_two(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Sum of full contractions after leaving two t_1 factors unmatched."""
    na, nb = a.count(1), b.count(1)
    result = 0
    if na >= 2:
        result += comb(na, 2) * operator(a[2:], (), b)
    if na and nb:
        result += na * nb * operator(a[1:], (), b[1:])
    if nb >= 2:
        result += comb(nb, 2) * operator(a, (), b[2:])
    return result


@lru_cache(None)
def source_gamma_two(x: tuple[int, ...], v: tuple[int, ...]) -> int:
    return sum(coefficient * gamma_two(term, v)
               for term, coefficient in dual_polynomial(x))


def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(None)
def cumulative_angle(x: tuple[int, ...], shift: int, energy: int,
                     output: tuple[int, ...], input_: tuple[int, ...],
                     row_cut: bool) -> int:
    h = sum(x)
    result = 0
    for term, coefficient in dual_polynomial(x):
        length = h + sum(term)
        if row_cut:
            k = energy - shift
            weight = (-1) ** k * (
                energy * choose(length - 1, k - 1)
                - length * choose(length - 2, k - 2)) if k >= 0 else 0
        else:
            k = 2 * h - energy - shift - 1
            weight = -(-1) ** k * (
                (2 * h - energy) * choose(length - 1, k)
                - length * choose(length - 2, k - 1)) if k >= 0 else 0
        result += coefficient * weight * operator(term, output, input_)
    return result


def cumulative_from_channels(h: int, shift: int, row_cut: bool) -> list[list[int]]:
    basis = partitions(h)
    d = len(basis)
    size = d * d
    result = [[0] * size for _ in range(size)]
    top = 2 * h + shift if row_cut else 2 * h - shift - 1
    for energy in range(max(-1, top) + 1):
        for middle in partitions(energy):
            for xi, x in enumerate(basis):
                for vi, v in enumerate(basis):
                    scale = factorial(h) // partition_denominator(v)
                    for yi, y in enumerate(basis):
                        for ui, u in enumerate(basis):
                            if row_cut:
                                value = cumulative_angle(x, shift, energy, middle, y, True)
                                value *= scale * operator(v, u, middle)
                            else:
                                value = cumulative_angle(x, shift, energy, u, middle, False)
                                value *= scale * operator(v, middle, y)
                            result[xi * d + yi][vi * d + ui] += value
    return result


def rank_fraction(matrix) -> int:
    values = [[Fraction(value) for value in row] for row in matrix]
    if not values:
        return 0
    rank = 0
    for column in range(len(values[0])):
        pivot = next((row for row in range(rank, len(values)) if values[row][column]), None)
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        divisor = values[rank][column]
        values[rank] = [value / divisor for value in values[rank]]
        for row in range(rank + 1, len(values)):
            multiplier = values[row][column]
            if multiplier:
                values[row] = [value - multiplier * other
                               for value, other in zip(values[row], values[rank])]
        rank += 1
        if rank == len(values):
            break
    return rank


def inverse(matrix) -> list[list[Fraction]]:
    size = len(matrix)
    values = [[Fraction(value) for value in row]
              + [Fraction(i == j) for j in range(size)]
              for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if values[row][column])
        values[column], values[pivot] = values[pivot], values[column]
        divisor = values[column][column]
        values[column] = [value / divisor for value in values[column]]
        for row in range(size):
            if row != column and values[row][column]:
                multiplier = values[row][column]
                values[row] = [value - multiplier * other
                               for value, other in zip(values[row], values[column])]
    return [row[size:] for row in values]


def pivot_minor(matrix) -> tuple[list[int], list[int]]:
    values = [[Fraction(value) for value in row] for row in matrix]
    labels = list(range(len(values)))
    columns = []
    rank = 0
    for column in range(len(values[0])):
        pivot = next((row for row in range(rank, len(values)) if values[row][column]), None)
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        labels[rank], labels[pivot] = labels[pivot], labels[rank]
        divisor = values[rank][column]
        values[rank] = [value / divisor for value in values[rank]]
        for row in range(rank + 1, len(values)):
            multiplier = values[row][column]
            if multiplier:
                values[row] = [value - multiplier * other
                               for value, other in zip(values[row], values[rank])]
        columns.append(column)
        rank += 1
    return labels[:rank], columns


def active_profile(h: int, profile: list[list[list[int]]], gamma: list[list[int]]) -> dict:
    """Normalize an actual nonsingular moment block, then measure both cuts."""
    rows, columns = pivot_minor(gamma)
    rank = len(rows)
    d = len(partitions(h))
    block_inverse = inverse([[gamma[i][j] for j in columns] for i in rows])
    pieces = []
    for matrix in profile:
        piece = [[sum(block_inverse[x][k] * matrix[rows[k] * d + y][columns[v] * d + u]
                      for k in range(rank)) / (-4 * h)
                  for v in range(rank) for u in range(rank)]
                 for x in range(rank) for y in range(rank)]
        pieces.append(piece)
    size = rank * rank
    assert all(sum(matrix[i][j] for matrix in pieces) == int(i == j)
               for i in range(size) for j in range(size))
    records = []
    for shift, matrix in zip(range(-2 * h, 2 * h), pieces):
        records.append({"s": shift, "row_rank_Q": rank_fraction(matrix),
                        "column_rank_Q": rank_fraction(column_matrix(matrix, rank))})
    return {"h": h, "dimension": rank, "operator_rows": rows, "operator_columns": columns,
            "state_labels": list(range(rank)), "pieces": records,
            "normalization": "left inverse of exact Gamma_2 minor, then scalar -1/(4h)",
            "crossing_entries_checked": size * size}


def verify(max_profile_weight: int = 4, max_moment_weight: int = 12) -> dict:
    if not 1 <= max_profile_weight <= 6 or not max_profile_weight <= max_moment_weight <= 14:
        raise ValueError("profile weights 1..6 and moment weights through 14 are supported")
    moment_records = []
    pairing_entries = annihilation_checks = commutator_checks = 0
    prime = 1000003
    for h in range(1, max_moment_weight + 1):
        basis = partitions(h)
        size = len(basis)
        matrix = [[Fraction(0) for _ in basis] for _ in basis]
        positions = {mu: i for i, mu in enumerate(basis)}
        for column, mu in enumerate(basis):
            image = dict(k_basis(mu))
            assert lower(lower(image)) == {}
            annihilation_checks += 1
            for output, value in image.items():
                matrix[positions[output]][column] = value
            for index in range(min(h + 1, 5)):
                vector = {mu: Fraction(1)}
                commutator = lower(virasoro_minus(index, vector))
                add_scaled(commutator, virasoro_minus(index, lower(vector)), -1)
                expected = lower(vector) if index == 0 else virasoro_minus(index - 1, vector)
                assert commutator == {key: (index + 1) * value for key, value in expected.items()}
                commutator_checks += 1
        if h <= 9:
            for i, x in enumerate(basis):
                for j, v in enumerate(basis):
                    assert source_gamma_two(x, v) == partition_denominator(x) * matrix[i][j]
                    pairing_entries += 1
        modular = [[value.numerator * pow(value.denominator, -1, prime) % prime
                    for value in row] for row in matrix]
        expected_rank = len(partitions(h)) - len(partitions(h - 2))
        assert rank_mod(modular, prime) == expected_rank
        eigenspaces = {}
        for value in sorted({0, h, 2 - h}):
            eigenspaces[str(value)] = size - rank_mod(
                [[entry - value * (i == j) for j, entry in enumerate(row)]
                 for i, row in enumerate(modular)], prime)
        assert sum(eigenspaces.values()) == size
        if h >= 3:
            assert eigenspaces[str(h)] == len(partitions(h)) - len(partitions(h - 1))
            assert eigenspaces[str(2 - h)] == len(partitions(h - 1)) - len(partitions(h - 2))
            assert eigenspaces["0"] == len(partitions(h - 2))
        moment_records.append({"h": h, "dimension": size, "rank_mod_prime": expected_rank,
                               "eigenspace_dimensions_mod_prime": eigenspaces})

    profile_records = []
    active_records = []
    moment_entries = hermite_entries = channel_entries = 0
    for h in range(1, max_profile_weight + 1):
        basis = partitions(h)
        d = len(basis)
        size = d * d
        shifts = range(-2 * h, 2 * h + 1)
        matrices = [piece_matrix(h, shift) for shift in shifts]
        gamma = [[source_gamma_two(x, v) * factorial(h) // partition_denominator(v)
                  for v in basis] for x in basis]
        cumulative = []
        running = [[0] * size for _ in range(size)]
        for shift, matrix in zip(shifts, matrices):
            for i in range(size):
                for j in range(size):
                    running[i][j] += shift * matrix[i][j]
            if shift < 2 * h:
                cumulative.append([row[:] for row in running])
        assert not any(any(row) for row in running)
        assert cumulative[2 * h - 1] == cumulative[2 * h]
        for i in range(size):
            for j in range(size):
                values = [matrix[i][j] for matrix in matrices]
                constant = factorial(h) * int(i == j)
                assert sum(values) == constant
                assert sum(shift * value for shift, value in zip(shifts, values)) == 0
                second = sum(shift * (shift - 1) // 2 * value for shift, value in zip(shifts, values))
                x, y = divmod(i, d)
                v, u = divmod(j, d)
                assert second == 2 * h * gamma[x][v] * int(y == u)
                assert sum(matrix[i][j] for matrix in cumulative) == -2 * second
                moment_entries += 1
                residual = values[:]
                residual[2 * h] -= constant
                quotient = []
                for degree in range(4 * h - 1):
                    quotient.append(sum((degree - r + 1) * residual[r] for r in range(degree + 1)))
                for degree in range(4 * h + 1):
                    reconstructed = sum(multiplier * quotient[degree - offset]
                                        for offset, multiplier in ((0, 1), (1, -2), (2, 1))
                                        if 0 <= degree - offset < len(quotient))
                    assert reconstructed == residual[degree]
                    hermite_entries += 1
                assert quotient[2 * h - 1] == -cumulative[2 * h - 1][i][j]
                assert sum(quotient) == second
        # Every shift, both cuts; these are independently collapsed channel formulas.
        if h <= 3:
            for shift, matrix in zip(range(-2 * h, 2 * h), cumulative):
                for row_cut in (False, True):
                    assert cumulative_from_channels(h, shift, row_cut) == matrix
                    channel_entries += size * size
        central = cumulative[2 * h - 1]
        profile_records.append({
            "h": h, "dimension": d, "central_row_rank_Q": rank_fraction(central),
            "central_column_rank_Q": rank_fraction(column_matrix(central, d)),
            "central_declared_channel_dimension": sum(len(partitions(j)) for j in range(1, 2 * h)),
        })
        active_records.append(active_profile(h, cumulative, gamma))

    # Negative control: truncating the Virasoro telescoping sum at k=1 fails.
    bad = {(3,): Fraction(3)}
    add_scaled(bad, virasoro_minus(1, lower({(3,): Fraction(1)})), -1)
    assert lower(lower(bad))
    # A universal second-moment zero relation is false over Q.
    assert gamma_two((1,), (1,)) == 1
    # Omitting the unmatched-on-one-side terms changes the scalar moment.
    assert gamma_two((1, 1), (2,)) == 0
    assert gamma_two((1, 1), ()) == 1
    return {
        "scope": "exact finite verification supporting separate all-weight ordinary proofs",
        "prime_for_large_moment_matrices": prime,
        "moment_operator_records": moment_records,
        "exact_Hall_pairing_entries": pairing_entries,
        "exact_D_squared_K_checks": annihilation_checks,
        "exact_Virasoro_commutator_checks": commutator_checks,
        "exact_complete_profile_moment_entries": moment_entries,
        "exact_Hermite_reconstruction_entries": hermite_entries,
        "exact_two_cut_collapsed_channel_entries": channel_entries,
        "original_label_central_cumulatives": profile_records,
        "normalized_active_crossing_profiles": active_records,
        "negative_controls": ["truncated Virasoro telescoping sum", "false second-moment zero relation",
                              "omitted unmatched-on-one-side contractions"],
        "new_factoring_exponent_claimed": False,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-profile-weight", type=int, default=4)
    parser.add_argument("--max-moment-weight", type=int, default=12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.max_profile_weight, args.max_moment_weight)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")
